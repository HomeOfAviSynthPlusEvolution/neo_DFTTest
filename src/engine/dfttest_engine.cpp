#include "engine/dfttest_engine.hpp"

#include <avisynth.h>

#include <dualsynth/global_lock.hpp>
#include <dualsynth/video_filter.hpp>

#include "dft_common.h"
#include "engine/dfttest_config.hpp"
#include "engine/pixel_plane.hpp"
#include "executor/dft_executor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace neo_dfttest {

namespace detail {

template <class T>
AlignedBuffer<T> make_aligned_buffer(std::size_t count, const char* name) {
  try {
    return AlignedBuffer<T>(count);
  } catch (const std::bad_alloc&) {
    throw std::runtime_error(std::string("malloc failure (") + name + ")");
  }
}

} // namespace detail

namespace {

std::unique_ptr<fft::Backend> create_fft_backend(std::string_view name) {
  if (name == "fftw") {
    return fft::create_fftw_backend();
  }
  if (name == "pocketfft") {
    return fft::create_pocketfft_backend();
  }
  if (name == "vkfft" || name == "vkfft-vulkan") {
#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)
    return fft::create_vkfft_vulkan_backend();
#else
    throw std::runtime_error("fft_backend 'vkfft-vulkan' requires a build with NEO_DFTTEST_ENABLE_VKFFT_VULKAN=ON");
#endif
  }
  throw std::runtime_error("unsupported FFT backend");
}

int logical_cpu_threads() noexcept {
  const unsigned int hardware_threads = std::thread::hardware_concurrency();
  return hardware_threads > 0 ? static_cast<int>(hardware_threads) : 4;
}

int conservative_cpu_budget() noexcept {
  return std::clamp(logical_cpu_threads() / 2, 1, 16);
}

int auto_worker_threads(int host_concurrency) noexcept {
  return std::clamp(conservative_cpu_budget() / std::max(1, host_concurrency), 1, 16);
}

} // namespace

class DftTestEngine::Impl {
public:
  Impl() = default;
  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;

  ~Impl() {
    if (fft_ && fft_->loaded()) {
      ds::HostGlobalLockGuard lock("fftw", host_locks_);
      state_.fft.forward = fft::Plan();
      state_.fft.inverse = fft::Plan();
      fft_->unload();
    }
  }

  void initialize(
    const ds::VideoInputInfo& input,
    const ds::ParamValues* params,
    ds::HostGlobalLockCallbacks host_locks
  ) {
    input_ = input;
    host_locks_ = host_locks;
    const ds::ParamValues empty_params{};
    const ds::ParamValues& values = params ? *params : empty_params;

    state_.format.bits_per_sample = ds::bits_per_sample(input.format.sample_format);
    state_.format.bytes_per_sample = ds::bytes_per_sample(input.format.sample_format);
    state_.format.integer = input.format.sample_format != ds::SampleFormat::Float32;
    state_.format.num_planes = input.format.plane_count;
    state_.format.width = input.width;
    state_.format.height = input.height;
    state_.format.subsampling_h = input.format.subsampling_h;
    state_.format.subsampling_w = input.format.subsampling_w;

    config_ = DfttestConfig::read(values, state_);
    fft_ = create_fft_backend(config_.fft_backend);
    fft_->load();
    state_.fft.backend = fft_.get();
    apply_thread_policy(1);
    if (config_.fft_threads > 1 && state_.fft.backend->has_threading()) {
      state_.fft.backend->set_thread_count(config_.fft_threads);
    }
    config_.validate(input, state_);
    config_.configure_planes(values, state_);
    configure_batch_policy();

    state_.kernels = selectFunctions(
      static_cast<unsigned>(config_.ftype),
      static_cast<unsigned>(config_.opt),
      state_.format,
      state_.block
    );
    executor_ = create_cpu_dft_executor(static_cast<unsigned>(config_.opt), state_.format);

    if (state_.format.integer) {
      state_.sample.multiplier = static_cast<float>(1 << (state_.format.bits_per_sample - 8));
      state_.sample.divisor = 1.0f / state_.sample.multiplier;
      state_.sample.peak = (1 << state_.format.bits_per_sample) - 1;
    }

    if (config_.ftype != 0) {
      state_.block.f0_beta = 1.0f;
    }

    configure_geometry();
    create_fft_plans();
    initialize_sigma_profile();
    prepare_noise_points();
    resize_thread_storage(state_.block.worker_threads * config_.fft_threads * 16);
  }

  void request_frames(ds::VideoRequestContext& context) const {
    if (state_.block.temporal_size == 1) {
      context.request_frame(0, context.output_frame);
    } else {
      const int start = std::max(context.output_frame - state_.block.temporal_size / 2, 0);
      const int stop = std::min(context.output_frame + state_.block.temporal_size / 2, input_.num_frames - 1);
      for (int frame = start; frame <= stop; ++frame) {
        context.request_frame(0, frame);
      }
    }

    if (noise_points_.empty() || noise_profile_ready_) {
      return;
    }

    for (const NPInfo& point : noise_points_) {
      for (int z = 0; z < state_.block.temporal_size; ++z) {
        context.request_frame(0, point.fn + z);
      }
    }
  }

  void process_frame(int n, ds::VideoFrameProvider& provider, ds::MutableVideoFrameView dst) {
    build_noise_profile_if_needed(provider);

    ThreadSlot slot(*this);
    auto current = provider.get(0, n);
    if (!current.has_value()) {
      throw std::runtime_error(current.error().message);
    }

    if (state_.block.temporal_size == 1) {
      process_spatial_frame(slot.id(), current.value().frame, dst);
      return;
    }

    process_temporal_frame(slot.id(), n, current.value().frame, provider, dst);
  }

  int cache_hints(int cachehints, int frame_range, int default_response) {
    if (cachehints == CACHE_INFORM_NUM_THREADS) {
      const int n_threads = std::max(frame_range, 0);
      std::lock_guard<std::mutex> lock(thread_check_mutex_);
      apply_thread_policy(std::max(n_threads, 1));
      if (n_threads > static_cast<int>(thread_id_store_.size())) {
        resize_thread_storage_unlocked(n_threads);
      }
      return 0;
    }

    return default_response;
  }

private:
  class ThreadSlot {
  public:
    explicit ThreadSlot(Impl& owner)
      : owner_(owner),
        id_(owner.acquire_thread_slot()) {}

    ThreadSlot(const ThreadSlot&) = delete;
    ThreadSlot& operator=(const ThreadSlot&) = delete;

    ~ThreadSlot() {
      owner_.release_thread_slot(id_);
    }

    unsigned int id() const {
      return id_;
    }

  private:
    Impl& owner_;
    unsigned int id_;
  };

  void apply_thread_policy(int host_concurrency) {
    if (config_.fft_threads_auto) {
      config_.fft_threads = 1;
    }

    if (config_.worker_threads_auto) {
      state_.block.worker_threads =
        (!config_.fft_threads_auto && config_.fft_threads > 1)
          ? 1
          : auto_worker_threads(host_concurrency);
    }

    state_.block.worker_threads = std::clamp(state_.block.worker_threads, 1, 16);
    config_.fft_threads = std::max(config_.fft_threads, 1);
  }

  void configure_geometry() {
    state_.derived.block_area = state_.block.spatial_size * state_.block.spatial_size;
    state_.derived.block_volume = state_.derived.block_area * state_.block.temporal_size;
    state_.derived.complex_count = (state_.block.spatial_size / 2 + 1) * state_.block.spatial_size * state_.block.temporal_size;
    state_.derived.coefficient_count = state_.derived.complex_count * 2;
    state_.derived.transform_type =
      config_.tmode * 4 + (state_.block.temporal_size > 1 ? 2 : 0) + config_.smode;
    state_.derived.spatial_center = state_.block.spatial_size / 2;
    state_.derived.custom_f0_beta = std::abs(state_.block.f0_beta - 1.0f) >= 0.00005f;
    state_.derived.step = (state_.derived.transform_type & 1) ? state_.block.spatial_size - state_.block.spatial_overlap : 1;

    for (int plane = 0; plane < state_.format.num_planes; plane++) {
      const int shift_w = (plane == 1 || plane == 2) ? state_.format.subsampling_w : 0;
      const int shift_h = (plane == 1 || plane == 2) ? state_.format.subsampling_h : 0;
      state_.planes.width[plane] = state_.format.width >> shift_w;
      state_.planes.height[plane] = state_.format.height >> shift_h;
      const int width = state_.planes.width[plane];
      const int height = state_.planes.height[plane];

      if (config_.smode == 0) {
        const int ae = (state_.block.spatial_size >> 1) << 1;
        state_.planes.pad_width[plane] = width + ae;
        state_.planes.pad_height[plane] = height + ae;
        state_.planes.e_height[plane] = height;
      } else {
        const int ae = std::max(state_.block.spatial_size - state_.block.spatial_overlap, state_.block.spatial_overlap) * 2;
        state_.planes.pad_width[plane] = width + EXTRA(width, state_.block.spatial_size) + ae;
        state_.planes.pad_height[plane] = height + EXTRA(height, state_.block.spatial_size) + ae;
        state_.planes.e_height[plane] =
          (state_.planes.pad_height[plane] - state_.block.spatial_overlap) / (state_.block.spatial_size - state_.block.spatial_overlap) *
          (state_.block.spatial_size - state_.block.spatial_overlap);
      }

      state_.planes.pad_stride[plane] =
        ((state_.planes.pad_width[plane] * state_.format.bytes_per_sample - 1) | (FRAME_ALIGN - 1)) + 1;
      state_.planes.pad_block_size[plane] = state_.planes.pad_stride[plane] * state_.planes.pad_height[plane];
      state_.planes.e_stride[plane] = ((state_.planes.pad_width[plane] * static_cast<int>(sizeof(float)) - 1) | (FRAME_ALIGN - 1)) + 1;
    }
  }

  void configure_batch_policy() {
    state_.batch_policy = make_cpu_dft_batch_policy(state_.block, state_.fft.backend->capabilities());
  }

  void create_fft_plans() {
    state_.coefficients.window =
      detail::make_aligned_buffer<float>(state_.derived.block_volume + 7, "hw");
    createWindow(
      DFTMutableFloatSpan{state_.coefficients.window.data(), state_.derived.block_volume},
      config_.tmode,
      config_.smode,
      state_.block,
      state_.derived
    );

    const int real_stride = dft_scratch_real_stride(state_.derived);
    const int complex_stride = dft_scratch_complex_stride(state_.derived);
    const fft::BatchLayout fft_layout{
      dft_fft_batch_capacity(state_.batch_policy),
      real_stride,
      complex_stride
    };

    auto dftgr = detail::make_aligned_buffer<float>(real_stride, "dftgr");
    state_.coefficients.window_dft =
      detail::make_aligned_buffer<fft::Complex>(complex_stride, "dftgc");

    {
      ds::HostGlobalLockGuard fftw_lock("fftw", host_locks_);
      const fft::TransformShape shape = state_.block.temporal_size > 1
        ? fft::TransformShape{3, {state_.block.temporal_size, state_.block.spatial_size, state_.block.spatial_size}}
        : fft::TransformShape{2, {state_.block.spatial_size, state_.block.spatial_size, 1}};
      state_.fft.forward = state_.fft.backend->make_plan(
        fft::TransformDirection::r2c,
        shape,
        fft_layout,
        dftgr.data(),
        state_.coefficients.window_dft.data(),
        fft::kPatientDestroyInputPlanOptions
      );
      state_.fft.inverse = state_.fft.backend->make_plan(
        fft::TransformDirection::c2r,
        shape,
        fft_layout,
        dftgr.data(),
        state_.coefficients.window_dft.data(),
        fft::kPatientDestroyInputPlanOptions
      );
    }

    float wscale = 0.0f;
    const float* hw_t = state_.coefficients.window.data();
    float* dftgr_t = dftgr.data();
    for (int s = 0; s < state_.block.temporal_size; s++) {
      for (int y = 0; y < state_.block.spatial_size; y++) {
        for (int x = 0; x < state_.block.spatial_size; x++) {
          dftgr_t[x] = 255.0f * hw_t[x];
          wscale += hw_t[x] * hw_t[x];
        }
        hw_t += state_.block.spatial_size;
        dftgr_t += state_.block.spatial_size;
      }
    }

    state_.fft.backend->submit_r2c(
      state_.fft.forward,
      fft::single_r2c_batch(dftgr.data(), state_.coefficients.window_dft.data(), real_stride, complex_stride)
    ).wait();
    wscale_ = 1.0f / wscale;
  }

  void initialize_sigma_profile() {
    const float wscalef = (config_.ftype < 2) ? wscale_ : 1.0f;

    state_.coefficients.sigmas =
      detail::make_aligned_buffer<float>(state_.derived.coefficient_count + 7, "sigmas");
    state_.coefficients.sigmas2 =
      detail::make_aligned_buffer<float>(state_.derived.coefficient_count + 7, "sigmas2");
    state_.coefficients.pmins =
      detail::make_aligned_buffer<float>(state_.derived.coefficient_count + 7, "pmins");
    state_.coefficients.pmaxs =
      detail::make_aligned_buffer<float>(state_.derived.coefficient_count + 7, "pmaxs");

    if (!config_.slocation.empty() || !config_.ssx.empty() || !config_.ssy.empty() || !config_.sst.empty()) {
      initialize_spatially_varying_sigmas(wscalef);
    } else {
      for (int i = 0; i < state_.derived.coefficient_count; i++) {
        state_.coefficients.sigmas[i] = config_.sigma / wscalef;
      }
    }

    for (int i = 0; i < state_.derived.coefficient_count; i++) {
      state_.coefficients.sigmas2[i] = config_.sigma2 / wscalef;
      state_.coefficients.pmins[i] = config_.pmin / wscale_;
      state_.coefficients.pmaxs[i] = config_.pmax / wscale_;
    }
  }

  void initialize_spatially_varying_sigmas(float wscalef) {
    int ndim = 3;
    if (state_.block.temporal_size == 1) {
      ndim -= 1;
    }
    if (state_.block.spatial_size == 1) {
      ndim -= 2;
    }

    const float ndiv = 1.0f / static_cast<float>(ndim);
    std::vector<float> tdata;
    std::vector<float> sydata;
    std::vector<float> sxdata;

    if (!config_.slocation.empty()) {
      tdata = parseSigmaLocation(config_.slocation, config_.sigma, config_.ssystem ? 1.0f : ndiv);
      sydata = parseSigmaLocation(config_.slocation, config_.sigma, config_.ssystem ? 1.0f : ndiv);
      sxdata = parseSigmaLocation(config_.slocation, config_.sigma, config_.ssystem ? 1.0f : ndiv);
    } else {
      tdata = parseSigmaLocation(config_.sst, config_.sigma, ndiv);
      sydata = parseSigmaLocation(config_.ssy, config_.sigma, ndiv);
      sxdata = parseSigmaLocation(config_.ssx, config_.sigma, ndiv);
    }

    const int tcnt = static_cast<int>(tdata.size() / 2);
    const int sycnt = static_cast<int>(sydata.size() / 2);
    const int sxcnt = static_cast<int>(sxdata.size() / 2);

    const int cpx = state_.block.spatial_size / 2 + 1;
    float pft = 0.0f;
    float pfy = 0.0f;
    float pfx = 0.0f;

    for (int z = 0; z < state_.block.temporal_size; z++) {
      const float tval = getSVal(z, state_.block.temporal_size, DFTConstFloatSpan{tdata.data(), tcnt * 2}, pft);
      for (int y = 0; y < state_.block.spatial_size; y++) {
        const float syval = getSVal(y, state_.block.spatial_size, DFTConstFloatSpan{sydata.data(), sycnt * 2}, pfy);
        for (int x = 0; x < cpx; x++) {
          const float sxval = getSVal(x, state_.block.spatial_size, DFTConstFloatSpan{sxdata.data(), sxcnt * 2}, pfx);
          float val = 0.0f;

          if (config_.ssystem) {
            const float dw = std::sqrt((pft * pft + pfy * pfy + pfx * pfx) / static_cast<float>(ndim));
            val = interp(dw, DFTConstFloatSpan{tdata.data(), tcnt * 2});
          } else {
            val = tval * syval * sxval;
          }

          const int pos = ((z * state_.block.spatial_size + y) * cpx + x) * 2;
          state_.coefficients.sigmas[pos] = state_.coefficients.sigmas[pos + 1] = val / wscalef;
        }
      }
    }
  }

  void prepare_noise_points() {
    if (config_.nlocation.empty() || config_.ftype >= 2) {
      noise_profile_ready_ = true;
      return;
    }

    for (std::size_t i = 0; i < config_.nlocation.size(); i += 4) {
      const int fn = config_.nlocation[i + 0];
      const int plane = config_.nlocation[i + 1];
      const int y = config_.nlocation[i + 2];
      const int x = config_.nlocation[i + 3];

      if (fn < 0 || fn > input_.num_frames - state_.block.temporal_size) {
        throw std::runtime_error("invalid frame number in nlocation (" + std::to_string(fn) + ")");
      }
      if (plane < 0 || plane >= state_.format.num_planes) {
        throw std::runtime_error("invalid plane number in nlocation (" + std::to_string(plane) + ")");
      }

      const int height = state_.format.height >> (plane > 0 ? state_.format.subsampling_h : 0);
      if (y < 0 || y > height - state_.block.spatial_size) {
        throw std::runtime_error("invalid y pos in nlocation (" + std::to_string(y) + ")");
      }

      const int width = state_.format.width >> (plane > 0 ? state_.format.subsampling_w : 0);
      if (x < 0 || x > width - state_.block.spatial_size) {
        throw std::runtime_error("invalid x pos in nlocation (" + std::to_string(x) + ")");
      }
      if (noise_points_.size() >= 500) {
        throw std::runtime_error("maximum number of entries in nlocation is 500");
      }

      noise_points_.push_back(NPInfo{fn, plane, y, x});
    }
  }

  void build_noise_profile_if_needed(ds::VideoFrameProvider& provider) {
    if (noise_profile_ready_) {
      return;
    }

    std::lock_guard<std::mutex> lock(noise_profile_mutex_);
    if (noise_profile_ready_) {
      return;
    }

    std::memset(
      state_.coefficients.sigmas.data(),
      0,
      static_cast<std::size_t>(state_.derived.coefficient_count) * sizeof(float)
    );

    auto hw2 = detail::make_aligned_buffer<float>(state_.derived.block_volume + 7, "hw2");
    createWindow(
      DFTMutableFloatSpan{hw2.data(), state_.derived.block_volume},
      0,
      0,
      state_.block,
      state_.derived
    );

    const int real_stride = dft_scratch_real_stride(state_.derived);
    const int complex_stride = dft_scratch_complex_stride(state_.derived);
    auto dftr = detail::make_aligned_buffer<float>(real_stride, "dftr");
    auto dftgc2 = detail::make_aligned_buffer<fft::Complex>(complex_stride, "dftgc2");

    float wscale2 = 0.0f;
    int w = 0;
    for (int s = 0; s < state_.block.temporal_size; s++) {
      for (int y = 0; y < state_.block.spatial_size; y++) {
        for (int x = 0; x < state_.block.spatial_size; x++, w++) {
          dftr[w] = 255.0f * hw2[w];
          wscale2 += hw2[w] * hw2[w];
        }
      }
    }
    wscale2 = 1.0f / wscale2;
    state_.fft.backend->submit_r2c(
      state_.fft.forward,
      fft::single_r2c_batch(dftr.data(), dftgc2.data(), real_stride, complex_stride)
    ).wait();

    auto dftc = detail::make_aligned_buffer<fft::Complex>(complex_stride, "dftc");
    auto dftc2 = detail::make_aligned_buffer<fft::Complex>(complex_stride, "dftc2");

    for (const NPInfo& point : noise_points_) {
      for (int z = 0; z < state_.block.temporal_size; z++) {
        auto frame = provider.get(0, point.fn + z);
        if (!frame.has_value()) {
          throw std::runtime_error(frame.error().message);
        }

        const auto plane = engine::read_plane(frame.value().frame, point.b, state_);
        const int stride_elements = plane.stride_elements(state_.format.bytes_per_sample);

        if (state_.format.bytes_per_sample == 1) {
          const auto* srcp = plane.typed_at<std::uint8_t>(point.x, point.y, state_.format.bytes_per_sample);
          load_windowed_block_scalar(
            DFTConstSampleBlock<std::uint8_t>{srcp, stride_elements, state_.block.spatial_size},
            DFTConstFloatSpan{hw2.data() + state_.derived.block_area * z, state_.derived.block_area},
            DFTMutableFloatSpan{dftr.data() + state_.derived.block_area * z, state_.derived.block_area},
            state_.sample.divisor
          );
        } else if (state_.format.bytes_per_sample == 2) {
          const auto* srcp = plane.typed_at<std::uint16_t>(point.x, point.y, state_.format.bytes_per_sample);
          load_windowed_block_scalar(
            DFTConstSampleBlock<std::uint16_t>{srcp, stride_elements, state_.block.spatial_size},
            DFTConstFloatSpan{hw2.data() + state_.derived.block_area * z, state_.derived.block_area},
            DFTMutableFloatSpan{dftr.data() + state_.derived.block_area * z, state_.derived.block_area},
            state_.sample.divisor
          );
        } else {
          const auto* srcp = plane.typed_at<float>(point.x, point.y, state_.format.bytes_per_sample);
          load_windowed_block_scalar(
            DFTConstSampleBlock<float>{srcp, stride_elements, state_.block.spatial_size},
            DFTConstFloatSpan{hw2.data() + state_.derived.block_area * z, state_.derived.block_area},
            DFTMutableFloatSpan{dftr.data() + state_.derived.block_area * z, state_.derived.block_area},
            state_.sample.divisor
          );
        }
      }

      state_.fft.backend->submit_r2c(
        state_.fft.forward,
        fft::single_r2c_batch(dftr.data(), dftc.data(), real_stride, complex_stride)
      ).wait();

      if (state_.block.zero_mean) {
        remove_mean_scalar(
          DFTMutableFloatSpan{complex_float_data(dftc.data()), state_.derived.coefficient_count},
          DFTConstFloatSpan{complex_float_data(dftgc2.data()), state_.derived.coefficient_count},
          DFTMutableFloatSpan{complex_float_data(dftc2.data()), state_.derived.coefficient_count}
        );
      }

      for (int h = 0; h < state_.derived.coefficient_count; h += 2) {
        const float* dftc_f = complex_float_data(dftc.data());
        const float psd = dftc_f[h] * dftc_f[h] + dftc_f[h + 1] * dftc_f[h + 1];
        state_.coefficients.sigmas[h] += psd;
        state_.coefficients.sigmas[h + 1] += psd;
      }
    }

    const float scale = 1.0f / static_cast<float>(noise_points_.size());
    for (int h = 0; h < state_.derived.coefficient_count; h++) {
      state_.coefficients.sigmas[h] *= scale * (wscale2 / wscale_) * config_.alpha;
    }

    noise_profile_ready_ = true;
  }

  void process_spatial_frame(
    unsigned int thread_id,
    const ds::VideoFrameView& src,
    ds::MutableVideoFrameView& dst
  ) {
    const DFTKernelContext kernel_context = make_kernel_context(state_);

    for (int plane = 0; plane < state_.format.num_planes; plane++) {
      const auto src_plane = engine::read_plane(src, plane, state_);
      const auto dst_plane = engine::write_plane(dst, plane, state_);

      if (state_.planes.process[plane] == 3) {
        auto pad = detail::make_aligned_buffer<unsigned char>(
          state_.planes.pad_block_size[plane] * state_.block.temporal_size,
          "pad0"
        );
        const DFTMutablePlaneBytes padded_plane{pad.data(), state_.planes.pad_stride[plane]};
        executor_->copy_pad(
          plane,
          DFTPlaneBytes{src_plane.data, src_plane.stride_bytes},
          padded_plane,
          kernel_context
        );
        executor_->process_spatial(
          thread_id,
          plane,
          DFTPlaneBytes{padded_plane.data, padded_plane.stride_bytes},
          DFTMutablePlaneBytes{dst_plane.data, dst_plane.stride_bytes},
          kernel_context
        );
      } else if (state_.planes.process[plane] == 2) {
        engine::copy_plane_rows(src_plane, dst_plane);
      }
    }
  }

  void process_temporal_frame(
    unsigned int thread_id,
    int n,
    const ds::VideoFrameView& current,
    ds::VideoFrameProvider& provider,
    ds::MutableVideoFrameView& dst
  ) {
    const int pos = state_.block.temporal_size / 2;
    const DFTKernelContext kernel_context = make_kernel_context(state_);

    for (int plane = 0; plane < state_.format.num_planes; plane++) {
      const auto src0_plane = engine::read_plane(current, plane, state_);
      const auto dst_plane = engine::write_plane(dst, plane, state_);

      if (state_.planes.process[plane] == 3) {
        auto pad0 = detail::make_aligned_buffer<unsigned char>(
          state_.planes.pad_block_size[plane] * state_.block.temporal_size,
          "pad0"
        );

        for (int i = 0; i < state_.block.temporal_size; i++) {
          const int fn = i + n - pos;
          const int fn_real = std::min(std::max(fn, 0), input_.num_frames - 1);
          auto src_frame = provider.get(0, fn_real);
          if (!src_frame.has_value()) {
            throw std::runtime_error(src_frame.error().message);
          }

          const auto src_plane = engine::read_plane(src_frame.value().frame, plane, state_);
          auto* pad = pad0.data() + state_.planes.pad_block_size[plane] * i;
          executor_->copy_pad(
            plane,
            DFTPlaneBytes{src_plane.data, src_plane.stride_bytes},
            DFTMutablePlaneBytes{pad, state_.planes.pad_stride[plane]},
            kernel_context
          );
        }

        executor_->process_temporal(
          thread_id,
          plane,
          DFTPlaneBytes{pad0.data(), state_.planes.pad_stride[plane]},
          DFTMutablePlaneBytes{dst_plane.data, dst_plane.stride_bytes},
          pos,
          kernel_context
        );
      } else if (state_.planes.process[plane] == 2) {
        engine::copy_plane_rows(src0_plane, dst_plane);
      }
    }
  }

  unsigned int acquire_thread_slot() {
    std::lock_guard<std::mutex> lock(thread_check_mutex_);
    auto it = std::find(thread_id_store_.begin(), thread_id_store_.end(), 0);
    if (it == thread_id_store_.end()) {
      throw std::runtime_error("neo_dfttest: max thread limit reached.");
    }

    const auto thread_id = static_cast<unsigned int>(std::distance(thread_id_store_.begin(), it));
    thread_id_store_[thread_id] = 1;
    ensure_thread_buffers_unlocked(thread_id);
    return thread_id;
  }

  void release_thread_slot(unsigned int thread_id) {
    std::lock_guard<std::mutex> lock(thread_check_mutex_);
    thread_id_store_[thread_id] = 0;
  }

  void resize_thread_storage(int count) {
    std::lock_guard<std::mutex> lock(thread_check_mutex_);
    resize_thread_storage_unlocked(count);
  }

  void resize_thread_storage_unlocked(int count) {
    thread_id_store_.resize(static_cast<std::size_t>(count), 0);
    state_.scratch.slots.resize(static_cast<std::size_t>(count));
  }

  void ensure_thread_buffers_unlocked(unsigned int thread_id) {
    DFTThreadScratchSlot& slot = state_.scratch.slots[thread_id];
    const int scratch_slots = dft_fft_scratch_slots(state_.batch_policy);
    if (!slot.ebuff) {
      slot.ebuff = detail::make_aligned_buffer<float>(
        static_cast<std::size_t>(state_.planes.e_stride[0]) * state_.planes.pad_height[0],
        "thread ebuff"
      );

      if (state_.block.dither_mode > 0) {
        slot.dither_buffer = detail::make_aligned_buffer<float>(
          static_cast<std::size_t>(2) * state_.format.width,
          "dither buffer"
        );
        if (state_.block.dither_mode > 1) {
          slot.rng = std::make_unique<std::mt19937>(std::random_device{}());
        }
      }
    }

    if (slot.fft_scratch_slots >= scratch_slots) {
      return;
    }

    slot.dftr = detail::make_aligned_buffer<float>(
      static_cast<std::size_t>(dft_scratch_real_stride(state_.derived)) * scratch_slots,
      "thread dftr"
    );
    slot.dftc = detail::make_aligned_buffer<fft::Complex>(
      static_cast<std::size_t>(dft_scratch_complex_stride(state_.derived)) * scratch_slots,
      "thread dftc"
    );
    slot.dftc2 = detail::make_aligned_buffer<fft::Complex>(
      static_cast<std::size_t>(dft_scratch_complex_stride(state_.derived)) * scratch_slots,
      "thread dftc2"
    );
    slot.fft_scratch_slots = scratch_slots;
  }

  ds::VideoInputInfo input_{};
  ds::HostGlobalLockCallbacks host_locks_{};
  DFTTestData state_{};
  DfttestConfig config_{};
  std::unique_ptr<fft::Backend> fft_;
  std::unique_ptr<DftExecutor> executor_;
  float wscale_ = 1.0f;
  std::vector<NPInfo> noise_points_;
  bool noise_profile_ready_ = false;
  std::mutex noise_profile_mutex_;
  std::mutex thread_check_mutex_;
  std::vector<int> thread_id_store_;
};

DftTestEngine::DftTestEngine()
  : impl_(std::make_unique<Impl>()) {}

DftTestEngine::~DftTestEngine() = default;

DftTestEngine::DftTestEngine(DftTestEngine&&) noexcept = default;

DftTestEngine& DftTestEngine::operator=(DftTestEngine&&) noexcept = default;

void DftTestEngine::initialize(
  const ds::VideoInputInfo& input,
  const ds::ParamValues* params,
  ds::HostGlobalLockCallbacks host_locks
) {
  impl_->initialize(input, params, host_locks);
}

void DftTestEngine::request_frames(ds::VideoRequestContext& context) const {
  impl_->request_frames(context);
}

void DftTestEngine::process_frame(
  int n,
  ds::VideoFrameProvider& provider,
  ds::MutableVideoFrameView dst
) {
  impl_->process_frame(n, provider, dst);
}

int DftTestEngine::cache_hints(int cachehints, int frame_range, int default_response) {
  return impl_->cache_hints(cachehints, frame_range, default_response);
}

} // namespace neo_dfttest
