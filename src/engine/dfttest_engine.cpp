/*
 * Copyright 2020 Xinyue Lu
 *
 * DFTTest engine implementation.
 */

#include "engine/dfttest_engine.hpp"

#include <avisynth.h>

#include <dualsynth/global_lock.hpp>
#include <dualsynth/video_filter.hpp>

#include "dft_common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <locale>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace neo_dfttest {

namespace detail {

template <class T>
T unwrap(ds::Result<T> result) {
  if (!result.has_value()) {
    throw std::runtime_error(result.error().message);
  }
  return std::move(result.value());
}

inline bool has_param(const ds::ParamValues& params, const std::string& name) {
  return std::ranges::any_of(params.entries, [&](const ds::ParamEntry& entry) {
    return entry.name == name;
  });
}

inline int read_int(const ds::ParamValues& params, const std::string& name, int default_value) {
  return unwrap(params.get_int(name, default_value));
}

inline float read_float(const ds::ParamValues& params, const std::string& name, float default_value) {
  return static_cast<float>(unwrap(params.get_double(name, default_value)));
}

inline bool read_bool(const ds::ParamValues& params, const std::string& name, bool default_value) {
  return unwrap(params.get_bool(name, default_value));
}

inline std::vector<int> read_int_array(
  const ds::ParamValues& params,
  const std::string& name,
  std::vector<std::int64_t> default_value = {}
) {
  const auto values = unwrap(params.get_int_array(name, std::move(default_value)));
  std::vector<int> output;
  output.reserve(values.size());
  for (const std::int64_t value : values) {
    if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max()) {
      throw std::runtime_error("parameter '" + name + "' contains an integer outside int range");
    }
    output.push_back(static_cast<int>(value));
  }
  return output;
}

inline std::vector<float> read_float_array(
  const ds::ParamValues& params,
  const std::string& name,
  std::vector<double> default_value = {}
) {
  const auto values = unwrap(params.get_double_array(name, std::move(default_value)));
  std::vector<float> output;
  output.reserve(values.size());
  for (const double value : values) {
    output.push_back(static_cast<float>(value));
  }
  return output;
}

inline void framecpy(
  unsigned char* dst_ptr,
  int dst_stride,
  const unsigned char* src_ptr,
  int src_stride,
  int width_byte,
  int height
) {
  if (src_stride == dst_stride) {
    std::memcpy(dst_ptr, src_ptr, static_cast<std::size_t>(dst_stride) * static_cast<std::size_t>(height));
    return;
  }

  for (int h = 0; h < height; h++) {
    std::memcpy(dst_ptr, src_ptr, static_cast<std::size_t>(width_byte));
    dst_ptr += dst_stride;
    src_ptr += src_stride;
  }
}

struct AlignedDeleter {
  void operator()(void* p) const {
    _aligned_free(p);
  }
};

template <class T>
using AlignedPtr = std::unique_ptr<T, AlignedDeleter>;

inline unsigned char* writable_plane_data(ds::MutableVideoFrameView& frame, int plane) {
  return static_cast<unsigned char*>(frame.plane(plane).data);
}

inline const unsigned char* readable_plane_data(const ds::VideoFrameView& frame, int plane) {
  return static_cast<const unsigned char*>(frame.plane(plane).data);
}

inline int plane_stride(const ds::PlaneView& plane) {
  return static_cast<int>(plane.stride_bytes);
}

inline int plane_stride(const ds::MutablePlaneView& plane) {
  return static_cast<int>(plane.stride_bytes);
}

} // namespace detail

class DftTestEngine::Impl {
public:
  Impl() = default;
  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;

  ~Impl() {
    _aligned_free(state_.coefficients.window);
    _aligned_free(state_.coefficients.window_dft);
    _aligned_free(state_.coefficients.sigmas);
    _aligned_free(state_.coefficients.sigmas2);
    _aligned_free(state_.coefficients.pmins);
    _aligned_free(state_.coefficients.pmaxs);

    for (auto&& buf : state_.scratch.ebuff) {
      _aligned_free(buf);
    }
    for (auto&& buf : state_.scratch.dftr) {
      _aligned_free(buf);
    }
    for (auto&& buf : state_.scratch.dftc) {
      _aligned_free(buf);
    }
    for (auto&& buf : state_.scratch.dftc2) {
      _aligned_free(buf);
    }
    for (auto&& buf : state_.scratch.dither_buffers) {
      _aligned_free(buf);
    }

    if (fft_ && fft_->loaded()) {
      ds::HostGlobalLockGuard lock("fftw", host_locks_);
      fft_->destroy_plan(state_.fft.forward);
      fft_->destroy_plan(state_.fft.inverse);
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

    fft_ = fft::create_fftw_backend();
    fft_->load();
    state_.fft.backend = fft_.get();
    state_.format.bits_per_sample = ds::bits_per_sample(input.format.sample_format);
    state_.format.bytes_per_sample = ds::bytes_per_sample(input.format.sample_format);
    state_.format.integer = input.format.sample_format != ds::SampleFormat::Float32;
    state_.format.num_planes = input.format.plane_count;
    state_.format.width = input.width;
    state_.format.height = input.height;
    state_.format.subsampling_h = input.format.subsampling_h;
    state_.format.subsampling_w = input.format.subsampling_w;

    read_parameters(values);
    validate_parameters(input);
    configure_planes(values);

#ifndef ENABLE_PAR
    state_.block.worker_threads = 1;
#endif

    selectFunctions(static_cast<unsigned>(ftype_), static_cast<unsigned>(opt_), &state_);

    if (state_.format.integer) {
      state_.sample.multiplier = static_cast<float>(1 << (state_.format.bits_per_sample - 8));
      state_.sample.divisor = 1.0f / state_.sample.multiplier;
      state_.sample.peak = (1 << state_.format.bits_per_sample) - 1;
    }

    if (ftype_ != 0) {
      state_.block.f0_beta = 1.0f;
    }

    configure_geometry();
    create_fft_plans();
    initialize_sigma_profile();
    prepare_noise_points();
    resize_thread_storage(state_.block.worker_threads * fft_threads_ * 16);
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

  void read_parameters(const ds::ParamValues& values) {
    ftype_ = detail::read_int(values, "ftype", 0);
    sigma_ = detail::read_float(values, "sigma", 8.0f);
    sigma2_ = detail::read_float(values, "sigma2", 8.0f);
    pmin_ = detail::read_float(values, "pmin", 0.0f);
    pmax_ = detail::read_float(values, "pmax", 500.0f);
    smode_ = detail::read_int(values, "smode", 1);
    tmode_ = detail::read_int(values, "tmode", 0);
    opt_ = detail::read_int(values, "opt", 0);

    state_.block.spatial_size = detail::read_int(values, "sbsize", state_.block.spatial_size);
    state_.block.spatial_overlap = detail::read_int(values, "sosize", state_.block.spatial_overlap);
    state_.block.temporal_size = detail::read_int(values, "tbsize", state_.block.temporal_size);
    state_.block.temporal_overlap = detail::read_int(values, "tosize", state_.block.temporal_overlap);
    state_.block.spatial_window = detail::read_int(values, "swin", state_.block.spatial_window);
    state_.block.temporal_window = detail::read_int(values, "twin", state_.block.temporal_window);
    state_.block.spatial_beta = detail::read_float(values, "sbeta", state_.block.spatial_beta);
    state_.block.temporal_beta = detail::read_float(values, "tbeta", state_.block.temporal_beta);
    state_.block.zero_mean = detail::read_bool(values, "zmean", state_.block.zero_mean);
    state_.block.f0_beta = detail::read_float(values, "f0beta", state_.block.f0_beta);
    state_.block.worker_threads = detail::read_int(values, "threads", state_.block.worker_threads);
    state_.block.dither_mode = detail::read_int(values, "dither", state_.block.dither_mode);
    fft_threads_ = detail::read_int(values, "fft_threads", fft_threads_);
    if (fft_threads_ < 1) {
      fft_threads_ = 1;
    }

    if (smode_ == 0) {
      state_.block.spatial_overlap = 0;
    }
    if (tmode_ == 0) {
      state_.block.temporal_overlap = 0;
    }

    if (fft_threads_ > 1 && state_.fft.backend->has_threading()) {
      state_.fft.backend->set_thread_count(fft_threads_);
    }

    nlocation_ = detail::read_int_array(values, "nlocation");
    alpha_ = detail::read_float(values, "alpha", ftype_ == 0 ? 5.0f : 7.0f);
    slocation_ = detail::read_float_array(values, "slocation");
    ssx_ = detail::read_float_array(values, "ssx");
    ssy_ = detail::read_float_array(values, "ssy");
    sst_ = detail::read_float_array(values, "sst");
    ssystem_ = detail::read_int(values, "ssystem", 0);

    if (state_.block.worker_threads <= 0) {
      state_.block.worker_threads = 4;
    }
    if (state_.block.worker_threads > 16) {
      state_.block.worker_threads = 16;
    }
  }

  void configure_planes(const ds::ParamValues& values) {
    std::fill(std::begin(state_.planes.process), std::end(state_.planes.process), 2);

    if (detail::has_param(values, "planes")) {
      const auto planes = detail::read_int_array(values, "planes");
      for (const int plane : planes) {
        if (plane < 0 || plane >= state_.format.num_planes) {
          throw std::runtime_error("plane index out of range");
        }
        state_.planes.process[plane] = 3;
      }
      return;
    }

    state_.planes.process[0] = 3;
    state_.planes.process[1] = 3;
    state_.planes.process[2] = 3;
    state_.planes.process[3] = 2;
    state_.planes.process[0] = detail::read_int(values, "y", state_.planes.process[0]);
    state_.planes.process[1] = detail::read_int(values, "u", state_.planes.process[1]);
    state_.planes.process[2] = detail::read_int(values, "v", state_.planes.process[2]);
    state_.planes.process[3] = detail::read_int(values, "a", state_.planes.process[3]);
  }

  void validate_parameters(const ds::VideoInputInfo& input) const {
    if (input.width <= 0 || input.height <= 0) {
      throw std::runtime_error("only constant format input supported");
    }

    if (
      (state_.format.integer && state_.format.bits_per_sample > 16) ||
      (!state_.format.integer && state_.format.bits_per_sample != 32)
    ) {
      throw std::runtime_error("only 8-16 bit integer and 32 bit float input supported");
    }

    if (ftype_ < 0 || ftype_ > 4) {
      throw std::runtime_error("ftype must be 0, 1, 2, 3, or 4");
    }
    if (state_.block.spatial_size < 1) {
      throw std::runtime_error("sbsize must be greater than or equal to 1");
    }
    if (smode_ < 0 || smode_ > 1) {
      throw std::runtime_error("smode must be 0 or 1");
    }
    if (smode_ == 0 && !(state_.block.spatial_size & 1)) {
      throw std::runtime_error("sbsize must be odd when using smode=0");
    }
    if (state_.block.spatial_overlap < 0 || state_.block.spatial_overlap >= state_.block.spatial_size) {
      throw std::runtime_error("sosize must be between 0 and sbsize-1 (inclusive)");
    }
    if (state_.block.spatial_overlap > state_.block.spatial_size / 2 && state_.block.spatial_size % (state_.block.spatial_size - state_.block.spatial_overlap) != 0) {
      throw std::runtime_error(
        "spatial overlap greater than 50% requires that sbsize-sosize is a divisor of sbsize"
      );
    }
    if (state_.block.temporal_size < 1 || state_.block.temporal_size > 15) {
      throw std::runtime_error("tbsize must be between 1 and 15 (inclusive)");
    }
    if (tmode_ != 0) {
      throw std::runtime_error("tmode must be 0. tmode=1 is not implemented");
    }
    if (tmode_ == 0 && !(state_.block.temporal_size & 1)) {
      throw std::runtime_error("tbsize must be odd when using tmode=0");
    }
    if (state_.block.temporal_overlap < 0 || state_.block.temporal_overlap >= state_.block.temporal_size) {
      throw std::runtime_error("tosize must be between 0 and tbsize-1 (inclusive)");
    }
    if (state_.block.temporal_overlap > state_.block.temporal_size / 2 && state_.block.temporal_size % (state_.block.temporal_size - state_.block.temporal_overlap) != 0) {
      throw std::runtime_error(
        "temporal overlap greater than 50% requires that tbsize-tosize is a divisor of tbsize"
      );
    }
    if (state_.block.temporal_size > input.num_frames) {
      throw std::runtime_error("tbsize must be less than or equal to the number of frames in the clip");
    }
    if (state_.block.spatial_window < 0 || state_.block.spatial_window > 11) {
      throw std::runtime_error("swin must be between 0 and 11 (inclusive)");
    }
    if (state_.block.temporal_window < 0 || state_.block.temporal_window > 11) {
      throw std::runtime_error("twin must be between 0 and 11 (inclusive)");
    }
    if (nlocation_.size() & 3U) {
      throw std::runtime_error("the number of elements in nlocation must be a multiple of 4");
    }
    if (alpha_ <= 0.0f) {
      throw std::runtime_error("alpha must be greater than 0.0");
    }
    if (slocation_.size() & 1U) {
      throw std::runtime_error("the number of elements in slocation must be even");
    }
    if (ssx_.size() & 1U) {
      throw std::runtime_error("the number of elements in ssx must be even");
    }
    if (ssy_.size() & 1U) {
      throw std::runtime_error("the number of elements in ssy must be even");
    }
    if (sst_.size() & 1U) {
      throw std::runtime_error("the number of elements in sst must be even");
    }
    if (ssystem_ < 0 || ssystem_ > 1) {
      throw std::runtime_error("ssystem must be 0 or 1");
    }
    if (opt_ != 8 && (opt_ < 0 || opt_ > 3)) {
      throw std::runtime_error("opt must be 0, 1, 2, 3, or 8");
    }
  }

  void configure_geometry() {
    state_.derived.block_area = state_.block.spatial_size * state_.block.spatial_size;
    state_.derived.block_volume = state_.derived.block_area * state_.block.temporal_size;
    state_.derived.complex_count = (state_.block.spatial_size / 2 + 1) * state_.block.spatial_size * state_.block.temporal_size;
    state_.derived.coefficient_count = state_.derived.complex_count * 2;
    state_.derived.transform_type = tmode_ * 4 + (state_.block.temporal_size > 1 ? 2 : 0) + smode_;
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

      if (smode_ == 0) {
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
      state_.planes.e_batch_size[plane] =
        ((state_.planes.e_height[plane] - 1) / state_.block.worker_threads / state_.derived.step + 1) * state_.derived.step;
    }
  }

  void create_fft_plans() {
    state_.coefficients.window = static_cast<float*>(_aligned_malloc((state_.derived.block_volume + 7) * sizeof(float), FRAME_ALIGN));
    if (!state_.coefficients.window) {
      throw std::runtime_error("malloc failure (hw)");
    }
    createWindow(state_.coefficients.window, tmode_, smode_, &state_);

    float* dftgr = static_cast<float*>(_aligned_malloc((state_.derived.block_volume + 7) * sizeof(float), FRAME_ALIGN));
    detail::AlignedPtr<float> dftgr_smart(dftgr);
    state_.coefficients.window_dft = static_cast<fft::Complex*>(
      _aligned_malloc((state_.derived.complex_count + 7) * sizeof(fft::Complex), FRAME_ALIGN)
    );
    if (!dftgr || !state_.coefficients.window_dft) {
      throw std::runtime_error("malloc failure (dftgr/dftgc)");
    }

    {
      ds::HostGlobalLockGuard fftw_lock("fftw", host_locks_);
      if (state_.block.temporal_size > 1) {
        state_.fft.forward = state_.fft.backend->plan_r2c_3d(
          state_.block.temporal_size,
          state_.block.spatial_size,
          state_.block.spatial_size,
          dftgr,
          state_.coefficients.window_dft,
          fft::kPatientDestroyInputPlanFlags
        );
        state_.fft.inverse = state_.fft.backend->plan_c2r_3d(
          state_.block.temporal_size,
          state_.block.spatial_size,
          state_.block.spatial_size,
          state_.coefficients.window_dft,
          dftgr,
          fft::kPatientDestroyInputPlanFlags
        );
      } else {
        state_.fft.forward = state_.fft.backend->plan_r2c_2d(
          state_.block.spatial_size,
          state_.block.spatial_size,
          dftgr,
          state_.coefficients.window_dft,
          fft::kPatientDestroyInputPlanFlags
        );
        state_.fft.inverse = state_.fft.backend->plan_c2r_2d(
          state_.block.spatial_size,
          state_.block.spatial_size,
          state_.coefficients.window_dft,
          dftgr,
          fft::kPatientDestroyInputPlanFlags
        );
      }
    }

    float wscale = 0.0f;
    const float* hw_t = state_.coefficients.window;
    float* dftgr_t = dftgr;
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

    state_.fft.backend->execute_r2c(state_.fft.forward, dftgr, state_.coefficients.window_dft);
    wscale_ = 1.0f / wscale;
  }

  void initialize_sigma_profile() {
    const float wscalef = (ftype_ < 2) ? wscale_ : 1.0f;

    state_.coefficients.sigmas = static_cast<float*>(_aligned_malloc((state_.derived.coefficient_count + 7) * sizeof(float), FRAME_ALIGN));
    state_.coefficients.sigmas2 = static_cast<float*>(_aligned_malloc((state_.derived.coefficient_count + 7) * sizeof(float), FRAME_ALIGN));
    state_.coefficients.pmins = static_cast<float*>(_aligned_malloc((state_.derived.coefficient_count + 7) * sizeof(float), FRAME_ALIGN));
    state_.coefficients.pmaxs = static_cast<float*>(_aligned_malloc((state_.derived.coefficient_count + 7) * sizeof(float), FRAME_ALIGN));
    if (!state_.coefficients.sigmas || !state_.coefficients.sigmas2 || !state_.coefficients.pmins || !state_.coefficients.pmaxs) {
      throw std::runtime_error("malloc failure (sigmas/sigmas2/pmins/pmaxs)");
    }

    if (!slocation_.empty() || !ssx_.empty() || !ssy_.empty() || !sst_.empty()) {
      initialize_spatially_varying_sigmas(wscalef);
    } else {
      for (int i = 0; i < state_.derived.coefficient_count; i++) {
        state_.coefficients.sigmas[i] = sigma_ / wscalef;
      }
    }

    for (int i = 0; i < state_.derived.coefficient_count; i++) {
      state_.coefficients.sigmas2[i] = sigma2_ / wscalef;
      state_.coefficients.pmins[i] = pmin_ / wscale_;
      state_.coefficients.pmaxs[i] = pmax_ / wscale_;
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
    int tcnt = 0;
    int sycnt = 0;
    int sxcnt = 0;
    float* tdata = nullptr;
    float* sydata = nullptr;
    float* sxdata = nullptr;

    if (!slocation_.empty()) {
      tdata = parseSigmaLocation(slocation_, tcnt, sigma_, ssystem_ ? 1.0f : ndiv);
      sydata = parseSigmaLocation(slocation_, sycnt, sigma_, ssystem_ ? 1.0f : ndiv);
      sxdata = parseSigmaLocation(slocation_, sxcnt, sigma_, ssystem_ ? 1.0f : ndiv);
    } else {
      tdata = parseSigmaLocation(sst_, tcnt, sigma_, ndiv);
      sydata = parseSigmaLocation(ssy_, sycnt, sigma_, ndiv);
      sxdata = parseSigmaLocation(ssx_, sxcnt, sigma_, ndiv);
    }

    std::unique_ptr<float[]> t_smart(tdata);
    std::unique_ptr<float[]> sx_smart(sxdata);
    std::unique_ptr<float[]> sy_smart(sydata);

    const int cpx = state_.block.spatial_size / 2 + 1;
    float pft = 0.0f;
    float pfy = 0.0f;
    float pfx = 0.0f;

    for (int z = 0; z < state_.block.temporal_size; z++) {
      const float tval = getSVal(z, state_.block.temporal_size, tdata, tcnt, pft);
      for (int y = 0; y < state_.block.spatial_size; y++) {
        const float syval = getSVal(y, state_.block.spatial_size, sydata, sycnt, pfy);
        for (int x = 0; x < cpx; x++) {
          const float sxval = getSVal(x, state_.block.spatial_size, sxdata, sxcnt, pfx);
          float val = 0.0f;

          if (ssystem_) {
            const float dw = std::sqrt((pft * pft + pfy * pfy + pfx * pfx) / static_cast<float>(ndim));
            val = interp(dw, tdata, tcnt);
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
    if (nlocation_.empty() || ftype_ >= 2) {
      noise_profile_ready_ = true;
      return;
    }

    for (std::size_t i = 0; i < nlocation_.size(); i += 4) {
      const int fn = nlocation_[i + 0];
      const int plane = nlocation_[i + 1];
      const int y = nlocation_[i + 2];
      const int x = nlocation_[i + 3];

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

    std::memset(state_.coefficients.sigmas, 0, static_cast<std::size_t>(state_.derived.coefficient_count) * sizeof(float));

    float* hw2 = static_cast<float*>(_aligned_malloc((state_.derived.block_volume + 7) * sizeof(float), FRAME_ALIGN));
    if (!hw2) {
      throw std::runtime_error("malloc failure (hw2)");
    }
    detail::AlignedPtr<float> hw2_smart(hw2);
    createWindow(hw2, 0, 0, &state_);

    float* dftr = static_cast<float*>(_aligned_malloc((state_.derived.block_volume + 7) * sizeof(float), FRAME_ALIGN));
    auto* dftgc2 = static_cast<fft::Complex*>(
      _aligned_malloc((state_.derived.complex_count + 7) * sizeof(fft::Complex), FRAME_ALIGN)
    );
    if (!dftr || !dftgc2) {
      throw std::runtime_error("malloc failure (dftr/dftgc2)");
    }
    detail::AlignedPtr<float> dftr_smart(dftr);
    detail::AlignedPtr<fft::Complex> dftgc2_smart(dftgc2);

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
    state_.fft.backend->execute_r2c(state_.fft.forward, dftr, dftgc2);

    auto* dftc = static_cast<fft::Complex*>(
      _aligned_malloc((state_.derived.complex_count + 7) * sizeof(fft::Complex), FRAME_ALIGN)
    );
    auto* dftc2 = static_cast<fft::Complex*>(
      _aligned_malloc((state_.derived.complex_count + 7) * sizeof(fft::Complex), FRAME_ALIGN)
    );
    if (!dftc || !dftc2) {
      throw std::runtime_error("malloc failure (dftc/dftc2)");
    }
    detail::AlignedPtr<fft::Complex> dftc_smart(dftc);
    detail::AlignedPtr<fft::Complex> dftc2_smart(dftc2);

    for (const NPInfo& point : noise_points_) {
      for (int z = 0; z < state_.block.temporal_size; z++) {
        auto frame = provider.get(0, point.fn + z);
        if (!frame.has_value()) {
          throw std::runtime_error(frame.error().message);
        }

        const ds::PlaneView& plane = frame.value().frame.plane(point.b);
        const int stride_elements = static_cast<int>(plane.stride_bytes) / state_.format.bytes_per_sample;

        if (state_.format.bytes_per_sample == 1) {
          const auto* srcp =
            static_cast<const std::uint8_t*>(plane.data) + stride_elements * point.y + point.x;
          proc0_c(srcp, hw2 + state_.derived.block_area * z, dftr + state_.derived.block_area * z, stride_elements, state_.block.spatial_size, state_.sample.divisor);
        } else if (state_.format.bytes_per_sample == 2) {
          const auto* srcp =
            static_cast<const std::uint16_t*>(plane.data) + stride_elements * point.y + point.x;
          proc0_c(srcp, hw2 + state_.derived.block_area * z, dftr + state_.derived.block_area * z, stride_elements, state_.block.spatial_size, state_.sample.divisor);
        } else {
          const auto* srcp =
            static_cast<const float*>(plane.data) + stride_elements * point.y + point.x;
          proc0_c(srcp, hw2 + state_.derived.block_area * z, dftr + state_.derived.block_area * z, stride_elements, state_.block.spatial_size, state_.sample.divisor);
        }
      }

      state_.fft.backend->execute_r2c(state_.fft.forward, dftr, dftc);

      if (state_.block.zero_mean) {
        removeMean_c(
          reinterpret_cast<float*>(dftc),
          reinterpret_cast<const float*>(dftgc2),
          state_.derived.coefficient_count,
          reinterpret_cast<float*>(dftc2)
        );
      }

      for (int h = 0; h < state_.derived.coefficient_count; h += 2) {
        const float* dftc_f = reinterpret_cast<float*>(dftc);
        const float psd = dftc_f[h] * dftc_f[h] + dftc_f[h + 1] * dftc_f[h + 1];
        state_.coefficients.sigmas[h] += psd;
        state_.coefficients.sigmas[h + 1] += psd;
      }
    }

    const float scale = 1.0f / static_cast<float>(noise_points_.size());
    for (int h = 0; h < state_.derived.coefficient_count; h++) {
      state_.coefficients.sigmas[h] *= scale * (wscale2 / wscale_) * alpha_;
    }

    noise_profile_ready_ = true;
  }

  void process_spatial_frame(
    unsigned int thread_id,
    const ds::VideoFrameView& src,
    ds::MutableVideoFrameView& dst
  ) {
    for (int plane = 0; plane < state_.format.num_planes; plane++) {
      const int height = state_.planes.height[plane];
      const int width = state_.planes.width[plane];
      const auto& src_plane = src.plane(plane);
      const auto& dst_plane = dst.plane(plane);
      const int src_stride = detail::plane_stride(src_plane);
      const int dst_stride = detail::plane_stride(dst_plane);
      const auto* src_ptr = detail::readable_plane_data(src, plane);
      auto* dst_ptr = detail::writable_plane_data(dst, plane);

      if (state_.planes.process[plane] == 3) {
        auto* pad = static_cast<unsigned char*>(_aligned_malloc(state_.planes.pad_block_size[plane] * state_.block.temporal_size, FRAME_ALIGN));
        if (!pad) {
          throw std::runtime_error("pad0 allocation failed.");
        }
        detail::AlignedPtr<unsigned char> pad0_smart(pad);
        state_.kernels.copy_pad(plane, src_ptr, src_stride, pad, &state_);
        state_.kernels.process_spatial(thread_id, plane, pad, dst_ptr, dst_stride, &state_);
      } else if (state_.planes.process[plane] == 2) {
        detail::framecpy(
          dst_ptr,
          dst_stride,
          src_ptr,
          src_stride,
          width * state_.format.bytes_per_sample,
          height
        );
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

    for (int plane = 0; plane < state_.format.num_planes; plane++) {
      const int height = state_.planes.height[plane];
      const int width = state_.planes.width[plane];
      const auto& src0_plane = current.plane(plane);
      const auto& dst_plane = dst.plane(plane);
      const int src0_stride = detail::plane_stride(src0_plane);
      const int dst_stride = detail::plane_stride(dst_plane);
      const auto* src0_ptr = detail::readable_plane_data(current, plane);
      auto* dst_ptr = detail::writable_plane_data(dst, plane);

      if (state_.planes.process[plane] == 3) {
        auto* pad0 = static_cast<unsigned char*>(_aligned_malloc(state_.planes.pad_block_size[plane] * state_.block.temporal_size, FRAME_ALIGN));
        if (!pad0) {
          throw std::runtime_error("pad0 allocation failed.");
        }
        detail::AlignedPtr<unsigned char> pad0_smart(pad0);

        for (int i = 0; i < state_.block.temporal_size; i++) {
          const int fn = i + n - pos;
          const int fn_real = std::min(std::max(fn, 0), input_.num_frames - 1);
          auto src_frame = provider.get(0, fn_real);
          if (!src_frame.has_value()) {
            throw std::runtime_error(src_frame.error().message);
          }

          const auto& src_plane = src_frame.value().frame.plane(plane);
          const int src_stride = detail::plane_stride(src_plane);
          const auto* src_ptr = detail::readable_plane_data(src_frame.value().frame, plane);
          auto* pad = pad0 + state_.planes.pad_block_size[plane] * i;
          state_.kernels.copy_pad(plane, src_ptr, src_stride, pad, &state_);
        }

        state_.kernels.process_temporal(thread_id, plane, pad0, dst_ptr, dst_stride, pos, &state_);
      } else if (state_.planes.process[plane] == 2) {
        detail::framecpy(
          dst_ptr,
          dst_stride,
          src0_ptr,
          src0_stride,
          width * state_.format.bytes_per_sample,
          height
        );
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
    state_.scratch.ebuff.resize(static_cast<std::size_t>(count), nullptr);
    state_.scratch.dftr.resize(static_cast<std::size_t>(count), nullptr);
    state_.scratch.dftc.resize(static_cast<std::size_t>(count), nullptr);
    state_.scratch.dftc2.resize(static_cast<std::size_t>(count), nullptr);
    state_.scratch.rngs.resize(static_cast<std::size_t>(count));
    state_.scratch.dither_buffers.resize(static_cast<std::size_t>(count), nullptr);
  }

  void ensure_thread_buffers_unlocked(unsigned int thread_id) {
    if (state_.scratch.ebuff[thread_id]) {
      return;
    }

    state_.scratch.ebuff[thread_id] = static_cast<float*>(
      _aligned_malloc(sizeof(float) * state_.planes.e_stride[0] * state_.planes.pad_height[0], FRAME_ALIGN)
    );
    state_.scratch.dftr[thread_id] = static_cast<float*>(
      _aligned_malloc(sizeof(float) * (((state_.derived.block_volume + 7) | 15) + 1) * state_.block.worker_threads, FRAME_ALIGN)
    );
    state_.scratch.dftc[thread_id] = static_cast<fft::Complex*>(
      _aligned_malloc(sizeof(fft::Complex) * (((state_.derived.complex_count + 7) | 15) + 1) * state_.block.worker_threads, FRAME_ALIGN)
    );
    state_.scratch.dftc2[thread_id] = static_cast<fft::Complex*>(
      _aligned_malloc(sizeof(fft::Complex) * (((state_.derived.complex_count + 7) | 15) + 1) * state_.block.worker_threads, FRAME_ALIGN)
    );

    if (!state_.scratch.ebuff[thread_id] || !state_.scratch.dftr[thread_id] || !state_.scratch.dftc[thread_id] || !state_.scratch.dftc2[thread_id]) {
      throw std::runtime_error("thread buffer allocation failed.");
    }

    if (state_.block.dither_mode > 0) {
      state_.scratch.dither_buffers[thread_id] = static_cast<float*>(
        _aligned_malloc(sizeof(float) * 2 * state_.format.width, FRAME_ALIGN)
      );
      if (!state_.scratch.dither_buffers[thread_id]) {
        throw std::runtime_error("dither buffer allocation failed.");
      }
      if (state_.block.dither_mode > 1) {
        state_.scratch.rngs[thread_id] = std::make_unique<std::mt19937>(std::random_device{}());
      }
    }
  }

  ds::VideoInputInfo input_{};
  ds::HostGlobalLockCallbacks host_locks_{};
  DFTTestData state_{};
  std::unique_ptr<fft::Backend> fft_;
  int fft_threads_ = 2;
  int ftype_ = 0;
  int smode_ = 1;
  int tmode_ = 0;
  int opt_ = 0;
  int ssystem_ = 0;
  float sigma_ = 8.0f;
  float sigma2_ = 8.0f;
  float pmin_ = 0.0f;
  float pmax_ = 500.0f;
  float alpha_ = 5.0f;
  float wscale_ = 1.0f;
  std::vector<int> nlocation_;
  std::vector<float> slocation_;
  std::vector<float> ssx_;
  std::vector<float> ssy_;
  std::vector<float> sst_;
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
