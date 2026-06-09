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
    _aligned_free(ep_.hw);
    _aligned_free(ep_.dftgc);
    _aligned_free(ep_.sigmas);
    _aligned_free(ep_.sigmas2);
    _aligned_free(ep_.pmins);
    _aligned_free(ep_.pmaxs);

    for (auto&& buf : ep_.ebuff) {
      _aligned_free(buf);
    }
    for (auto&& buf : ep_.dftr) {
      _aligned_free(buf);
    }
    for (auto&& buf : ep_.dftc) {
      _aligned_free(buf);
    }
    for (auto&& buf : ep_.dftc2) {
      _aligned_free(buf);
    }
    for (auto&& buf : ep_.d_buffs) {
      _aligned_free(buf);
    }

    if (fft_.library) {
      ds::HostGlobalLockGuard lock("fftw", host_locks_);
      if (ep_.ft) {
        fft_.fftwf_destroy_plan(ep_.ft);
      }
      if (ep_.fti) {
        fft_.fftwf_destroy_plan(ep_.fti);
      }
      fft_.free();
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

    fft_.load();
    ep_.fft = &fft_;
    ep_.vi_bitsPerSample = ds::bits_per_sample(input.format.sample_format);
    ep_.vi_bytesPerSample = ds::bytes_per_sample(input.format.sample_format);
    ep_.vi_integer = input.format.sample_format != ds::SampleFormat::Float32;
    ep_.vi_numPlanes = input.format.plane_count;
    ep_.vi_width = input.width;
    ep_.vi_height = input.height;
    ep_.vi_subSamplingH = input.format.subsampling_h;
    ep_.vi_subSamplingW = input.format.subsampling_w;

    read_parameters(values);
    validate_parameters(input);
    configure_planes(values);

#ifndef ENABLE_PAR
    ep_.threads = 1;
#endif

    selectFunctions(static_cast<unsigned>(ftype_), static_cast<unsigned>(opt_), &ep_);

    if (ep_.vi_integer) {
      ep_.multiplier = static_cast<float>(1 << (ep_.vi_bitsPerSample - 8));
      ep_.divisor = 1.0f / ep_.multiplier;
      ep_.peak = (1 << ep_.vi_bitsPerSample) - 1;
    }

    if (ftype_ != 0) {
      ep_.f0beta = 1.0f;
    }

    configure_geometry();
    create_fft_plans();
    initialize_sigma_profile();
    prepare_noise_points();
    resize_thread_storage(ep_.threads * fft_threads_ * 16);
  }

  void request_frames(ds::VideoRequestContext& context) const {
    if (ep_.tbsize == 1) {
      context.request_frame(0, context.output_frame);
    } else {
      const int start = std::max(context.output_frame - ep_.tbsize / 2, 0);
      const int stop = std::min(context.output_frame + ep_.tbsize / 2, input_.num_frames - 1);
      for (int frame = start; frame <= stop; ++frame) {
        context.request_frame(0, frame);
      }
    }

    if (noise_points_.empty() || noise_profile_ready_) {
      return;
    }

    for (const NPInfo& point : noise_points_) {
      for (int z = 0; z < ep_.tbsize; ++z) {
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

    if (ep_.tbsize == 1) {
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

    ep_.sbsize = detail::read_int(values, "sbsize", ep_.sbsize);
    ep_.sosize = detail::read_int(values, "sosize", ep_.sosize);
    ep_.tbsize = detail::read_int(values, "tbsize", ep_.tbsize);
    ep_.tosize = detail::read_int(values, "tosize", ep_.tosize);
    ep_.swin = detail::read_int(values, "swin", ep_.swin);
    ep_.twin = detail::read_int(values, "twin", ep_.twin);
    ep_.sbeta = detail::read_float(values, "sbeta", ep_.sbeta);
    ep_.tbeta = detail::read_float(values, "tbeta", ep_.tbeta);
    ep_.zmean = detail::read_bool(values, "zmean", ep_.zmean);
    ep_.f0beta = detail::read_float(values, "f0beta", ep_.f0beta);
    ep_.threads = detail::read_int(values, "threads", ep_.threads);
    ep_.dither = detail::read_int(values, "dither", ep_.dither);
    fft_threads_ = detail::read_int(values, "fft_threads", fft_threads_);
    if (fft_threads_ < 1) {
      fft_threads_ = 1;
    }

    if (smode_ == 0) {
      ep_.sosize = 0;
    }
    if (tmode_ == 0) {
      ep_.tosize = 0;
    }

    if (fft_threads_ > 1 && fft_.has_threading()) {
      fft_.fftwf_init_threads();
      fft_.fftwf_plan_with_nthreads(fft_threads_);
    }

    nlocation_ = detail::read_int_array(values, "nlocation");
    alpha_ = detail::read_float(values, "alpha", ftype_ == 0 ? 5.0f : 7.0f);
    slocation_ = detail::read_float_array(values, "slocation");
    ssx_ = detail::read_float_array(values, "ssx");
    ssy_ = detail::read_float_array(values, "ssy");
    sst_ = detail::read_float_array(values, "sst");
    ssystem_ = detail::read_int(values, "ssystem", 0);

    if (ep_.threads <= 0) {
      ep_.threads = 4;
    }
    if (ep_.threads > 16) {
      ep_.threads = 16;
    }
  }

  void configure_planes(const ds::ParamValues& values) {
    std::fill(std::begin(ep_.process), std::end(ep_.process), 2);

    if (detail::has_param(values, "planes")) {
      const auto planes = detail::read_int_array(values, "planes");
      for (const int plane : planes) {
        if (plane < 0 || plane >= ep_.vi_numPlanes) {
          throw std::runtime_error("plane index out of range");
        }
        ep_.process[plane] = 3;
      }
      return;
    }

    ep_.process[0] = 3;
    ep_.process[1] = 3;
    ep_.process[2] = 3;
    ep_.process[3] = 2;
    ep_.process[0] = detail::read_int(values, "y", ep_.process[0]);
    ep_.process[1] = detail::read_int(values, "u", ep_.process[1]);
    ep_.process[2] = detail::read_int(values, "v", ep_.process[2]);
    ep_.process[3] = detail::read_int(values, "a", ep_.process[3]);
  }

  void validate_parameters(const ds::VideoInputInfo& input) const {
    if (input.width <= 0 || input.height <= 0) {
      throw std::runtime_error("only constant format input supported");
    }

    if (
      (ep_.vi_integer && ep_.vi_bitsPerSample > 16) ||
      (!ep_.vi_integer && ep_.vi_bitsPerSample != 32)
    ) {
      throw std::runtime_error("only 8-16 bit integer and 32 bit float input supported");
    }

    if (ftype_ < 0 || ftype_ > 4) {
      throw std::runtime_error("ftype must be 0, 1, 2, 3, or 4");
    }
    if (ep_.sbsize < 1) {
      throw std::runtime_error("sbsize must be greater than or equal to 1");
    }
    if (smode_ < 0 || smode_ > 1) {
      throw std::runtime_error("smode must be 0 or 1");
    }
    if (smode_ == 0 && !(ep_.sbsize & 1)) {
      throw std::runtime_error("sbsize must be odd when using smode=0");
    }
    if (ep_.sosize < 0 || ep_.sosize >= ep_.sbsize) {
      throw std::runtime_error("sosize must be between 0 and sbsize-1 (inclusive)");
    }
    if (ep_.sosize > ep_.sbsize / 2 && ep_.sbsize % (ep_.sbsize - ep_.sosize) != 0) {
      throw std::runtime_error(
        "spatial overlap greater than 50% requires that sbsize-sosize is a divisor of sbsize"
      );
    }
    if (ep_.tbsize < 1 || ep_.tbsize > 15) {
      throw std::runtime_error("tbsize must be between 1 and 15 (inclusive)");
    }
    if (tmode_ != 0) {
      throw std::runtime_error("tmode must be 0. tmode=1 is not implemented");
    }
    if (tmode_ == 0 && !(ep_.tbsize & 1)) {
      throw std::runtime_error("tbsize must be odd when using tmode=0");
    }
    if (ep_.tosize < 0 || ep_.tosize >= ep_.tbsize) {
      throw std::runtime_error("tosize must be between 0 and tbsize-1 (inclusive)");
    }
    if (ep_.tosize > ep_.tbsize / 2 && ep_.tbsize % (ep_.tbsize - ep_.tosize) != 0) {
      throw std::runtime_error(
        "temporal overlap greater than 50% requires that tbsize-tosize is a divisor of tbsize"
      );
    }
    if (ep_.tbsize > input.num_frames) {
      throw std::runtime_error("tbsize must be less than or equal to the number of frames in the clip");
    }
    if (ep_.swin < 0 || ep_.swin > 11) {
      throw std::runtime_error("swin must be between 0 and 11 (inclusive)");
    }
    if (ep_.twin < 0 || ep_.twin > 11) {
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
    ep_.barea = ep_.sbsize * ep_.sbsize;
    ep_.bvolume = ep_.barea * ep_.tbsize;
    ep_.ccnt = (ep_.sbsize / 2 + 1) * ep_.sbsize * ep_.tbsize;
    ep_.ccnt2 = ep_.ccnt * 2;
    ep_.type = tmode_ * 4 + (ep_.tbsize > 1 ? 2 : 0) + smode_;
    ep_.sbd1 = ep_.sbsize / 2;
    ep_.uf0b = std::abs(ep_.f0beta - 1.0f) >= 0.00005f;
    ep_.inc = (ep_.type & 1) ? ep_.sbsize - ep_.sosize : 1;

    for (int plane = 0; plane < ep_.vi_numPlanes; plane++) {
      const int shift_w = (plane == 1 || plane == 2) ? ep_.vi_subSamplingW : 0;
      const int shift_h = (plane == 1 || plane == 2) ? ep_.vi_subSamplingH : 0;
      ep_.planeWidth[plane] = ep_.vi_width >> shift_w;
      ep_.planeHeight[plane] = ep_.vi_height >> shift_h;
      const int width = ep_.planeWidth[plane];
      const int height = ep_.planeHeight[plane];

      if (smode_ == 0) {
        const int ae = (ep_.sbsize >> 1) << 1;
        ep_.padWidth[plane] = width + ae;
        ep_.padHeight[plane] = height + ae;
        ep_.eHeight[plane] = height;
      } else {
        const int ae = std::max(ep_.sbsize - ep_.sosize, ep_.sosize) * 2;
        ep_.padWidth[plane] = width + EXTRA(width, ep_.sbsize) + ae;
        ep_.padHeight[plane] = height + EXTRA(height, ep_.sbsize) + ae;
        ep_.eHeight[plane] =
          (ep_.padHeight[plane] - ep_.sosize) / (ep_.sbsize - ep_.sosize) *
          (ep_.sbsize - ep_.sosize);
      }

      ep_.padStride[plane] =
        ((ep_.padWidth[plane] * ep_.vi_bytesPerSample - 1) | (FRAME_ALIGN - 1)) + 1;
      ep_.padBlockSize[plane] = ep_.padStride[plane] * ep_.padHeight[plane];
      ep_.eStride[plane] = ((ep_.padWidth[plane] * static_cast<int>(sizeof(float)) - 1) | (FRAME_ALIGN - 1)) + 1;
      ep_.eBatchSize[plane] =
        ((ep_.eHeight[plane] - 1) / ep_.threads / ep_.inc + 1) * ep_.inc;
    }
  }

  void create_fft_plans() {
    ep_.hw = static_cast<float*>(_aligned_malloc((ep_.bvolume + 7) * sizeof(float), FRAME_ALIGN));
    if (!ep_.hw) {
      throw std::runtime_error("malloc failure (hw)");
    }
    createWindow(ep_.hw, tmode_, smode_, &ep_);

    float* dftgr = static_cast<float*>(_aligned_malloc((ep_.bvolume + 7) * sizeof(float), FRAME_ALIGN));
    detail::AlignedPtr<float> dftgr_smart(dftgr);
    ep_.dftgc = static_cast<fftwf_complex*>(
      _aligned_malloc((ep_.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN)
    );
    if (!dftgr || !ep_.dftgc) {
      throw std::runtime_error("malloc failure (dftgr/dftgc)");
    }

    {
      ds::HostGlobalLockGuard fftw_lock("fftw", host_locks_);
      if (ep_.tbsize > 1) {
        ep_.ft = fft_.fftwf_plan_dft_r2c_3d(
          ep_.tbsize,
          ep_.sbsize,
          ep_.sbsize,
          dftgr,
          ep_.dftgc,
          FFTW_PATIENT | FFTW_DESTROY_INPUT
        );
        ep_.fti = fft_.fftwf_plan_dft_c2r_3d(
          ep_.tbsize,
          ep_.sbsize,
          ep_.sbsize,
          ep_.dftgc,
          dftgr,
          FFTW_PATIENT | FFTW_DESTROY_INPUT
        );
      } else {
        ep_.ft = fft_.fftwf_plan_dft_r2c_2d(
          ep_.sbsize,
          ep_.sbsize,
          dftgr,
          ep_.dftgc,
          FFTW_PATIENT | FFTW_DESTROY_INPUT
        );
        ep_.fti = fft_.fftwf_plan_dft_c2r_2d(
          ep_.sbsize,
          ep_.sbsize,
          ep_.dftgc,
          dftgr,
          FFTW_PATIENT | FFTW_DESTROY_INPUT
        );
      }
    }

    float wscale = 0.0f;
    const float* hw_t = ep_.hw;
    float* dftgr_t = dftgr;
    for (int s = 0; s < ep_.tbsize; s++) {
      for (int y = 0; y < ep_.sbsize; y++) {
        for (int x = 0; x < ep_.sbsize; x++) {
          dftgr_t[x] = 255.0f * hw_t[x];
          wscale += hw_t[x] * hw_t[x];
        }
        hw_t += ep_.sbsize;
        dftgr_t += ep_.sbsize;
      }
    }

    fft_.fftwf_execute_dft_r2c(ep_.ft, dftgr, ep_.dftgc);
    wscale_ = 1.0f / wscale;
  }

  void initialize_sigma_profile() {
    const float wscalef = (ftype_ < 2) ? wscale_ : 1.0f;

    ep_.sigmas = static_cast<float*>(_aligned_malloc((ep_.ccnt2 + 7) * sizeof(float), FRAME_ALIGN));
    ep_.sigmas2 = static_cast<float*>(_aligned_malloc((ep_.ccnt2 + 7) * sizeof(float), FRAME_ALIGN));
    ep_.pmins = static_cast<float*>(_aligned_malloc((ep_.ccnt2 + 7) * sizeof(float), FRAME_ALIGN));
    ep_.pmaxs = static_cast<float*>(_aligned_malloc((ep_.ccnt2 + 7) * sizeof(float), FRAME_ALIGN));
    if (!ep_.sigmas || !ep_.sigmas2 || !ep_.pmins || !ep_.pmaxs) {
      throw std::runtime_error("malloc failure (sigmas/sigmas2/pmins/pmaxs)");
    }

    if (!slocation_.empty() || !ssx_.empty() || !ssy_.empty() || !sst_.empty()) {
      initialize_spatially_varying_sigmas(wscalef);
    } else {
      for (int i = 0; i < ep_.ccnt2; i++) {
        ep_.sigmas[i] = sigma_ / wscalef;
      }
    }

    for (int i = 0; i < ep_.ccnt2; i++) {
      ep_.sigmas2[i] = sigma2_ / wscalef;
      ep_.pmins[i] = pmin_ / wscale_;
      ep_.pmaxs[i] = pmax_ / wscale_;
    }
  }

  void initialize_spatially_varying_sigmas(float wscalef) {
    int ndim = 3;
    if (ep_.tbsize == 1) {
      ndim -= 1;
    }
    if (ep_.sbsize == 1) {
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

    const int cpx = ep_.sbsize / 2 + 1;
    float pft = 0.0f;
    float pfy = 0.0f;
    float pfx = 0.0f;

    for (int z = 0; z < ep_.tbsize; z++) {
      const float tval = getSVal(z, ep_.tbsize, tdata, tcnt, pft);
      for (int y = 0; y < ep_.sbsize; y++) {
        const float syval = getSVal(y, ep_.sbsize, sydata, sycnt, pfy);
        for (int x = 0; x < cpx; x++) {
          const float sxval = getSVal(x, ep_.sbsize, sxdata, sxcnt, pfx);
          float val = 0.0f;

          if (ssystem_) {
            const float dw = std::sqrt((pft * pft + pfy * pfy + pfx * pfx) / static_cast<float>(ndim));
            val = interp(dw, tdata, tcnt);
          } else {
            val = tval * syval * sxval;
          }

          const int pos = ((z * ep_.sbsize + y) * cpx + x) * 2;
          ep_.sigmas[pos] = ep_.sigmas[pos + 1] = val / wscalef;
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

      if (fn < 0 || fn > input_.num_frames - ep_.tbsize) {
        throw std::runtime_error("invalid frame number in nlocation (" + std::to_string(fn) + ")");
      }
      if (plane < 0 || plane >= ep_.vi_numPlanes) {
        throw std::runtime_error("invalid plane number in nlocation (" + std::to_string(plane) + ")");
      }

      const int height = ep_.vi_height >> (plane > 0 ? ep_.vi_subSamplingH : 0);
      if (y < 0 || y > height - ep_.sbsize) {
        throw std::runtime_error("invalid y pos in nlocation (" + std::to_string(y) + ")");
      }

      const int width = ep_.vi_width >> (plane > 0 ? ep_.vi_subSamplingW : 0);
      if (x < 0 || x > width - ep_.sbsize) {
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

    std::memset(ep_.sigmas, 0, static_cast<std::size_t>(ep_.ccnt2) * sizeof(float));

    float* hw2 = static_cast<float*>(_aligned_malloc((ep_.bvolume + 7) * sizeof(float), FRAME_ALIGN));
    if (!hw2) {
      throw std::runtime_error("malloc failure (hw2)");
    }
    detail::AlignedPtr<float> hw2_smart(hw2);
    createWindow(hw2, 0, 0, &ep_);

    float* dftr = static_cast<float*>(_aligned_malloc((ep_.bvolume + 7) * sizeof(float), FRAME_ALIGN));
    auto* dftgc2 = static_cast<fftwf_complex*>(
      _aligned_malloc((ep_.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN)
    );
    if (!dftr || !dftgc2) {
      throw std::runtime_error("malloc failure (dftr/dftgc2)");
    }
    detail::AlignedPtr<float> dftr_smart(dftr);
    detail::AlignedPtr<fftwf_complex> dftgc2_smart(dftgc2);

    float wscale2 = 0.0f;
    int w = 0;
    for (int s = 0; s < ep_.tbsize; s++) {
      for (int y = 0; y < ep_.sbsize; y++) {
        for (int x = 0; x < ep_.sbsize; x++, w++) {
          dftr[w] = 255.0f * hw2[w];
          wscale2 += hw2[w] * hw2[w];
        }
      }
    }
    wscale2 = 1.0f / wscale2;
    fft_.fftwf_execute_dft_r2c(ep_.ft, dftr, dftgc2);

    auto* dftc = static_cast<fftwf_complex*>(
      _aligned_malloc((ep_.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN)
    );
    auto* dftc2 = static_cast<fftwf_complex*>(
      _aligned_malloc((ep_.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN)
    );
    if (!dftc || !dftc2) {
      throw std::runtime_error("malloc failure (dftc/dftc2)");
    }
    detail::AlignedPtr<fftwf_complex> dftc_smart(dftc);
    detail::AlignedPtr<fftwf_complex> dftc2_smart(dftc2);

    for (const NPInfo& point : noise_points_) {
      for (int z = 0; z < ep_.tbsize; z++) {
        auto frame = provider.get(0, point.fn + z);
        if (!frame.has_value()) {
          throw std::runtime_error(frame.error().message);
        }

        const ds::PlaneView& plane = frame.value().frame.plane(point.b);
        const int stride_elements = static_cast<int>(plane.stride_bytes) / ep_.vi_bytesPerSample;

        if (ep_.vi_bytesPerSample == 1) {
          const auto* srcp =
            static_cast<const std::uint8_t*>(plane.data) + stride_elements * point.y + point.x;
          proc0_c(srcp, hw2 + ep_.barea * z, dftr + ep_.barea * z, stride_elements, ep_.sbsize, ep_.divisor);
        } else if (ep_.vi_bytesPerSample == 2) {
          const auto* srcp =
            static_cast<const std::uint16_t*>(plane.data) + stride_elements * point.y + point.x;
          proc0_c(srcp, hw2 + ep_.barea * z, dftr + ep_.barea * z, stride_elements, ep_.sbsize, ep_.divisor);
        } else {
          const auto* srcp =
            static_cast<const float*>(plane.data) + stride_elements * point.y + point.x;
          proc0_c(srcp, hw2 + ep_.barea * z, dftr + ep_.barea * z, stride_elements, ep_.sbsize, ep_.divisor);
        }
      }

      fft_.fftwf_execute_dft_r2c(ep_.ft, dftr, dftc);

      if (ep_.zmean) {
        removeMean_c(
          reinterpret_cast<float*>(dftc),
          reinterpret_cast<const float*>(dftgc2),
          ep_.ccnt2,
          reinterpret_cast<float*>(dftc2)
        );
      }

      for (int h = 0; h < ep_.ccnt2; h += 2) {
        const float* dftc_f = reinterpret_cast<float*>(dftc);
        const float psd = dftc_f[h] * dftc_f[h] + dftc_f[h + 1] * dftc_f[h + 1];
        ep_.sigmas[h] += psd;
        ep_.sigmas[h + 1] += psd;
      }
    }

    const float scale = 1.0f / static_cast<float>(noise_points_.size());
    for (int h = 0; h < ep_.ccnt2; h++) {
      ep_.sigmas[h] *= scale * (wscale2 / wscale_) * alpha_;
    }

    noise_profile_ready_ = true;
  }

  void process_spatial_frame(
    unsigned int thread_id,
    const ds::VideoFrameView& src,
    ds::MutableVideoFrameView& dst
  ) {
    for (int plane = 0; plane < ep_.vi_numPlanes; plane++) {
      const int height = ep_.planeHeight[plane];
      const int width = ep_.planeWidth[plane];
      const auto& src_plane = src.plane(plane);
      const auto& dst_plane = dst.plane(plane);
      const int src_stride = detail::plane_stride(src_plane);
      const int dst_stride = detail::plane_stride(dst_plane);
      const auto* src_ptr = detail::readable_plane_data(src, plane);
      auto* dst_ptr = detail::writable_plane_data(dst, plane);

      if (ep_.process[plane] == 3) {
        auto* pad = static_cast<unsigned char*>(_aligned_malloc(ep_.padBlockSize[plane] * ep_.tbsize, FRAME_ALIGN));
        if (!pad) {
          throw std::runtime_error("pad0 allocation failed.");
        }
        detail::AlignedPtr<unsigned char> pad0_smart(pad);
        ep_.copyPad(plane, src_ptr, src_stride, pad, &ep_);
        ep_.func_0(thread_id, plane, pad, dst_ptr, dst_stride, &ep_);
      } else if (ep_.process[plane] == 2) {
        detail::framecpy(
          dst_ptr,
          dst_stride,
          src_ptr,
          src_stride,
          width * ep_.vi_bytesPerSample,
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
    const int pos = ep_.tbsize / 2;

    for (int plane = 0; plane < ep_.vi_numPlanes; plane++) {
      const int height = ep_.planeHeight[plane];
      const int width = ep_.planeWidth[plane];
      const auto& src0_plane = current.plane(plane);
      const auto& dst_plane = dst.plane(plane);
      const int src0_stride = detail::plane_stride(src0_plane);
      const int dst_stride = detail::plane_stride(dst_plane);
      const auto* src0_ptr = detail::readable_plane_data(current, plane);
      auto* dst_ptr = detail::writable_plane_data(dst, plane);

      if (ep_.process[plane] == 3) {
        auto* pad0 = static_cast<unsigned char*>(_aligned_malloc(ep_.padBlockSize[plane] * ep_.tbsize, FRAME_ALIGN));
        if (!pad0) {
          throw std::runtime_error("pad0 allocation failed.");
        }
        detail::AlignedPtr<unsigned char> pad0_smart(pad0);

        for (int i = 0; i < ep_.tbsize; i++) {
          const int fn = i + n - pos;
          const int fn_real = std::min(std::max(fn, 0), input_.num_frames - 1);
          auto src_frame = provider.get(0, fn_real);
          if (!src_frame.has_value()) {
            throw std::runtime_error(src_frame.error().message);
          }

          const auto& src_plane = src_frame.value().frame.plane(plane);
          const int src_stride = detail::plane_stride(src_plane);
          const auto* src_ptr = detail::readable_plane_data(src_frame.value().frame, plane);
          auto* pad = pad0 + ep_.padBlockSize[plane] * i;
          ep_.copyPad(plane, src_ptr, src_stride, pad, &ep_);
        }

        ep_.func_1(thread_id, plane, pad0, dst_ptr, dst_stride, pos, &ep_);
      } else if (ep_.process[plane] == 2) {
        detail::framecpy(
          dst_ptr,
          dst_stride,
          src0_ptr,
          src0_stride,
          width * ep_.vi_bytesPerSample,
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
    ep_.ebuff.resize(static_cast<std::size_t>(count), nullptr);
    ep_.dftr.resize(static_cast<std::size_t>(count), nullptr);
    ep_.dftc.resize(static_cast<std::size_t>(count), nullptr);
    ep_.dftc2.resize(static_cast<std::size_t>(count), nullptr);
    ep_.rngs.resize(static_cast<std::size_t>(count));
    ep_.d_buffs.resize(static_cast<std::size_t>(count), nullptr);
  }

  void ensure_thread_buffers_unlocked(unsigned int thread_id) {
    if (ep_.ebuff[thread_id]) {
      return;
    }

    ep_.ebuff[thread_id] = static_cast<float*>(
      _aligned_malloc(sizeof(float) * ep_.eStride[0] * ep_.padHeight[0], FRAME_ALIGN)
    );
    ep_.dftr[thread_id] = static_cast<float*>(
      _aligned_malloc(sizeof(float) * (((ep_.bvolume + 7) | 15) + 1) * ep_.threads, FRAME_ALIGN)
    );
    ep_.dftc[thread_id] = static_cast<fftwf_complex*>(
      _aligned_malloc(sizeof(fftwf_complex) * (((ep_.ccnt + 7) | 15) + 1) * ep_.threads, FRAME_ALIGN)
    );
    ep_.dftc2[thread_id] = static_cast<fftwf_complex*>(
      _aligned_malloc(sizeof(fftwf_complex) * (((ep_.ccnt + 7) | 15) + 1) * ep_.threads, FRAME_ALIGN)
    );

    if (!ep_.ebuff[thread_id] || !ep_.dftr[thread_id] || !ep_.dftc[thread_id] || !ep_.dftc2[thread_id]) {
      throw std::runtime_error("thread buffer allocation failed.");
    }

    if (ep_.dither > 0) {
      ep_.d_buffs[thread_id] = static_cast<float*>(
        _aligned_malloc(sizeof(float) * 2 * ep_.vi_width, FRAME_ALIGN)
      );
      if (!ep_.d_buffs[thread_id]) {
        throw std::runtime_error("dither buffer allocation failed.");
      }
      if (ep_.dither > 1) {
        ep_.rngs[thread_id] = std::make_unique<std::mt19937>(std::random_device{}());
      }
    }
  }

  ds::VideoInputInfo input_{};
  ds::HostGlobalLockCallbacks host_locks_{};
  DFTTestData ep_{};
  FFTFunctionPointers fft_{};
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
