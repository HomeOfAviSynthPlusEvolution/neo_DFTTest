#include "fft/fft_backend.hpp"

#include "fftwlite.h"

#include <limits>
#include <stdexcept>

namespace neo_dfttest::fft {
namespace {

class FftwPlan final : public Plan::State {
public:
  FftwPlan(FFTFunctionPointers& api, fftwf_plan plan) noexcept
    : api_(api), plan_(plan) {}

  FftwPlan(const FftwPlan&) = delete;
  FftwPlan& operator=(const FftwPlan&) = delete;

  ~FftwPlan() override {
    if (plan_) {
      api_.fftwf_destroy_plan(plan_);
    }
  }

  [[nodiscard]] fftwf_plan get() const noexcept {
    return plan_;
  }

private:
  FFTFunctionPointers& api_;
  fftwf_plan plan_ {nullptr};
};

fftwf_complex* fftw_complex(Complex* data) noexcept {
  return reinterpret_cast<fftwf_complex*>(data);
}

unsigned fftw_flags(PlanOptions options) noexcept {
  unsigned flags = 0;
  if (options.patient) {
    flags |= FFTW_PATIENT;
  }
  if (options.destroy_input) {
    flags |= FFTW_DESTROY_INPUT;
  }
  return flags;
}

const FftwPlan& fftw_plan(const Plan& plan) noexcept {
  return static_cast<const FftwPlan&>(plan.state());
}

class FftwBackend final : public Backend {
public:
  FftwBackend() = default;
  FftwBackend(const FftwBackend&) = delete;
  FftwBackend& operator=(const FftwBackend&) = delete;

  ~FftwBackend() override {
    unload();
  }

  void load() override {
    try {
      api_.load();
    } catch (const char* error) {
      throw std::runtime_error(error);
    }
  }

  void unload() noexcept override {
    api_.free();
  }

  bool loaded() const noexcept override {
    return api_.library != nullptr;
  }

  bool has_threading() const noexcept override {
    return api_.has_threading();
  }

  void set_thread_count(int nthreads) override {
    api_.fftwf_init_threads();
    api_.fftwf_plan_with_nthreads(nthreads);
  }

  BackendCapabilities capabilities() const noexcept override {
    return BackendCapabilities{
      true,
      false,
      false,
      1,
      std::numeric_limits<int>::max(),
    };
  }

  Plan make_plan(
    TransformDirection direction,
    TransformShape shape,
    Real* real_workspace,
    Complex* complex_workspace,
    PlanOptions options
  ) override {
    const unsigned flags = fftw_flags(options);
    fftwf_plan plan = nullptr;
    if (shape.rank == 2) {
      if (direction == TransformDirection::r2c) {
        plan = api_.fftwf_plan_dft_r2c_2d(shape.extents[0], shape.extents[1], real_workspace, fftw_complex(complex_workspace), flags);
      } else {
        plan = api_.fftwf_plan_dft_c2r_2d(shape.extents[0], shape.extents[1], fftw_complex(complex_workspace), real_workspace, flags);
      }
    } else if (shape.rank == 3) {
      if (direction == TransformDirection::r2c) {
        plan = api_.fftwf_plan_dft_r2c_3d(shape.extents[0], shape.extents[1], shape.extents[2], real_workspace, fftw_complex(complex_workspace), flags);
      } else {
        plan = api_.fftwf_plan_dft_c2r_3d(shape.extents[0], shape.extents[1], shape.extents[2], fftw_complex(complex_workspace), real_workspace, flags);
      }
    } else {
      throw std::runtime_error("unsupported FFT rank");
    }

    if (!plan) {
      throw std::runtime_error("failed to create FFTW plan");
    }
    return Plan(std::make_unique<FftwPlan>(api_, plan));
  }

  Completion submit_r2c(const Plan& plan, R2CBatch batch, SubmitOptions) const noexcept override {
    const fftwf_plan native_plan = fftw_plan(plan).get();
    auto* in = batch.input.data;
    auto* out = batch.output.data;
    for (int index = 0; index < batch.count; ++index) {
      api_.fftwf_execute_dft_r2c(native_plan, in, fftw_complex(out));
      in += batch.input.stride_elements;
      out += batch.output.stride_elements;
    }
    return Completion::completed();
  }

  Completion submit_c2r(const Plan& plan, C2RBatch batch, SubmitOptions) const noexcept override {
    const fftwf_plan native_plan = fftw_plan(plan).get();
    auto* in = batch.input.data;
    auto* out = batch.output.data;
    for (int index = 0; index < batch.count; ++index) {
      api_.fftwf_execute_dft_c2r(native_plan, fftw_complex(in), out);
      in += batch.input.stride_elements;
      out += batch.output.stride_elements;
    }
    return Completion::completed();
  }

private:
  FFTFunctionPointers api_{};
};

} // namespace

std::unique_ptr<Backend> create_fftw_backend() {
  return std::make_unique<FftwBackend>();
}

} // namespace neo_dfttest::fft
