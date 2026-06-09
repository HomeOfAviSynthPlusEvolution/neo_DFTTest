#include "fft/fft_backend.hpp"

#include <stdexcept>

namespace neo_dfttest::fft {
namespace {

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

  Plan plan_r2c_2d(int n0, int n1, Real* in, Complex* out, unsigned flags) override {
    return api_.fftwf_plan_dft_r2c_2d(n0, n1, in, out, flags);
  }

  Plan plan_c2r_2d(int n0, int n1, Complex* in, Real* out, unsigned flags) override {
    return api_.fftwf_plan_dft_c2r_2d(n0, n1, in, out, flags);
  }

  Plan plan_r2c_3d(int n0, int n1, int n2, Real* in, Complex* out, unsigned flags) override {
    return api_.fftwf_plan_dft_r2c_3d(n0, n1, n2, in, out, flags);
  }

  Plan plan_c2r_3d(int n0, int n1, int n2, Complex* in, Real* out, unsigned flags) override {
    return api_.fftwf_plan_dft_c2r_3d(n0, n1, n2, in, out, flags);
  }

  void execute_r2c(Plan plan, Real* in, Complex* out) const noexcept override {
    api_.fftwf_execute_dft_r2c(plan, in, out);
  }

  void execute_c2r(Plan plan, Complex* in, Real* out) const noexcept override {
    api_.fftwf_execute_dft_c2r(plan, in, out);
  }

  void destroy_plan(Plan plan) noexcept override {
    if (plan) {
      api_.fftwf_destroy_plan(plan);
    }
  }

private:
  FFTFunctionPointers api_{};
};

} // namespace

std::unique_ptr<Backend> create_fftw_backend() {
  return std::make_unique<FftwBackend>();
}

} // namespace neo_dfttest::fft
