#ifndef NEO_DFTTEST_FFT_BACKEND_HPP
#define NEO_DFTTEST_FFT_BACKEND_HPP

#include <cstddef>
#include <memory>

#include "fftwlite.h"

namespace neo_dfttest::fft {

using Real = float;
using Complex = fftwf_complex;
using Plan = fftwf_plan;

constexpr unsigned kPatientDestroyInputPlanFlags = FFTW_PATIENT | FFTW_DESTROY_INPUT;

enum class MemoryDomain {
  host,
  device,
  unified,
};

enum class ExecutionMode {
  synchronous,
  may_defer,
};

struct SubmitOptions {
  ExecutionMode mode {ExecutionMode::may_defer};
  void* stream {nullptr};
};

struct BackendCapabilities {
  bool batched_transforms {false};
  bool asynchronous_submission {false};
  bool device_memory {false};
  int preferred_batch_size {1};
  int max_batch_size {1};
};

struct RealBatchView {
  Real* data {nullptr};
  std::ptrdiff_t stride_elements {0};
  MemoryDomain domain {MemoryDomain::host};
};

struct ComplexBatchView {
  Complex* data {nullptr};
  std::ptrdiff_t stride_elements {0};
  MemoryDomain domain {MemoryDomain::host};
};

struct R2CBatch {
  RealBatchView input;
  ComplexBatchView output;
  int count {1};
};

struct C2RBatch {
  ComplexBatchView input;
  RealBatchView output;
  int count {1};
};

class Completion {
public:
  Completion(const Completion&) = delete;
  Completion& operator=(const Completion&) = delete;
  Completion(Completion&&) noexcept = default;
  Completion& operator=(Completion&&) noexcept = default;

  static Completion completed() noexcept {
    return Completion();
  }

  void wait() noexcept {}
  [[nodiscard]] bool ready() const noexcept { return true; }

private:
  Completion() = default;
};

class Backend {
public:
  virtual ~Backend() = default;

  virtual void load() = 0;
  virtual void unload() noexcept = 0;
  virtual bool loaded() const noexcept = 0;

  virtual bool has_threading() const noexcept = 0;
  virtual void set_thread_count(int nthreads) = 0;
  virtual BackendCapabilities capabilities() const noexcept = 0;

  virtual Plan plan_r2c_2d(int n0, int n1, Real* in, Complex* out, unsigned flags) = 0;
  virtual Plan plan_c2r_2d(int n0, int n1, Complex* in, Real* out, unsigned flags) = 0;
  virtual Plan plan_r2c_3d(int n0, int n1, int n2, Real* in, Complex* out, unsigned flags) = 0;
  virtual Plan plan_c2r_3d(int n0, int n1, int n2, Complex* in, Real* out, unsigned flags) = 0;

  virtual Completion submit_r2c(Plan plan, R2CBatch batch, SubmitOptions options = {}) const noexcept = 0;
  virtual Completion submit_c2r(Plan plan, C2RBatch batch, SubmitOptions options = {}) const noexcept = 0;
  virtual void destroy_plan(Plan plan) noexcept = 0;
};

inline R2CBatch single_r2c_batch(Real* input, Complex* output, std::ptrdiff_t real_stride, std::ptrdiff_t complex_stride) noexcept {
  return R2CBatch{
    RealBatchView{input, real_stride, MemoryDomain::host},
    ComplexBatchView{output, complex_stride, MemoryDomain::host},
    1,
  };
}

inline C2RBatch single_c2r_batch(Complex* input, Real* output, std::ptrdiff_t complex_stride, std::ptrdiff_t real_stride) noexcept {
  return C2RBatch{
    ComplexBatchView{input, complex_stride, MemoryDomain::host},
    RealBatchView{output, real_stride, MemoryDomain::host},
    1,
  };
}

std::unique_ptr<Backend> create_fftw_backend();

} // namespace neo_dfttest::fft

#endif
