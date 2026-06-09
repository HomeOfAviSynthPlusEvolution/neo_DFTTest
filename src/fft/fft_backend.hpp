#ifndef NEO_DFTTEST_FFT_BACKEND_HPP
#define NEO_DFTTEST_FFT_BACKEND_HPP

#include <array>
#include <complex>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace neo_dfttest::fft {

using Real = float;
using Complex = std::complex<float>;

static_assert(std::is_standard_layout_v<Complex>);
static_assert(sizeof(Complex) == sizeof(float) * 2);
static_assert(alignof(Complex) == alignof(float));

enum class TransformDirection {
  r2c,
  c2r,
};

struct TransformShape {
  int rank {2};
  std::array<int, 3> extents {1, 1, 1};
};

struct PlanOptions {
  bool patient {true};
  bool destroy_input {true};
};

constexpr PlanOptions kPatientDestroyInputPlanOptions {};

struct BatchLayout {
  int max_batch {1};
  std::ptrdiff_t real_stride_elements {0};
  std::ptrdiff_t complex_stride_elements {0};
};

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
  class State {
  public:
    virtual ~State() = default;
    virtual void wait() noexcept = 0;
    [[nodiscard]] virtual bool ready() const noexcept = 0;
  };

  Completion(const Completion&) = delete;
  Completion& operator=(const Completion&) = delete;
  Completion(Completion&&) noexcept = default;
  Completion& operator=(Completion&&) noexcept = default;

  explicit Completion(std::unique_ptr<State> state) noexcept
    : state_(std::move(state)) {}

  static Completion completed() noexcept {
    return Completion();
  }

  void wait() noexcept {
    if (state_) {
      state_->wait();
    }
  }

  [[nodiscard]] bool ready() const noexcept {
    return !state_ || state_->ready();
  }

private:
  Completion() = default;
  std::unique_ptr<State> state_;
};

class Plan {
public:
  class State {
  public:
    virtual ~State() = default;
  };

  Plan() = default;
  Plan(std::unique_ptr<State> state, BatchLayout layout) noexcept
    : state_(std::move(state)), layout_(layout) {}

  Plan(const Plan&) = delete;
  Plan& operator=(const Plan&) = delete;
  Plan(Plan&&) noexcept = default;
  Plan& operator=(Plan&&) noexcept = default;

  [[nodiscard]] bool valid() const noexcept {
    return state_ != nullptr;
  }

  [[nodiscard]] const State& state() const noexcept {
    return *state_;
  }

  [[nodiscard]] BatchLayout layout() const noexcept {
    return layout_;
  }

private:
  std::unique_ptr<State> state_;
  BatchLayout layout_;
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

  virtual Plan make_plan(
    TransformDirection direction,
    TransformShape shape,
    BatchLayout layout,
    Real* real_workspace,
    Complex* complex_workspace,
    PlanOptions options = {}
  ) = 0;

  virtual Completion submit_r2c(const Plan& plan, R2CBatch batch, SubmitOptions options = {}) const = 0;
  virtual Completion submit_c2r(const Plan& plan, C2RBatch batch, SubmitOptions options = {}) const = 0;
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
std::unique_ptr<Backend> create_pocketfft_backend();

} // namespace neo_dfttest::fft

#endif
