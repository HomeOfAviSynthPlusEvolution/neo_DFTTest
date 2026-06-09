#include "fft/fft_backend.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "pocketfft_hdronly.h"

namespace neo_dfttest::fft {
namespace {

class PocketfftPlan final : public Plan::State {
public:
  PocketfftPlan(TransformDirection direction, TransformShape shape)
    : direction_(direction),
      shape_(make_shape(shape)),
      complex_shape_(shape_),
      axes_(make_axes(shape_)),
      real_strides_(make_strides<Real>(shape_)) {
    if (shape.rank < 1 || shape.rank > 3) {
      throw std::runtime_error("unsupported pocketfft rank");
    }
    complex_shape_.back() = complex_shape_.back() / 2 + 1;
    complex_strides_ = make_strides<Complex>(complex_shape_);
  }

  [[nodiscard]] TransformDirection direction() const noexcept {
    return direction_;
  }

  [[nodiscard]] const pocketfft::shape_t& shape() const noexcept {
    return shape_;
  }

  [[nodiscard]] const pocketfft::shape_t& axes() const noexcept {
    return axes_;
  }

  [[nodiscard]] const pocketfft::stride_t& real_strides() const noexcept {
    return real_strides_;
  }

  [[nodiscard]] const pocketfft::stride_t& complex_strides() const noexcept {
    return complex_strides_;
  }

private:
  static pocketfft::shape_t make_shape(TransformShape shape) {
    if (shape.rank < 1 || shape.rank > 3) {
      throw std::runtime_error("unsupported pocketfft rank");
    }

    pocketfft::shape_t output;
    output.reserve(static_cast<std::size_t>(shape.rank));
    for (int index = 0; index < shape.rank; ++index) {
      if (shape.extents[static_cast<std::size_t>(index)] <= 0) {
        throw std::runtime_error("invalid FFT shape extent");
      }
      output.push_back(static_cast<std::size_t>(shape.extents[static_cast<std::size_t>(index)]));
    }
    return output;
  }

  static pocketfft::shape_t make_axes(const pocketfft::shape_t& shape) {
    pocketfft::shape_t axes(shape.size());
    std::iota(axes.begin(), axes.end(), std::size_t{0});
    return axes;
  }

  template<typename T>
  static pocketfft::stride_t make_strides(const pocketfft::shape_t& shape) {
    pocketfft::stride_t strides(shape.size());
    std::ptrdiff_t stride = static_cast<std::ptrdiff_t>(sizeof(T));
    for (auto index = static_cast<std::ptrdiff_t>(shape.size()) - 1; index >= 0; --index) {
      strides[static_cast<std::size_t>(index)] = stride;
      stride *= static_cast<std::ptrdiff_t>(shape[static_cast<std::size_t>(index)]);
    }
    return strides;
  }

  TransformDirection direction_;
  pocketfft::shape_t shape_;
  pocketfft::shape_t complex_shape_;
  pocketfft::shape_t axes_;
  pocketfft::stride_t real_strides_;
  pocketfft::stride_t complex_strides_;
};

const PocketfftPlan& pocketfft_plan(const Plan& plan) noexcept {
  return static_cast<const PocketfftPlan&>(plan.state());
}

class PocketfftBackend final : public Backend {
public:
  void load() override {
    loaded_ = true;
  }

  void unload() noexcept override {
    loaded_ = false;
  }

  bool loaded() const noexcept override {
    return loaded_;
  }

  bool has_threading() const noexcept override {
    return true;
  }

  void set_thread_count(int nthreads) override {
    thread_count_ = static_cast<std::size_t>(std::max(1, nthreads));
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
    Real*,
    Complex*,
    PlanOptions
  ) override {
    return Plan(std::make_unique<PocketfftPlan>(direction, shape));
  }

  Completion submit_r2c(const Plan& plan, R2CBatch batch, SubmitOptions) const override {
    const auto& native = pocketfft_plan(plan);
    if (native.direction() != TransformDirection::r2c) {
      throw std::runtime_error("pocketfft r2c submitted with a non-r2c plan");
    }

    auto* in = batch.input.data;
    auto* out = batch.output.data;
    for (int index = 0; index < batch.count; ++index) {
      pocketfft::r2c<float>(
        native.shape(),
        native.real_strides(),
        native.complex_strides(),
        native.axes(),
        pocketfft::FORWARD,
        in,
        out,
        1.0f,
        thread_count_
      );
      in += batch.input.stride_elements;
      out += batch.output.stride_elements;
    }
    return Completion::completed();
  }

  Completion submit_c2r(const Plan& plan, C2RBatch batch, SubmitOptions) const override {
    const auto& native = pocketfft_plan(plan);
    if (native.direction() != TransformDirection::c2r) {
      throw std::runtime_error("pocketfft c2r submitted with a non-c2r plan");
    }

    auto* in = batch.input.data;
    auto* out = batch.output.data;
    for (int index = 0; index < batch.count; ++index) {
      pocketfft::c2r<float>(
        native.shape(),
        native.complex_strides(),
        native.real_strides(),
        native.axes(),
        pocketfft::BACKWARD,
        in,
        out,
        1.0f,
        thread_count_
      );
      in += batch.input.stride_elements;
      out += batch.output.stride_elements;
    }
    return Completion::completed();
  }

private:
  bool loaded_ {false};
  std::size_t thread_count_ {1};
};

} // namespace

std::unique_ptr<Backend> create_pocketfft_backend() {
  return std::make_unique<PocketfftBackend>();
}

} // namespace neo_dfttest::fft
