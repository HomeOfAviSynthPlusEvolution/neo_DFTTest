#include "executor/dft_executor.hpp"

namespace neo_dfttest {

class VulkanDftExecutor final : public DftExecutor {
public:
  VulkanDftExecutor(unsigned opt, DFTClipFormat format)
    : fallback_(create_cpu_dft_executor(opt, format)) {}

  [[nodiscard]] DftExecutorCapabilities capabilities() const noexcept override {
    return DftExecutorCapabilities{
      fft::MemoryDomain::device,
      true,
      true,
      true,
      true
    };
  }

  [[nodiscard]] DftMemoryPlan memory_plan() const noexcept override {
    return dft_vulkan_memory_plan(true);
  }

  [[nodiscard]] DFTBatchPolicy make_batch_policy(
    const DFTBlockSettings& block,
    const fft::BackendCapabilities& fft_capabilities
  ) const noexcept override {
    return fallback_->make_batch_policy(block, fft_capabilities);
  }

  void copy_pad(DftCopyPadRequest request) override {
    fallback_->copy_pad(request);
  }

  void process_frame(DftFrameProcessRequest request) override {
    fallback_->process_frame(request);
  }

  void process(DftProcessRequest request) override {
    fallback_->process(request);
  }

  void process_spatial(DftProcessSpatialRequest request) override {
    fallback_->process_spatial(request);
  }

  void process_temporal(DftProcessTemporalRequest request) override {
    fallback_->process_temporal(request);
  }

private:
  std::unique_ptr<DftExecutor> fallback_;
};

std::unique_ptr<DftExecutor> create_vulkan_dft_executor(unsigned opt, DFTClipFormat format) {
  return std::make_unique<VulkanDftExecutor>(opt, format);
}

} // namespace neo_dfttest
