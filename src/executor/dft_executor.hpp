#pragma once

#include "dft_common.h"
#include "executor/dft_memory_plan.hpp"
#include "workspace/dft_workspace.hpp"

#include <array>
#include <memory>

namespace neo_dfttest {

struct DftExecutorCapabilities {
  fft::MemoryDomain working_memory_domain {fft::MemoryDomain::host};
  bool supports_device_coefficients {false};
  bool supports_asynchronous_batches {false};
  bool supports_host_copy_pad {true};
  bool uses_host_fallback {false};
};

struct DftCopyPadRequest {
  int plane {0};
  DFTPlaneBytes source;
  DFTMutablePlaneBytes destination;
  const DFTKernelContext& context;
};

struct DftProcessSpatialRequest {
  DftWorkspaceLease workspace;
  int plane {0};
  DFTPlaneBytes source;
  DFTMutablePlaneBytes destination;
  const DFTKernelContext& context;
};

struct DftProcessTemporalRequest {
  DftWorkspaceLease workspace;
  int plane {0};
  DFTPlaneBytes source;
  DFTMutablePlaneBytes destination;
  int temporal_position {0};
  const DFTKernelContext& context;
};

enum class DftProcessMode {
  spatial,
  temporal,
};

struct DftProcessRequest {
  DftWorkspaceLease workspace;
  int plane {0};
  DftProcessMode mode {DftProcessMode::spatial};
  std::array<DFTPlaneBytes, kMaxDftTemporalFrames> sources {};
  int source_count {0};
  DFTMutablePlaneBytes destination;
  int temporal_position {0};
  const DFTKernelContext& context;
};

struct DftFramePlaneRequest {
  std::array<DFTPlaneBytes, kMaxDftTemporalFrames> sources {};
  int source_count {0};
  DFTMutablePlaneBytes destination;
};

struct DftFrameProcessRequest {
  DftWorkspaceLease workspace;
  DftProcessMode mode {DftProcessMode::spatial};
  std::array<DftFramePlaneRequest, 4> planes {};
  int plane_count {0};
  int temporal_position {0};
  const DFTKernelContext& context;
};

class DftExecutor {
public:
  virtual ~DftExecutor() = default;

  [[nodiscard]] virtual DftExecutorCapabilities capabilities() const noexcept = 0;
  [[nodiscard]] virtual DftMemoryPlan memory_plan() const noexcept = 0;
  [[nodiscard]] virtual DFTBatchPolicy make_batch_policy(
    const DFTBlockSettings& block,
    const fft::BackendCapabilities& fft_capabilities
  ) const noexcept = 0;

  virtual void copy_pad(DftCopyPadRequest request) = 0;
  virtual void process_frame(DftFrameProcessRequest request) = 0;
  virtual void process(DftProcessRequest request) = 0;
  virtual void process_spatial(DftProcessSpatialRequest request) = 0;
  virtual void process_temporal(DftProcessTemporalRequest request) = 0;
};

std::unique_ptr<DftExecutor> create_cpu_dft_executor(unsigned opt, DFTClipFormat format);
#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)
std::unique_ptr<DftExecutor> create_vulkan_dft_executor(
  unsigned opt,
  DFTClipFormat format,
  fft::VulkanRuntime* runtime
);
#endif

} // namespace neo_dfttest
