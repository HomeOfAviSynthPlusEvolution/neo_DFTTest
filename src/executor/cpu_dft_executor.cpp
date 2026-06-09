#include "executor/dft_executor.hpp"

#include "core.h"

#include <stdexcept>

namespace neo_dfttest {

class CpuDftExecutor final : public DftExecutor {
public:
  explicit CpuDftExecutor(DFTCpuProcessDispatch dispatch) noexcept
    : process_spatial_(dispatch.process_spatial),
      process_temporal_(dispatch.process_temporal) {}

  void process_spatial(
    unsigned int thread_id,
    int plane,
    DFTPlaneBytes src,
    DFTMutablePlaneBytes dst,
    const DFTKernelContext& context
  ) override {
    if (!process_spatial_) {
      throw std::runtime_error("CPU DFT executor has no spatial processor");
    }
    process_spatial_(thread_id, plane, src, dst, context);
  }

  void process_temporal(
    unsigned int thread_id,
    int plane,
    DFTPlaneBytes src,
    DFTMutablePlaneBytes dst,
    int temporal_position,
    const DFTKernelContext& context
  ) override {
    if (!process_temporal_) {
      throw std::runtime_error("CPU DFT executor has no temporal processor");
    }
    process_temporal_(thread_id, plane, src, dst, temporal_position, context);
  }

private:
  DFTProcessSpatialFunction process_spatial_ {nullptr};
  DFTProcessTemporalFunction process_temporal_ {nullptr};
};

std::unique_ptr<DftExecutor> create_cpu_dft_executor(unsigned opt, DFTClipFormat format) {
  return std::make_unique<CpuDftExecutor>(select_cpu_process_dispatch(opt, format));
}

} // namespace neo_dfttest
