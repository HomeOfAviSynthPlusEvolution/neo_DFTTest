#include "executor/dft_executor.hpp"

#include "core.h"

#include <stdexcept>

namespace neo_dfttest {

class CpuDftExecutor final : public DftExecutor {
public:
  CpuDftExecutor(DFTCopyPadFunction copy_pad, DFTCpuProcessDispatch dispatch) noexcept
    : copy_pad_(copy_pad),
      process_spatial_(dispatch.process_spatial),
      process_temporal_(dispatch.process_temporal) {}

  [[nodiscard]] DftExecutorCapabilities capabilities() const noexcept override {
    return DftExecutorCapabilities{
      fft::MemoryDomain::host,
      false,
      false,
      true
    };
  }

  [[nodiscard]] DFTBatchPolicy make_batch_policy(
    const DFTBlockSettings& block,
    const fft::BackendCapabilities& fft_capabilities
  ) const noexcept override {
    return make_cpu_dft_batch_policy(block, fft_capabilities);
  }

  void copy_pad(DftCopyPadRequest request) override {
    if (!copy_pad_) {
      throw std::runtime_error("CPU DFT executor has no copy-pad processor");
    }
    copy_pad_(request.plane, request.source, request.destination, request.context);
  }

  void process_spatial(DftProcessSpatialRequest request) override {
    if (!process_spatial_) {
      throw std::runtime_error("CPU DFT executor has no spatial processor");
    }
    process_spatial_(
      request.thread_id,
      request.plane,
      request.source,
      request.destination,
      request.context
    );
  }

  void process_temporal(DftProcessTemporalRequest request) override {
    if (!process_temporal_) {
      throw std::runtime_error("CPU DFT executor has no temporal processor");
    }
    process_temporal_(
      request.thread_id,
      request.plane,
      request.source,
      request.destination,
      request.temporal_position,
      request.context
    );
  }

private:
  DFTCopyPadFunction copy_pad_ {nullptr};
  DFTProcessSpatialFunction process_spatial_ {nullptr};
  DFTProcessTemporalFunction process_temporal_ {nullptr};
};

std::unique_ptr<DftExecutor> create_cpu_dft_executor(unsigned opt, DFTClipFormat format) {
  return std::make_unique<CpuDftExecutor>(
    select_cpu_copy_pad(format),
    select_cpu_process_dispatch(opt, format)
  );
}

} // namespace neo_dfttest
