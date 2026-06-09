#ifndef NEO_DFTTEST_DFT_EXECUTOR_HPP
#define NEO_DFTTEST_DFT_EXECUTOR_HPP

#include "dft_common.h"

#include <memory>

namespace neo_dfttest {

class DftExecutor {
public:
  virtual ~DftExecutor() = default;

  virtual void process_spatial(
    unsigned int thread_id,
    int plane,
    DFTPlaneBytes src,
    DFTMutablePlaneBytes dst,
    const DFTKernelContext& context
  ) = 0;

  virtual void process_temporal(
    unsigned int thread_id,
    int plane,
    DFTPlaneBytes src,
    DFTMutablePlaneBytes dst,
    int temporal_position,
    const DFTKernelContext& context
  ) = 0;
};

std::unique_ptr<DftExecutor> create_cpu_dft_executor(DFTKernelDispatch dispatch);

} // namespace neo_dfttest

#endif
