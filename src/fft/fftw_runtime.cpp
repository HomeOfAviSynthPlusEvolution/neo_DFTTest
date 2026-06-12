#include "fft/fftw_runtime.hpp"

#include <utility>

namespace neo_dfttest::fft {

std::shared_ptr<const FftwRuntime> FftwRuntime::load() {
  return std::shared_ptr<const FftwRuntime>(new FftwRuntime(DynamicLibrary::open_fftw()));
}

FftwRuntime::FftwRuntime(DynamicLibrary library) : library_(std::move(library)) {
  load_symbols();
}

bool FftwRuntime::loaded() const noexcept {
  return library_.loaded();
}

bool FftwRuntime::has_threading() const noexcept {
  return loaded() && fftwf_init_threads && fftwf_plan_with_nthreads;
}

void FftwRuntime::load_symbols() {
  fftwf_destroy_plan = library_.required_symbol<DestroyPlanProc>("fftwf_destroy_plan");
  fftwf_execute_dft_r2c = library_.required_symbol<ExecuteDftR2CProc>("fftwf_execute_dft_r2c");
  fftwf_execute_dft_c2r = library_.required_symbol<ExecuteDftC2RProc>("fftwf_execute_dft_c2r");
  fftwf_plan_dft_r2c_2d = library_.required_symbol<PlanDftR2C2DProc>("fftwf_plan_dft_r2c_2d");
  fftwf_plan_dft_c2r_2d = library_.required_symbol<PlanDftC2R2DProc>("fftwf_plan_dft_c2r_2d");
  fftwf_plan_dft_r2c_3d = library_.required_symbol<PlanDftR2C3DProc>("fftwf_plan_dft_r2c_3d");
  fftwf_plan_dft_c2r_3d = library_.required_symbol<PlanDftC2R3DProc>("fftwf_plan_dft_c2r_3d");
  fftwf_init_threads = library_.optional_symbol<InitThreadsProc>("fftwf_init_threads");
  fftwf_plan_with_nthreads = library_.optional_symbol<PlanWithNThreadsProc>("fftwf_plan_with_nthreads");
}

} // namespace neo_dfttest::fft
