#pragma once

#include "fft/dynamic_library.hpp"
#include "fft/fftw_abi.hpp"

#include <memory>

namespace neo_dfttest::fft {

class FftwRuntime {
public:
  using DestroyPlanProc = void (*)(fftwf_plan plan);
  using ExecuteDftR2CProc = void (*)(fftwf_plan plan, float* realdata, fftwf_complex* fftsrc);
  using ExecuteDftC2RProc = void (*)(fftwf_plan plan, fftwf_complex* fftsrc, float* realdata);
  using PlanDftR2C2DProc = fftwf_plan (*)(int n0, int n1, float* in, fftwf_complex* out, unsigned flags);
  using PlanDftC2R2DProc = fftwf_plan (*)(int n0, int n1, fftwf_complex* in, float* out, unsigned flags);
  using PlanDftR2C3DProc = fftwf_plan (*)(int n0, int n1, int n2, float* in, fftwf_complex* out, unsigned flags);
  using PlanDftC2R3DProc = fftwf_plan (*)(int n0, int n1, int n2, fftwf_complex* in, float* out, unsigned flags);
  using InitThreadsProc = int (*)();
  using PlanWithNThreadsProc = void (*)(int nthreads);

  static std::shared_ptr<const FftwRuntime> load();

  [[nodiscard]] bool loaded() const noexcept;
  [[nodiscard]] bool has_threading() const noexcept;

  DestroyPlanProc fftwf_destroy_plan {nullptr};
  ExecuteDftR2CProc fftwf_execute_dft_r2c {nullptr};
  ExecuteDftC2RProc fftwf_execute_dft_c2r {nullptr};
  PlanDftR2C2DProc fftwf_plan_dft_r2c_2d {nullptr};
  PlanDftC2R2DProc fftwf_plan_dft_c2r_2d {nullptr};
  PlanDftR2C3DProc fftwf_plan_dft_r2c_3d {nullptr};
  PlanDftC2R3DProc fftwf_plan_dft_c2r_3d {nullptr};
  InitThreadsProc fftwf_init_threads {nullptr};
  PlanWithNThreadsProc fftwf_plan_with_nthreads {nullptr};

private:
  explicit FftwRuntime(DynamicLibrary library);

  void load_symbols();

  DynamicLibrary library_;
};

} // namespace neo_dfttest::fft
