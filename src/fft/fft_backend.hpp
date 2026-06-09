#ifndef NEO_DFTTEST_FFT_BACKEND_HPP
#define NEO_DFTTEST_FFT_BACKEND_HPP

#include <memory>

#include "fftwlite.h"

namespace neo_dfttest::fft {

using Real = float;
using Complex = fftwf_complex;
using Plan = fftwf_plan;

constexpr unsigned kPatientDestroyInputPlanFlags = FFTW_PATIENT | FFTW_DESTROY_INPUT;

class Backend {
public:
  virtual ~Backend() = default;

  virtual void load() = 0;
  virtual void unload() noexcept = 0;
  virtual bool loaded() const noexcept = 0;

  virtual bool has_threading() const noexcept = 0;
  virtual void set_thread_count(int nthreads) = 0;

  virtual Plan plan_r2c_2d(int n0, int n1, Real* in, Complex* out, unsigned flags) = 0;
  virtual Plan plan_c2r_2d(int n0, int n1, Complex* in, Real* out, unsigned flags) = 0;
  virtual Plan plan_r2c_3d(int n0, int n1, int n2, Real* in, Complex* out, unsigned flags) = 0;
  virtual Plan plan_c2r_3d(int n0, int n1, int n2, Complex* in, Real* out, unsigned flags) = 0;

  virtual void execute_r2c(Plan plan, Real* in, Complex* out) const noexcept = 0;
  virtual void execute_c2r(Plan plan, Complex* in, Real* out) const noexcept = 0;
  virtual void destroy_plan(Plan plan) noexcept = 0;
};

std::unique_ptr<Backend> create_fftw_backend();

} // namespace neo_dfttest::fft

#endif
