#pragma once

#include <complex>

namespace neo_dfttest::fft {

using fftwf_complex = float[2];
struct fftwf_plan_s;
using fftwf_plan = fftwf_plan_s*;

inline constexpr unsigned FFTW_MEASURE = 0U;
inline constexpr unsigned FFTW_DESTROY_INPUT = 1U << 0;
inline constexpr unsigned FFTW_UNALIGNED = 1U << 1;
inline constexpr unsigned FFTW_CONSERVE_MEMORY = 1U << 2;
inline constexpr unsigned FFTW_EXHAUSTIVE = 1U << 3;
inline constexpr unsigned FFTW_PRESERVE_INPUT = 1U << 4;
inline constexpr unsigned FFTW_PATIENT = 1U << 5;
inline constexpr unsigned FFTW_ESTIMATE = 1U << 6;
inline constexpr unsigned FFTW_WISDOM_ONLY = 1U << 21;

static_assert(sizeof(std::complex<float>) == sizeof(fftwf_complex) &&
                alignof(std::complex<float>) == alignof(fftwf_complex),
              "std::complex<float> and fftwf_complex must have the identical ABI size and alignment.");

} // namespace neo_dfttest::fft
