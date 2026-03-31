#ifndef CORE_H
#define CORE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <random>
#include <type_traits>

#include "dft_common.h"

#ifdef VS_TARGET_CPU_X86
#include "vectorclass/vectorclass.h"

template<int type> void filter_sse2(float *, const float *, const int, const float *, const float *, const float *) noexcept;
template<int type> void filter_avx2(float *, const float *, const int, const float *, const float *, const float *) noexcept;

template<typename T> void func_0_sse2(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template<typename T> void func_0_avx2(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;

template<typename T> void func_1_sse2(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
template<typename T> void func_1_avx2(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
#endif

template<int type> void filter_c(float *, const float *, const int, const float *, const float *, const float *) noexcept;
template<typename T> void func_0_c(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template<typename T> void func_1_c(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;

#endif
