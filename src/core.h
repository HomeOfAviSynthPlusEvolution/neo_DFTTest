#ifndef CORE_H
#define CORE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>

#include "dft_common.h"
#include "MersenneTwister.h"

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

// Highway functions
namespace neo_dfttest {
    using FilterFunc = void (*)(float *, const float *, const int, const float *, const float *, const float *);
    FilterFunc GetHighwayFilter(int ftype, float f0beta);
    void GetHighwayFunc0(DFTTestData* d);
    void GetHighwayFunc1(DFTTestData* d);

    // Getters for internal testing
    using Proc0Func_u8 = void (*)(const uint8_t*, const float*, float*, int, int, float);
    Proc0Func_u8 GetHighwayProc0_u8();

    using Proc0Func_u16 = void (*)(const uint16_t*, const float*, float*, int, int, float);
    Proc0Func_u16 GetHighwayProc0_u16();

    using Proc0Func_f32 = void (*)(const float*, const float*, float*, int, int, float);
    Proc0Func_f32 GetHighwayProc0_f32();

    using Proc1Func = void (*)(const float*, const float*, float*, int, int);
    Proc1Func GetHighwayProc1();

    using RemoveMeanFunc = void (*)(float*, const float*, int, float*);
    RemoveMeanFunc GetHighwayRemoveMean();

    using AddMeanFunc = void (*)(float*, int, const float*);
    AddMeanFunc GetHighwayAddMean();

    using CastFunc_u8 = void (*)(const float *, uint8_t *, int, int, int, int, float, int);
    CastFunc_u8 GetHighwayCast_u8();

    using CastFunc_u16 = void (*)(const float *, uint16_t *, int, int, int, int, float, int);
    CastFunc_u16 GetHighwayCast_u16();

    using CastFunc_f32 = void (*)(const float *, float *, int, int, int, int, float, int);
    CastFunc_f32 GetHighwayCast_f32();

    FilterFunc GetHighwayFilter_Type0();
    FilterFunc GetHighwayFilter_Type1();
    FilterFunc GetHighwayFilter_Type2();
    FilterFunc GetHighwayFilter_Type3();
    FilterFunc GetHighwayFilter_Type4();
    FilterFunc GetHighwayFilter_Type5();
    FilterFunc GetHighwayFilter_Type6();
}

#endif
