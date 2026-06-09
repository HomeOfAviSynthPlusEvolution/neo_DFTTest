#ifndef CORE_H
#define CORE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <random>
#include <type_traits>

#include "dft_common.h"

template<int type> void filter_c(DFTFilterInput input) noexcept;
template<typename T> void func_0_c(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) noexcept;
template<typename T> void func_1_c(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, int temporal_position, const DFTKernelContext& context) noexcept;

// Highway functions
namespace neo_dfttest {
    DFTKernelDispatch GetHighwayDispatch(unsigned ftype, const DFTClipFormat& format, const DFTBlockSettings& block);
}

#endif
