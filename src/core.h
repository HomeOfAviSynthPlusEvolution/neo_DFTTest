#ifndef CORE_H
#define CORE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <optional>
#include <random>
#include <type_traits>
#include <utility>

#include "dft_common.h"

template<int type> void filter_scalar(DFTFilterInput input) noexcept;
template<typename T> void process_spatial_scalar(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context);
template<typename T> void process_temporal_scalar(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, int temporal_position, const DFTKernelContext& context);

// Highway functions
namespace neo_dfttest {
    DFTKernelDispatch make_highway_dispatch(DFTFilterPlan filter_plan, const DFTClipFormat& format);
}

#endif
