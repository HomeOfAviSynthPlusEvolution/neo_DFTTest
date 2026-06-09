#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <optional>
#include <random>
#include <type_traits>
#include <utility>

#include "dft_common.h"

using DFTCopyPadFunction = void (*)(int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) noexcept;
using DFTProcessSpatialFunction = void (*)(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context);
using DFTProcessTemporalFunction = void (*)(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, int temporal_position, const DFTKernelContext& context);

struct DFTCpuProcessDispatch {
    DFTProcessSpatialFunction process_spatial {nullptr};
    DFTProcessTemporalFunction process_temporal {nullptr};
};

template<int type> void filter_scalar(DFTFilterInput input) noexcept;
template<typename T> void process_spatial_scalar(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context);
template<typename T> void process_temporal_scalar(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, int temporal_position, const DFTKernelContext& context);

DFTCopyPadFunction select_cpu_copy_pad(const DFTClipFormat& format) noexcept;
DFTCpuProcessDispatch select_cpu_process_dispatch(unsigned opt, const DFTClipFormat& format) noexcept;

// Highway functions
namespace neo_dfttest {
    DFTCpuProcessDispatch make_highway_process_dispatch(const DFTClipFormat& format);
}
