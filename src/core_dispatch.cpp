#include "core.h"

template<typename T>
static void copyPad(int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) noexcept {
    int srcWidth = context.planes.width[plane];
    int srcHeight = context.planes.height[plane];
    int dstWidth = context.planes.pad_width[plane];
    int dstHeight = context.planes.pad_height[plane];
    int dstStrideBytes = dst.stride_bytes;
    int dstStride = dst.stride_bytes / sizeof(T);

    const int offy = (dstHeight - srcHeight) / 2;
    const int offx = (dstWidth - srcWidth) / 2;

    const unsigned char * scrp0 = src.data;
    unsigned char * dstp0 = dst.data + dstStrideBytes * offy + offx * sizeof(T);
    for (int h = 0; h < srcHeight; h++)
    {
        memcpy(dstp0, scrp0, srcWidth * sizeof(T));
        scrp0 += src.stride_bytes;
        dstp0 += dstStrideBytes;
    }
    
    T * dstp = reinterpret_cast<T *>(dst.data) + dstStride * offy;

    for (int y = offy; y < srcHeight + offy; y++) {
        int w = offx * 2;
        for (int x = 0; x < offx; x++, w--)
            dstp[x] = dstp[w];

        w = offx + srcWidth - 2;
        for (int x = offx + srcWidth; x < dstWidth; x++, w--)
            dstp[x] = dstp[w];

        dstp += dstStride;
    }

    int w = offy * 2;
    for (int y = 0; y < offy; y++, w--)
        memcpy(dst.data + dstStrideBytes * y, dst.data + dstStrideBytes * w, dstWidth * sizeof(T));

    w = offy + srcHeight - 2;
    for (int y = offy + srcHeight; y < dstHeight; y++, w--)
        memcpy(dst.data + dstStrideBytes * y, dst.data + dstStrideBytes * w, dstWidth * sizeof(T));
}

static DFTFilterPlan make_filter_plan(const unsigned ftype, const DFTBlockSettings& block) noexcept {
    if (ftype == 0) {
        if (std::abs(block.f0_beta - 1.0f) < 0.00005f)
            return DFTFilterPlan{DFTFilterKind::wiener, false};
        else if (std::abs(block.f0_beta - 0.5f) < 0.00005f)
            return DFTFilterPlan{DFTFilterKind::wiener_sqrt, true};
        else
            return DFTFilterPlan{DFTFilterKind::wiener_power, true};
    } else if (ftype == 1) {
        return DFTFilterPlan{DFTFilterKind::hard_threshold, false};
    } else if (ftype == 2) {
        return DFTFilterPlan{DFTFilterKind::multiplier, false};
    } else if (ftype == 3) {
        return DFTFilterPlan{DFTFilterKind::range_multiplier, false};
    }

    return DFTFilterPlan{DFTFilterKind::range_wiener, false};
}

static DFTFilterCoefficientsFunction select_scalar_filter(const DFTFilterPlan filter_plan) noexcept {
    switch (filter_plan.kind) {
        case DFTFilterKind::wiener:
            return filter_scalar<0>;
        case DFTFilterKind::hard_threshold:
            return filter_scalar<1>;
        case DFTFilterKind::multiplier:
            return filter_scalar<2>;
        case DFTFilterKind::range_multiplier:
            return filter_scalar<3>;
        case DFTFilterKind::range_wiener:
            return filter_scalar<4>;
        case DFTFilterKind::wiener_power:
            return filter_scalar<5>;
        case DFTFilterKind::wiener_sqrt:
            return filter_scalar<6>;
    }

    return filter_scalar<0>;
}

static DFTKernelDispatch make_scalar_dispatch(const DFTFilterPlan filter_plan, const DFTClipFormat& format) noexcept {
    DFTKernelDispatch kernels {};
    kernels.filter_coefficients = select_scalar_filter(filter_plan);
    kernels.filter_plan = filter_plan;

    if (format.bytes_per_sample == 1) {
        kernels.copy_pad = copyPad<uint8_t>;
    } else if (format.bytes_per_sample == 2) {
        kernels.copy_pad = copyPad<uint16_t>;
    } else {
        kernels.copy_pad = copyPad<float>;
    }

    return kernels;
}

static DFTCpuProcessDispatch make_scalar_process_dispatch(const DFTClipFormat& format) noexcept {
    if (format.bytes_per_sample == 1) {
        return DFTCpuProcessDispatch{process_spatial_scalar<uint8_t>, process_temporal_scalar<uint8_t>};
    }
    if (format.bytes_per_sample == 2) {
        return DFTCpuProcessDispatch{process_spatial_scalar<uint16_t>, process_temporal_scalar<uint16_t>};
    }

    return DFTCpuProcessDispatch{process_spatial_scalar<float>, process_temporal_scalar<float>};
}

DFTKernelDispatch selectFunctions(const unsigned ftype, const unsigned opt, const DFTClipFormat& format, const DFTBlockSettings& block) noexcept {
    const DFTFilterPlan filter_plan = make_filter_plan(ftype, block);
    DFTKernelDispatch kernels = make_scalar_dispatch(filter_plan, format);

    if (opt != 1) {
        DFTKernelDispatch highway = neo_dfttest::make_highway_dispatch(filter_plan, format);
        kernels.filter_coefficients = highway.filter_coefficients;
        kernels.filter_plan = highway.filter_plan;
    }

    return kernels;
}

DFTCpuProcessDispatch select_cpu_process_dispatch(const unsigned opt, const DFTClipFormat& format) noexcept {
    if (opt != 1) {
        return neo_dfttest::make_highway_process_dispatch(format);
    }

    return make_scalar_process_dispatch(format);
}
