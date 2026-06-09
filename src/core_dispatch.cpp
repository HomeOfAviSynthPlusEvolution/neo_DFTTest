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

static DFTFilterCoefficientsFunction select_scalar_filter(const unsigned ftype, const DFTBlockSettings& block) noexcept {
    if (ftype == 0) {
        if (std::abs(block.f0_beta - 1.0f) < 0.00005f)
            return filter_scalar<0>;
        else if (std::abs(block.f0_beta - 0.5f) < 0.00005f)
            return filter_scalar<6>;
        else
            return filter_scalar<5>;
    } else if (ftype == 1) {
        return filter_scalar<1>;
    } else if (ftype == 2) {
        return filter_scalar<2>;
    } else if (ftype == 3) {
        return filter_scalar<3>;
    }

    return filter_scalar<4>;
}

static DFTKernelDispatch make_scalar_dispatch(const unsigned ftype, const DFTClipFormat& format, const DFTBlockSettings& block) noexcept {
    DFTKernelDispatch kernels {};
    kernels.filter_coefficients = select_scalar_filter(ftype, block);

    if (format.bytes_per_sample == 1) {
        kernels.copy_pad = copyPad<uint8_t>;
        kernels.process_spatial = process_spatial_scalar<uint8_t>;
        kernels.process_temporal = process_temporal_scalar<uint8_t>;
    } else if (format.bytes_per_sample == 2) {
        kernels.copy_pad = copyPad<uint16_t>;
        kernels.process_spatial = process_spatial_scalar<uint16_t>;
        kernels.process_temporal = process_temporal_scalar<uint16_t>;
    } else {
        kernels.copy_pad = copyPad<float>;
        kernels.process_spatial = process_spatial_scalar<float>;
        kernels.process_temporal = process_temporal_scalar<float>;
    }

    return kernels;
}

DFTKernelDispatch selectFunctions(const unsigned ftype, const unsigned opt, const DFTClipFormat& format, const DFTBlockSettings& block) noexcept {
    DFTKernelDispatch kernels = make_scalar_dispatch(ftype, format, block);

    if (opt != 1) {
        DFTKernelDispatch highway = neo_dfttest::make_highway_dispatch(ftype, format, block);
        kernels.filter_coefficients = highway.filter_coefficients;
        kernels.process_spatial = highway.process_spatial;
        kernels.process_temporal = highway.process_temporal;
    }

    return kernels;
}
