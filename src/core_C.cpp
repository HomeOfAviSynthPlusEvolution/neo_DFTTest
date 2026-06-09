/*
**   VapourSynth port by HolyWu
**
**                    dfttest v1.8 for Avisynth 2.5.x
**
**   2D/3D frequency domain denoiser.
**
**   Copyright (C) 2007-2010 Kevin Stone
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "core.h"


template<typename T>
static inline void proc0(const T * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept;

template<>
inline void proc0(const uint8_t * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept {
    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v++)
            d[v] = s0[v] * s1[v];

        s0 += p0;
        s1 += p1;
        d += p1;
    }
}

template<>
inline void proc0(const uint16_t * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept {
    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v++)
            d[v] = s0[v] * divisor * s1[v];

        s0 += p0;
        s1 += p1;
        d += p1;
    }
}

template<>
inline void proc0(const float * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept {
    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v++)
            d[v] = s0[v] * 255.0f * s1[v];

        s0 += p0;
        s1 += p1;
        d += p1;
    }
}

static inline void proc1(const float * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1) noexcept {
    for (int u = 0; u < p0; u++) {
        for (int v = 0; v < p0; v++)
            d[v] += s0[v] * s1[v];

        s0 += p0;
        s1 += p0;
        d += p1;
    }
}

static inline void removeMean(float * VS_RESTRICT dftc, const float * dftgc, const int ccnt, float * VS_RESTRICT dftc2) noexcept {
    const float gf = dftc[0] / dftgc[0];

    for (int h = 0; h < ccnt; h += 2) {
        dftc2[h] = gf * dftgc[h];
        dftc2[h + 1] = gf * dftgc[h + 1];
        dftc[h] -= dftc2[h];
        dftc[h + 1] -= dftc2[h + 1];
    }
}

static inline void addMean(float * VS_RESTRICT dftc, const int ccnt, const float * dftc2) noexcept {
    for (int h = 0; h < ccnt; h += 2) {
        dftc[h] += dftc2[h];
        dftc[h + 1] += dftc2[h + 1];
    }
}

template<int type>
void filter_c(DFTFilterInput input) noexcept;

template<>
void filter_c<0>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        const float mult = std::max((psd - sigmas[h]) / (psd + 1e-15f), 0.0f);
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<>
void filter_c<1>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        if (psd < sigmas[h])
            dftc[h] = dftc[h + 1] = 0.0f;
    }
}

template<>
void filter_c<2>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        dftc[h] *= sigmas[h];
        dftc[h + 1] *= sigmas[h];
    }
}

template<>
void filter_c<3>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const float* pmin = input.pmins.data;
    const float* pmax = input.pmaxs.data;
    const float* sigmas2 = input.sigmas2.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];

        if (psd >= pmin[h] && psd <= pmax[h]) {
            dftc[h] *= sigmas[h];
            dftc[h + 1] *= sigmas[h];
        } else {
            dftc[h] *= sigmas2[h];
            dftc[h + 1] *= sigmas2[h];
        }
    }
}

template<>
void filter_c<4>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const float* pmin = input.pmins.data;
    const float* pmax = input.pmaxs.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1] + 1e-15f;
        const float mult = sigmas[h] * std::sqrt(psd * pmax[h] / ((psd + pmin[h]) * (psd + pmax[h])));
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<>
void filter_c<5>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const int ccnt = input.coefficients.size;
    const float beta = input.pmins.data[0];

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        const float mult = std::pow(std::max((psd - sigmas[h]) / (psd + 1e-15f), 0.0f), beta);
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<>
void filter_c<6>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        const float mult = std::sqrt(std::max((psd - sigmas[h]) / (psd + 1e-15f), 0.0f));
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<typename T>
void cast(const float * ebp, T * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak) noexcept;

template<>
void cast(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
          const float multiplier, const int peak) noexcept {
    for (int y = 0; y < dstHeight; y++) {
        for (int x = 0; x < dstWidth; x++)
            dstp[x] = std::min(std::max(static_cast<int>(ebp[x] + 0.5f), 0), 255);

        ebp += ebpStride;
        dstp += dstStride;
    }
}

template<>
void cast(const float * ebp, uint16_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
          const float multiplier, const int peak) noexcept {
    for (int y = 0; y < dstHeight; y++) {
        for (int x = 0; x < dstWidth; x++)
            dstp[x] = std::min(std::max(static_cast<int>(ebp[x] * multiplier + 0.5f), 0), peak);

        ebp += ebpStride;
        dstp += dstStride;
    }
}

template<>
void cast(const float * ebp, float * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
          const float multiplier, const int peak) noexcept {
    for (int y = 0; y < dstHeight; y++) {
        for (int x = 0; x < dstWidth; x++)
            dstp[x] = ebp[x] * (1.0f / 255.0f);

        ebp += ebpStride;
        dstp += dstStride;
    }
}

template<typename T>
static inline void dither(const float * ebp, T * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode, std::mt19937 &rng, float *dither_buff) noexcept {
    cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak);
}

template<>
inline void dither<uint8_t>(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode, std::mt19937 &rng, float *dither_buff) noexcept {
    float* dc = dither_buff;
    float* dn = dither_buff + dstWidth;
    const float scale = (dither_mode - 1) + 0.5f;
    const float off = scale * 0.5f;
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    memset(dc, 0, dstWidth * sizeof(float));
    for (int y = 0; y < dstHeight; ++y) {
        memset(dn, 0, dstWidth * sizeof(float));
        for (int x = 0; x < dstWidth; ++x) {
            const int v = dither_mode == 1 ?
                (int)(ebp[x] + dc[x] + 0.5f) :
                (int)(ebp[x] + dist(rng) * scale - off + dc[x] + 0.5f);
            dstp[x] = std::min(std::max(v, 0), 255);
            const float qerror = ebp[x] - dstp[x];
            if (x != 0)
                dn[x - 1] += qerror * 0.1875f;
            dn[x] += qerror * 0.3125f;
            if (x != dstWidth - 1) {
                dc[x + 1] += qerror * 0.4375f;
                dn[x + 1] += qerror * 0.0625f;
            }
        }
        ebp += ebpStride;
        dstp += dstStride;
        float* tn = dn;
        dn = dc;
        dc = tn;
    }
}

template<typename T>
void func_0_c(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) noexcept {
    float * ebuff = context.scratch.slots[thread_id].ebuff.data();
    float * dftr = context.scratch.slots[thread_id].dftr.data();
    auto* dftc = context.scratch.slots[thread_id].dftc.data();
    auto* dftc2 = context.scratch.slots[thread_id].dftc2.data();

    const int width = context.planes.pad_width[plane];
    const int height = context.planes.pad_height[plane];
    const int eheight = context.planes.e_height[plane];
    const int srcStride = context.planes.pad_stride[plane] / sizeof(T);
    const int ebpStride = context.planes.e_stride[plane];
    const T * srcp = reinterpret_cast<const T *>(src.data);
    float * ebpSaved = ebuff;

    memset(ebuff, 0, ebpStride * height * sizeof(float));

    for (int y = 0; y < eheight; y += context.derived.step) {
        for (int x = 0; x <= width - context.block.spatial_size; x += context.derived.step) {
            proc0(srcp + x, context.coefficients.window.data(), dftr, srcStride, context.block.spatial_size, context.sample.divisor);

            context.fft.backend->execute_r2c(context.fft.forward, dftr, dftc);
            if (context.block.zero_mean)
                removeMean(complex_float_data(dftc), complex_float_data(context.coefficients.window_dft.data()), context.derived.coefficient_count, complex_float_data(dftc2));

            context.filter_coefficients(make_filter_input(dftc, context));

            if (context.block.zero_mean)
                addMean(complex_float_data(dftc), context.derived.coefficient_count, complex_float_data(dftc2));
            context.fft.backend->execute_c2r(context.fft.inverse, dftc, dftr);

            if (context.derived.transform_type & 1) // spatial overlapping
                proc1(dftr, context.coefficients.window.data(), ebpSaved + x, context.block.spatial_size, ebpStride);
            else
                ebpSaved[x + context.derived.spatial_center * ebpStride + context.derived.spatial_center] = dftr[context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center] * context.coefficients.window.data()[context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center];
        }

        srcp += srcStride * context.derived.step;
        ebpSaved += ebpStride * context.derived.step;
    }

    int dstWidth = context.planes.width[plane];
    int dstHeight = context.planes.height[plane];
    int dstStride = dst.stride_bytes / sizeof(T);
    T * dstp = reinterpret_cast<T *>(dst.data);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (context.block.dither_mode > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak, context.block.dither_mode, *context.scratch.slots[thread_id].rng, context.scratch.slots[thread_id].dither_buffer.data());
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak);
}

template<typename T>
void func_1_c(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const int pos, const DFTKernelContext& context) noexcept {
    float * ebuff = context.scratch.slots[thread_id].ebuff.data();
    float * dftr = context.scratch.slots[thread_id].dftr.data();
    auto* dftc = context.scratch.slots[thread_id].dftc.data();
    auto* dftc2 = context.scratch.slots[thread_id].dftc2.data();

    const int width = context.planes.pad_width[plane];
    const int height = context.planes.pad_height[plane];
    const int eheight = context.planes.e_height[plane];
    const int srcStride = context.planes.pad_stride[plane] / sizeof(T);
    const int ebpStride = context.planes.e_stride[plane];

    const T * srcp[15] = {};
    for (int i = 0; i < context.block.temporal_size; i++)
        srcp[i] = reinterpret_cast<const T *>(src.data + context.planes.pad_block_size[plane] * i);

    memset(ebuff, 0, ebpStride * height * sizeof(float));

    for (int y = 0; y < eheight; y += context.derived.step) {
        for (int x = 0; x <= width - context.block.spatial_size; x += context.derived.step) {
            for (int z = 0; z < context.block.temporal_size; z++)
                proc0(srcp[z] + x, context.coefficients.window.data() + context.derived.block_area * z, dftr + context.derived.block_area * z, srcStride, context.block.spatial_size, context.sample.divisor);

            context.fft.backend->execute_r2c(context.fft.forward, dftr, dftc);
            if (context.block.zero_mean)
                removeMean(complex_float_data(dftc), complex_float_data(context.coefficients.window_dft.data()), context.derived.coefficient_count, complex_float_data(dftc2));

            context.filter_coefficients(make_filter_input(dftc, context));

            if (context.block.zero_mean)
                addMean(complex_float_data(dftc), context.derived.coefficient_count, complex_float_data(dftc2));
            context.fft.backend->execute_c2r(context.fft.inverse, dftc, dftr);

            if (context.derived.transform_type & 1) // spatial overlapping
                proc1(dftr + pos * context.derived.block_area, context.coefficients.window.data() + pos * context.derived.block_area, ebuff + y * ebpStride + x, context.block.spatial_size, ebpStride);
            else
                ebuff[(y + context.derived.spatial_center) * ebpStride + x + context.derived.spatial_center] = dftr[pos * context.derived.block_area + context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center] * context.coefficients.window.data()[pos * context.derived.block_area + context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center];
        }

        for (int q = 0; q < context.block.temporal_size; q++)
            srcp[q] += srcStride * context.derived.step;
    }

    int dstWidth = context.planes.width[plane];
    int dstHeight = context.planes.height[plane];
    int dstStride = dst.stride_bytes / sizeof(T);
    T * dstp = reinterpret_cast<T *>(dst.data);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (context.block.dither_mode > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak, context.block.dither_mode, *context.scratch.slots[thread_id].rng, context.scratch.slots[thread_id].dither_buffer.data());
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak);
}

template<typename T>
void proc0_c(const T * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept {
    proc0(s0, s1, d, p0, p1, divisor);
}
template void proc0_c<uint8_t>(const uint8_t * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept;
template void proc0_c<uint16_t>(const uint16_t * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept;
template void proc0_c<float>(const float * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept;

void removeMean_c(float * VS_RESTRICT dftc, const float * dftgc, const int ccnt, float * VS_RESTRICT dftc2) noexcept {
    removeMean(dftc, dftgc, ccnt, dftc2);
}

void dither_c(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode, std::mt19937 &rng, float *dither_buff) noexcept {
    dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak, dither_mode, rng, dither_buff);
}

template void func_0_c<uint8_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, const DFTKernelContext&) noexcept;
template void func_0_c<uint16_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, const DFTKernelContext&) noexcept;
template void func_0_c<float>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, const DFTKernelContext&) noexcept;
template void func_1_c<uint8_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, int, const DFTKernelContext&) noexcept;
template void func_1_c<uint16_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, int, const DFTKernelContext&) noexcept;
template void func_1_c<float>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, int, const DFTKernelContext&) noexcept;
