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
void filter_c(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept;

template<>
void filter_c<0>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        const float mult = std::max((psd - sigmas[h]) / (psd + 1e-15f), 0.0f);
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<>
void filter_c<1>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        if (psd < sigmas[h])
            dftc[h] = dftc[h + 1] = 0.0f;
    }
}

template<>
void filter_c<2>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 2) {
        dftc[h] *= sigmas[h];
        dftc[h + 1] *= sigmas[h];
    }
}

template<>
void filter_c<3>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
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
void filter_c<4>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1] + 1e-15f;
        const float mult = sigmas[h] * std::sqrt(psd * pmax[h] / ((psd + pmin[h]) * (psd + pmax[h])));
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<>
void filter_c<5>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    const float beta = pmin[0];

    for (int h = 0; h < ccnt; h += 2) {
        const float psd = dftc[h] * dftc[h] + dftc[h + 1] * dftc[h + 1];
        const float mult = std::pow(std::max((psd - sigmas[h]) / (psd + 1e-15f), 0.0f), beta);
        dftc[h] *= mult;
        dftc[h + 1] *= mult;
    }
}

template<>
void filter_c<6>(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
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
void func_0_c(unsigned int thread_id, int plane, const unsigned char * src_ptr, unsigned char * dst_ptr, int dst_stride_bytes, const DFTTestData * d) noexcept {
    float * ebuff = d->scratch.ebuff[thread_id];
    float * dftr = d->scratch.dftr[thread_id];
    fftwf_complex * dftc = d->scratch.dftc[thread_id];
    fftwf_complex * dftc2 = d->scratch.dftc2[thread_id];

    const int width = d->planes.pad_width[plane];
    const int height = d->planes.pad_height[plane];
    const int eheight = d->planes.e_height[plane];
    const int srcStride = d->planes.pad_stride[plane] / sizeof(T);
    const int ebpStride = d->planes.e_stride[plane];
    const T * srcp = reinterpret_cast<const T *>(src_ptr);
    float * ebpSaved = ebuff;

    memset(ebuff, 0, ebpStride * height * sizeof(float));

    for (int y = 0; y < eheight; y += d->derived.step) {
        for (int x = 0; x <= width - d->block.spatial_size; x += d->derived.step) {
            proc0(srcp + x, d->coefficients.window, dftr, srcStride, d->block.spatial_size, d->sample.divisor);

            d->fft.api->fftwf_execute_dft_r2c(d->fft.forward, dftr, dftc);
            if (d->block.zero_mean)
                removeMean(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(d->coefficients.window_dft), d->derived.coefficient_count, reinterpret_cast<float *>(dftc2));

            d->kernels.filter_coefficients(reinterpret_cast<float *>(dftc), d->coefficients.sigmas, d->derived.coefficient_count, d->derived.custom_f0_beta ? &d->block.f0_beta : d->coefficients.pmins, d->coefficients.pmaxs, d->coefficients.sigmas2);

            if (d->block.zero_mean)
                addMean(reinterpret_cast<float *>(dftc), d->derived.coefficient_count, reinterpret_cast<const float *>(dftc2));
            d->fft.api->fftwf_execute_dft_c2r(d->fft.inverse, dftc, dftr);

            if (d->derived.transform_type & 1) // spatial overlapping
                proc1(dftr, d->coefficients.window, ebpSaved + x, d->block.spatial_size, ebpStride);
            else
                ebpSaved[x + d->derived.spatial_center * ebpStride + d->derived.spatial_center] = dftr[d->derived.spatial_center * d->block.spatial_size + d->derived.spatial_center] * d->coefficients.window[d->derived.spatial_center * d->block.spatial_size + d->derived.spatial_center];
        }

        srcp += srcStride * d->derived.step;
        ebpSaved += ebpStride * d->derived.step;
    }

    int dstWidth = d->planes.width[plane];
    int dstHeight = d->planes.height[plane];
    int dstStride = dst_stride_bytes / sizeof(T);
    T * dstp = reinterpret_cast<T *>(dst_ptr);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (d->block.dither_mode > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->sample.multiplier, d->sample.peak, d->block.dither_mode, *d->scratch.rngs[thread_id], d->scratch.dither_buffers[thread_id]);
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->sample.multiplier, d->sample.peak);
}

template<typename T>
void func_1_c(unsigned int thread_id, int plane, const unsigned char * src_ptr, unsigned char * dst_ptr, int dst_stride_bytes, const int pos, const DFTTestData * d) noexcept {
    float * ebuff = d->scratch.ebuff[thread_id];
    float * dftr = d->scratch.dftr[thread_id];
    fftwf_complex * dftc = d->scratch.dftc[thread_id];
    fftwf_complex * dftc2 = d->scratch.dftc2[thread_id];

    const int width = d->planes.pad_width[plane];
    const int height = d->planes.pad_height[plane];
    const int eheight = d->planes.e_height[plane];
    const int srcStride = d->planes.pad_stride[plane] / sizeof(T);
    const int ebpStride = d->planes.e_stride[plane];

    const T * srcp[15] = {};
    for (int i = 0; i < d->block.temporal_size; i++)
        srcp[i] = reinterpret_cast<const T *>(src_ptr + d->planes.pad_block_size[plane] * i);

    memset(ebuff, 0, ebpStride * height * sizeof(float));

    for (int y = 0; y < eheight; y += d->derived.step) {
        for (int x = 0; x <= width - d->block.spatial_size; x += d->derived.step) {
            for (int z = 0; z < d->block.temporal_size; z++)
                proc0(srcp[z] + x, d->coefficients.window + d->derived.block_area * z, dftr + d->derived.block_area * z, srcStride, d->block.spatial_size, d->sample.divisor);

            d->fft.api->fftwf_execute_dft_r2c(d->fft.forward, dftr, dftc);
            if (d->block.zero_mean)
                removeMean(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(d->coefficients.window_dft), d->derived.coefficient_count, reinterpret_cast<float *>(dftc2));

            d->kernels.filter_coefficients(reinterpret_cast<float *>(dftc), d->coefficients.sigmas, d->derived.coefficient_count, d->derived.custom_f0_beta ? &d->block.f0_beta : d->coefficients.pmins, d->coefficients.pmaxs, d->coefficients.sigmas2);

            if (d->block.zero_mean)
                addMean(reinterpret_cast<float *>(dftc), d->derived.coefficient_count, reinterpret_cast<const float *>(dftc2));
            d->fft.api->fftwf_execute_dft_c2r(d->fft.inverse, dftc, dftr);

            if (d->derived.transform_type & 1) // spatial overlapping
                proc1(dftr + pos * d->derived.block_area, d->coefficients.window + pos * d->derived.block_area, ebuff + y * ebpStride + x, d->block.spatial_size, ebpStride);
            else
                ebuff[(y + d->derived.spatial_center) * ebpStride + x + d->derived.spatial_center] = dftr[pos * d->derived.block_area + d->derived.spatial_center * d->block.spatial_size + d->derived.spatial_center] * d->coefficients.window[pos * d->derived.block_area + d->derived.spatial_center * d->block.spatial_size + d->derived.spatial_center];
        }

        for (int q = 0; q < d->block.temporal_size; q++)
            srcp[q] += srcStride * d->derived.step;
    }

    int dstWidth = d->planes.width[plane];
    int dstHeight = d->planes.height[plane];
    int dstStride = dst_stride_bytes / sizeof(T);
    T * dstp = reinterpret_cast<T *>(dst_ptr);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (d->block.dither_mode > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->sample.multiplier, d->sample.peak, d->block.dither_mode, *d->scratch.rngs[thread_id], d->scratch.dither_buffers[thread_id]);
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->sample.multiplier, d->sample.peak);
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

template void func_0_c<uint8_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template void func_0_c<uint16_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template void func_0_c<float>(unsigned int, int, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template void func_1_c<uint8_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
template void func_1_c<uint16_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
template void func_1_c<float>(unsigned int, int, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
