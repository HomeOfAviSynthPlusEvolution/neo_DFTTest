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
#include "core_batch_pipeline.hpp"


template<typename T>
static inline void load_windowed_block(DFTConstSampleBlock<T> source, DFTConstFloatSpan window, DFTMutableFloatSpan destination, const float divisor) noexcept;

template<>
inline void load_windowed_block(DFTConstSampleBlock<uint8_t> source, DFTConstFloatSpan window, DFTMutableFloatSpan destination, const float divisor) noexcept {
    const uint8_t* s0 = source.data;
    const float* s1 = window.data;
    float* d = destination.data;
    for (int u = 0; u < source.size; u++) {
        for (int v = 0; v < source.size; v++)
            d[v] = s0[v] * s1[v];

        s0 += source.stride_elements;
        s1 += source.size;
        d += source.size;
    }
}

template<>
inline void load_windowed_block(DFTConstSampleBlock<uint16_t> source, DFTConstFloatSpan window, DFTMutableFloatSpan destination, const float divisor) noexcept {
    const uint16_t* s0 = source.data;
    const float* s1 = window.data;
    float* d = destination.data;
    for (int u = 0; u < source.size; u++) {
        for (int v = 0; v < source.size; v++)
            d[v] = s0[v] * divisor * s1[v];

        s0 += source.stride_elements;
        s1 += source.size;
        d += source.size;
    }
}

template<>
inline void load_windowed_block(DFTConstSampleBlock<float> source, DFTConstFloatSpan window, DFTMutableFloatSpan destination, const float divisor) noexcept {
    const float* s0 = source.data;
    const float* s1 = window.data;
    float* d = destination.data;
    for (int u = 0; u < source.size; u++) {
        for (int v = 0; v < source.size; v++)
            d[v] = s0[v] * 255.0f * s1[v];

        s0 += source.stride_elements;
        s1 += source.size;
        d += source.size;
    }
}

static inline void accumulate_overlap(const float * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1) noexcept {
    for (int u = 0; u < p0; u++) {
        for (int v = 0; v < p0; v++)
            d[v] += s0[v] * s1[v];

        s0 += p0;
        s1 += p0;
        d += p1;
    }
}

static inline void remove_mean(DFTMutableFloatSpan coefficients, DFTConstFloatSpan reference, DFTMutableFloatSpan removed) noexcept {
    float* VS_RESTRICT dftc = coefficients.data;
    const float* dftgc = reference.data;
    float* VS_RESTRICT dftc2 = removed.data;
    const int ccnt = coefficients.size;
    const float gf = dftc[0] / dftgc[0];

    for (int h = 0; h < ccnt; h += 2) {
        dftc2[h] = gf * dftgc[h];
        dftc2[h + 1] = gf * dftgc[h + 1];
        dftc[h] -= dftc2[h];
        dftc[h + 1] -= dftc2[h + 1];
    }
}

static inline void add_mean(DFTMutableFloatSpan coefficients, DFTConstFloatSpan removed) noexcept {
    float* VS_RESTRICT dftc = coefficients.data;
    const float* dftc2 = removed.data;
    const int ccnt = coefficients.size;
    for (int h = 0; h < ccnt; h += 2) {
        dftc[h] += dftc2[h];
        dftc[h + 1] += dftc2[h + 1];
    }
}

template<int type>
void filter_scalar(DFTFilterInput input) noexcept;

template<>
void filter_scalar<0>(DFTFilterInput input) noexcept {
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
void filter_scalar<1>(DFTFilterInput input) noexcept {
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
void filter_scalar<2>(DFTFilterInput input) noexcept {
    float* dftc = input.coefficients.data;
    const float* sigmas = input.sigmas.data;
    const int ccnt = input.coefficients.size;

    for (int h = 0; h < ccnt; h += 2) {
        dftc[h] *= sigmas[h];
        dftc[h + 1] *= sigmas[h];
    }
}

template<>
void filter_scalar<3>(DFTFilterInput input) noexcept {
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
void filter_scalar<4>(DFTFilterInput input) noexcept {
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
void filter_scalar<5>(DFTFilterInput input) noexcept {
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
void filter_scalar<6>(DFTFilterInput input) noexcept {
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

static void filter_coefficients_scalar(const DFTFilterPlan filter_plan, DFTFilterInput input) noexcept {
    switch (filter_plan.kind) {
        case DFTFilterKind::wiener:
            filter_scalar<0>(input);
            return;
        case DFTFilterKind::hard_threshold:
            filter_scalar<1>(input);
            return;
        case DFTFilterKind::multiplier:
            filter_scalar<2>(input);
            return;
        case DFTFilterKind::range_multiplier:
            filter_scalar<3>(input);
            return;
        case DFTFilterKind::range_wiener:
            filter_scalar<4>(input);
            return;
        case DFTFilterKind::wiener_power:
            filter_scalar<5>(input);
            return;
        case DFTFilterKind::wiener_sqrt:
            filter_scalar<6>(input);
            return;
    }

    filter_scalar<0>(input);
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
                 const float multiplier, const int peak, const int dither_mode, std::mt19937 *rng, float *dither_buff) noexcept {
    cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak);
}

template<>
inline void dither<uint8_t>(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode, std::mt19937 *rng, float *dither_buff) noexcept {
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
                (int)(ebp[x] + dist(*rng) * scale - off + dc[x] + 0.5f);
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
void process_spatial_scalar(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) {
    const auto workspace = neo_dfttest::dft_thread_workspace(context, thread_id);
    float* ebuff = workspace.accumulation;

    const int width = context.planes.pad_width[plane];
    const int height = context.planes.pad_height[plane];
    const int eheight = context.planes.e_height[plane];
    const auto source = dft_const_sample_plane<T>(src);
    const int srcStride = source.stride_elements;
    const int ebpStride = context.planes.e_stride[plane];
    const T * srcp = source.data;

    neo_dfttest::clear_accumulation(workspace, static_cast<std::size_t>(ebpStride) * height);

    neo_dfttest::run_dft_batch_pipeline(
        context,
        workspace.batch,
        width,
        eheight,
        [&](int y, int x, float* dftr) noexcept {
            load_windowed_block<T>(
                DFTConstSampleBlock<T>{srcp + y * srcStride + x, srcStride, context.block.spatial_size},
                DFTConstFloatSpan{context.coefficients.window.data, context.derived.block_area},
                DFTMutableFloatSpan{dftr, context.derived.block_area},
                context.sample.divisor
            );
        },
        [&](neo_dfttest::DFTCompletedBatch ready) {
            for (int index = 0; index < ready.batch.count; ++index) {
                auto* dftc = ready.coefficients.block(index);
                auto* dftc2 = ready.removed_mean.block(index);
                if (context.block.zero_mean)
                    remove_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(context.coefficients.window_dft.data), context.derived.coefficient_count},
                        DFTMutableFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );

                filter_coefficients_scalar(context.filter_plan, make_filter_input(dftc, context));

                if (context.block.zero_mean)
                    add_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );
            }

            neo_dfttest::DftFftOperations{context}.submit_inverse_and_wait(ready.coefficients, ready.real);

            float* output_row = ebuff + ready.y * ebpStride;
            for (int index = 0; index < ready.batch.count; ++index) {
                float* dftr = ready.real.block(index);
                const int block_x = dft_block_job(ready.batch, index).x;
                if (context.derived.transform_type & 1) // spatial overlapping
                    accumulate_overlap(dftr, context.coefficients.window.data, output_row + block_x, context.block.spatial_size, ebpStride);
                else
                    output_row[block_x + context.derived.spatial_center * ebpStride + context.derived.spatial_center] = dftr[context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center] * context.coefficients.window.data[context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center];
            }
        }
    );

    int dstWidth = context.planes.width[plane];
    int dstHeight = context.planes.height[plane];
    const auto destination = dft_mutable_sample_plane<T>(dst);
    int dstStride = destination.stride_elements;
    T * dstp = destination.data;
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (context.block.dither_mode > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak, context.block.dither_mode, workspace.dither_rng, workspace.dither_buffer);
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak);
}

template<typename T>
void process_temporal_scalar(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const int pos, const DFTKernelContext& context) {
    const auto workspace = neo_dfttest::dft_thread_workspace(context, thread_id);
    float* ebuff = workspace.accumulation;

    const int width = context.planes.pad_width[plane];
    const int height = context.planes.pad_height[plane];
    const int eheight = context.planes.e_height[plane];
    const auto first_source = dft_const_sample_plane<T>(src);
    const int srcStride = first_source.stride_elements;
    const int ebpStride = context.planes.e_stride[plane];

    const T * srcp[15] = {};
    for (int i = 0; i < context.block.temporal_size; i++) {
        const auto source = dft_const_sample_plane<T>(
            DFTPlaneBytes{src.data + context.planes.pad_block_size[plane] * i, src.stride_bytes}
        );
        srcp[i] = source.data;
    }

    neo_dfttest::clear_accumulation(workspace, static_cast<std::size_t>(ebpStride) * height);

    neo_dfttest::run_dft_batch_pipeline(
        context,
        workspace.batch,
        width,
        eheight,
        [&](int y, int x, float* dftr) noexcept {
            for (int z = 0; z < context.block.temporal_size; z++)
                load_windowed_block<T>(
                    DFTConstSampleBlock<T>{srcp[z] + y * srcStride + x, srcStride, context.block.spatial_size},
                    DFTConstFloatSpan{context.coefficients.window.data + context.derived.block_area * z, context.derived.block_area},
                    DFTMutableFloatSpan{dftr + context.derived.block_area * z, context.derived.block_area},
                    context.sample.divisor
                );
        },
        [&](neo_dfttest::DFTCompletedBatch ready) {
            for (int index = 0; index < ready.batch.count; ++index) {
                auto* dftc = ready.coefficients.block(index);
                auto* dftc2 = ready.removed_mean.block(index);
                if (context.block.zero_mean)
                    remove_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(context.coefficients.window_dft.data), context.derived.coefficient_count},
                        DFTMutableFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );

                filter_coefficients_scalar(context.filter_plan, make_filter_input(dftc, context));

                if (context.block.zero_mean)
                    add_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );
            }

            neo_dfttest::DftFftOperations{context}.submit_inverse_and_wait(ready.coefficients, ready.real);

            for (int index = 0; index < ready.batch.count; ++index) {
                float* dftr = ready.real.block(index);
                const int block_x = dft_block_job(ready.batch, index).x;
                if (context.derived.transform_type & 1) // spatial overlapping
                    accumulate_overlap(dftr + pos * context.derived.block_area, context.coefficients.window.data + pos * context.derived.block_area, ebuff + ready.y * ebpStride + block_x, context.block.spatial_size, ebpStride);
                else
                    ebuff[(ready.y + context.derived.spatial_center) * ebpStride + block_x + context.derived.spatial_center] = dftr[pos * context.derived.block_area + context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center] * context.coefficients.window.data[pos * context.derived.block_area + context.derived.spatial_center * context.block.spatial_size + context.derived.spatial_center];
            }
        }
    );

    int dstWidth = context.planes.width[plane];
    int dstHeight = context.planes.height[plane];
    const auto destination = dft_mutable_sample_plane<T>(dst);
    int dstStride = destination.stride_elements;
    T * dstp = destination.data;
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (context.block.dither_mode > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak, context.block.dither_mode, workspace.dither_rng, workspace.dither_buffer);
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, context.sample.multiplier, context.sample.peak);
}

template<typename T>
void load_windowed_block_scalar(DFTConstSampleBlock<T> source, DFTConstFloatSpan window, DFTMutableFloatSpan destination, float divisor) noexcept {
    load_windowed_block<T>(source, window, destination, divisor);
}
template void load_windowed_block_scalar<uint8_t>(DFTConstSampleBlock<uint8_t>, DFTConstFloatSpan, DFTMutableFloatSpan, float) noexcept;
template void load_windowed_block_scalar<uint16_t>(DFTConstSampleBlock<uint16_t>, DFTConstFloatSpan, DFTMutableFloatSpan, float) noexcept;
template void load_windowed_block_scalar<float>(DFTConstSampleBlock<float>, DFTConstFloatSpan, DFTMutableFloatSpan, float) noexcept;

void remove_mean_scalar(DFTMutableFloatSpan coefficients, DFTConstFloatSpan reference, DFTMutableFloatSpan removed) noexcept {
    remove_mean(coefficients, reference, removed);
}

void dither_u8_scalar(DFTConstFloatPlane source, DFTMutableBytePlane destination, DFTDitherContext context) noexcept {
    dither(
        source.data,
        destination.data,
        destination.width,
        destination.height,
        destination.stride_elements,
        source.stride_elements,
        context.multiplier,
        context.peak,
        context.mode,
        context.rng,
        context.buffer.data
    );
}

template void process_spatial_scalar<uint8_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, const DFTKernelContext&);
template void process_spatial_scalar<uint16_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, const DFTKernelContext&);
template void process_spatial_scalar<float>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, const DFTKernelContext&);
template void process_temporal_scalar<uint8_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, int, const DFTKernelContext&);
template void process_temporal_scalar<uint16_t>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, int, const DFTKernelContext&);
template void process_temporal_scalar<float>(unsigned int, int, DFTPlaneBytes, DFTMutablePlaneBytes, int, const DFTKernelContext&);
