#pragma once

#include "core_batch_pipeline.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace neo_dfttest {

template<typename T, typename Stages>
void run_cpu_spatial_dft_pipeline(
    DFTThreadWorkspaceView workspace,
    int plane,
    DFTPlaneBytes src,
    DFTMutablePlaneBytes dst,
    const DFTKernelContext& context,
    Stages stages
) {
    float* ebuff = workspace.accumulation;

    const int width = context.planes.pad_width[plane];
    const int height = context.planes.pad_height[plane];
    const int eheight = context.planes.e_height[plane];
    const auto source = dft_const_sample_plane<T>(src);
    const int srcStride = source.stride_elements;
    const int ebpStride = context.planes.e_stride[plane];
    const T* srcp = source.data;

    clear_accumulation(workspace, static_cast<std::size_t>(ebpStride) * height);

    run_dft_batch_pipeline(
        context,
        workspace.batch,
        width,
        eheight,
        [&](int y, int x, float* dftr) noexcept {
            stages.template load_windowed_block<T>(
                DFTConstSampleBlock<T>{srcp + y * srcStride + x, srcStride, context.block.spatial_size},
                DFTConstFloatSpan{context.coefficients.window.data, context.derived.block_area},
                DFTMutableFloatSpan{dftr, context.derived.block_area},
                context.sample.divisor
            );
        },
        [&](DFTCompletedBatch ready) {
            const auto filter_operation = dft_filter_stage_operation(ready.coefficients, context);
            for (int index = 0; index < ready.batch.count; ++index) {
                auto* dftc = ready.coefficients.block(index);
                auto* dftc2 = ready.removed_mean.block(index);
                if (context.block.zero_mean) {
                    stages.remove_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(context.coefficients.window_dft.data), context.derived.coefficient_count},
                        DFTMutableFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );
                }

                stages.apply_filter(filter_operation, index);

                if (context.block.zero_mean) {
                    stages.add_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );
                }
            }

            DftFftOperations{context}.submit_inverse_and_wait(ready.coefficients, ready.real);

            float* output_row = ebuff + ready.y * ebpStride;
            for (int index = 0; index < ready.batch.count; ++index) {
                float* dftr = ready.real.block(index);
                const int block_x = dft_block_job(ready.batch, index).x;
                stages.accumulate_inverse_block(
                    dft_inverse_accumulation_operation(
                        dftr,
                        context.coefficients.window.data,
                        output_row + block_x,
                        context,
                        ebpStride,
                        ready.real.domain
                    )
                );
            }
        }
    );

    stages.template write_output<T>(workspace, dst, plane, width, height, context);
}

template<typename T, typename Stages>
void run_cpu_temporal_dft_pipeline(
    DFTThreadWorkspaceView workspace,
    int plane,
    DFTPlaneBytes src,
    DFTMutablePlaneBytes dst,
    int temporal_position,
    const DFTKernelContext& context,
    Stages stages
) {
    float* ebuff = workspace.accumulation;

    const int width = context.planes.pad_width[plane];
    const int height = context.planes.pad_height[plane];
    const int eheight = context.planes.e_height[plane];
    const auto first_source = dft_const_sample_plane<T>(src);
    const int srcStride = first_source.stride_elements;
    const int ebpStride = context.planes.e_stride[plane];

    std::array<const T*, kMaxDftTemporalFrames> srcp {};
    for (int index = 0; index < context.block.temporal_size; ++index) {
        const auto source = dft_const_sample_plane<T>(
            DFTPlaneBytes{src.data + context.planes.pad_block_size[plane] * index, src.stride_bytes}
        );
        srcp[static_cast<std::size_t>(index)] = source.data;
    }

    clear_accumulation(workspace, static_cast<std::size_t>(ebpStride) * height);

    run_dft_batch_pipeline(
        context,
        workspace.batch,
        width,
        eheight,
        [&](int y, int x, float* dftr) noexcept {
            for (int z = 0; z < context.block.temporal_size; ++z) {
                stages.template load_windowed_block<T>(
                    DFTConstSampleBlock<T>{srcp[static_cast<std::size_t>(z)] + y * srcStride + x, srcStride, context.block.spatial_size},
                    DFTConstFloatSpan{context.coefficients.window.data + context.derived.block_area * z, context.derived.block_area},
                    DFTMutableFloatSpan{dftr + context.derived.block_area * z, context.derived.block_area},
                    context.sample.divisor
                );
            }
        },
        [&](DFTCompletedBatch ready) {
            const auto filter_operation = dft_filter_stage_operation(ready.coefficients, context);
            for (int index = 0; index < ready.batch.count; ++index) {
                auto* dftc = ready.coefficients.block(index);
                auto* dftc2 = ready.removed_mean.block(index);
                if (context.block.zero_mean) {
                    stages.remove_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(context.coefficients.window_dft.data), context.derived.coefficient_count},
                        DFTMutableFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );
                }

                stages.apply_filter(filter_operation, index);

                if (context.block.zero_mean) {
                    stages.add_mean(
                        DFTMutableFloatSpan{complex_float_data(dftc), context.derived.coefficient_count},
                        DFTConstFloatSpan{complex_float_data(dftc2), context.derived.coefficient_count}
                    );
                }
            }

            DftFftOperations{context}.submit_inverse_and_wait(ready.coefficients, ready.real);

            for (int index = 0; index < ready.batch.count; ++index) {
                float* dftr = ready.real.block(index);
                const int block_x = dft_block_job(ready.batch, index).x;
                const int temporal_offset = temporal_position * context.derived.block_area;
                stages.accumulate_inverse_block(
                    dft_inverse_accumulation_operation(
                        dftr + temporal_offset,
                        context.coefficients.window.data + temporal_offset,
                        ebuff + ready.y * ebpStride + block_x,
                        context,
                        ebpStride,
                        ready.real.domain
                    )
                );
            }
        }
    );

    stages.template write_output<T>(workspace, dst, plane, width, height, context);
}

} // namespace neo_dfttest
