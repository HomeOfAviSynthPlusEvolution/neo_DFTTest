#pragma once

#include "workspace/dft_workspace.hpp"

#include <optional>
#include <utility>

namespace neo_dfttest {

struct DFTCompletedBatch {
    const DFTBlockBatch& batch;
    int y {0};
    float* real {nullptr};
    fft::Complex* coefficients {nullptr};
    fft::Complex* removed_mean {nullptr};
    int real_stride {0};
    int complex_stride {0};
};

template<typename LoadBlock, typename CompleteBatch>
void run_dft_batch_pipeline(
    const DFTKernelContext& context,
    DFTBatchScratchBuffers scratch,
    const int width,
    const int eheight,
    LoadBlock&& load_block,
    CompleteBatch&& complete_batch
) {
    const int real_stride = dft_scratch_real_stride(context.derived);
    const int complex_stride = dft_scratch_complex_stride(context.derived);
    const int batch_capacity = dft_fft_batch_capacity(context.block, context.fft.backend->capabilities().max_batch_size);

    struct PendingBatch {
        DFTBlockBatch batch;
        int y {0};
        int slot {0};
        fft::Completion forward;
    };
    std::optional<PendingBatch> pending;

    auto real_slot = [&](const int slot) noexcept {
        return dft_real_batch_data(scratch.real, real_stride, slot * batch_capacity);
    };
    auto complex_slot = [&](fft::Complex* base, const int slot) noexcept {
        return dft_complex_batch_data(base, complex_stride, slot * batch_capacity);
    };
    auto consume_pending = [&](PendingBatch& pending_batch) {
        pending_batch.forward.wait();
        complete_batch(DFTCompletedBatch{
            pending_batch.batch,
            pending_batch.y,
            real_slot(pending_batch.slot),
            complex_slot(scratch.coefficients, pending_batch.slot),
            complex_slot(scratch.removed_mean, pending_batch.slot),
            real_stride,
            complex_stride
        });
    };

    int active_slot = 0;
    for (int y = 0; y < eheight; y += context.derived.step) {
        for (int x = 0; x <= width - context.block.spatial_size;) {
            float* active_real = real_slot(active_slot);
            auto* active_coefficients = complex_slot(scratch.coefficients, active_slot);
            DFTBlockBatch batch;

            for (; batch.count < batch_capacity && x <= width - context.block.spatial_size; ++batch.count, x += context.derived.step) {
                float* real_block = dft_real_batch_data(active_real, real_stride, batch.count);
                batch.x_offsets[static_cast<std::size_t>(batch.count)] = x;
                load_block(y, x, real_block);
            }

            auto forward = context.fft.backend->submit_r2c(
                context.fft.forward,
                fft::R2CBatch{
                    fft::RealBatchView{active_real, real_stride, fft::MemoryDomain::host},
                    fft::ComplexBatchView{active_coefficients, complex_stride, fft::MemoryDomain::host},
                    batch.count,
                }
            );

            if (pending) {
                consume_pending(*pending);
            }
            pending.emplace(PendingBatch{batch, y, active_slot, std::move(forward)});
            active_slot = 1 - active_slot;
        }
    }

    if (pending) {
        consume_pending(*pending);
    }
}

} // namespace neo_dfttest
