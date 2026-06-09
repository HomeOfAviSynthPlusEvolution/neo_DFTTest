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

class DftBatchPipeline {
public:
    DftBatchPipeline(
        const DFTKernelContext& context,
        DFTBatchScratchBuffers scratch,
        const int width,
        const int eheight
    ) noexcept
        : context_(context),
          scratch_(scratch),
          width_(width),
          eheight_(eheight),
          real_stride_(dft_scratch_real_stride(context.derived)),
          complex_stride_(dft_scratch_complex_stride(context.derived)),
          batch_capacity_(dft_fft_batch_capacity(context.batch_policy)) {}

    template<typename LoadBlock, typename CompleteBatch>
    void run(LoadBlock&& load_block, CompleteBatch&& complete_batch) {
        std::optional<PendingBatch> pending;
        int active_slot = 0;

        for (int y = 0; y < eheight_; y += context_.derived.step) {
            for (int x = 0; x <= width_ - context_.block.spatial_size;) {
                float* active_real = real_slot(active_slot);
                auto* active_coefficients = complex_slot(scratch_.coefficients, active_slot);
                DFTBlockBatch batch;

                for (; batch.count < batch_capacity_ && x <= width_ - context_.block.spatial_size; ++batch.count, x += context_.derived.step) {
                    float* real_block = dft_real_batch_data(active_real, real_stride_, batch.count);
                    batch.x_offsets[static_cast<std::size_t>(batch.count)] = x;
                    load_block(y, x, real_block);
                }

                auto forward = submit_forward(active_real, active_coefficients, batch.count);

                if (pending) {
                    consume_pending(*pending, complete_batch);
                }
                pending.emplace(PendingBatch{batch, y, active_slot, std::move(forward)});
                active_slot = 1 - active_slot;
            }
        }

        if (pending) {
            consume_pending(*pending, complete_batch);
        }
    }

private:
    struct PendingBatch {
        DFTBlockBatch batch;
        int y {0};
        int slot {0};
        fft::Completion forward;
    };

    float* real_slot(const int slot) const noexcept {
        return dft_real_batch_data(scratch_.real, real_stride_, slot * batch_capacity_);
    }

    fft::Complex* complex_slot(fft::Complex* base, const int slot) const noexcept {
        return dft_complex_batch_data(base, complex_stride_, slot * batch_capacity_);
    }

    fft::Completion submit_forward(float* active_real, fft::Complex* active_coefficients, const int count) const {
        return context_.fft.backend->submit_r2c(
            context_.fft.forward,
            fft::R2CBatch{
                fft::RealBatchView{active_real, real_stride_, fft::MemoryDomain::host},
                fft::ComplexBatchView{active_coefficients, complex_stride_, fft::MemoryDomain::host},
                count,
            }
        );
    }

    template<typename CompleteBatch>
    void consume_pending(PendingBatch& pending_batch, CompleteBatch& complete_batch) const {
        pending_batch.forward.wait();
        complete_batch(DFTCompletedBatch{
            pending_batch.batch,
            pending_batch.y,
            real_slot(pending_batch.slot),
            complex_slot(scratch_.coefficients, pending_batch.slot),
            complex_slot(scratch_.removed_mean, pending_batch.slot),
            real_stride_,
            complex_stride_
        });
    }

    const DFTKernelContext& context_;
    DFTBatchScratchBuffers scratch_;
    int width_ {0};
    int eheight_ {0};
    int real_stride_ {0};
    int complex_stride_ {0};
    int batch_capacity_ {0};
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
    DftBatchPipeline{context, scratch, width, eheight}.run(
        std::forward<LoadBlock>(load_block),
        std::forward<CompleteBatch>(complete_batch)
    );
}

} // namespace neo_dfttest
