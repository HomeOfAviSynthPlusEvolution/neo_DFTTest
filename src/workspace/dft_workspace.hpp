#pragma once

#include "dft_common.h"

#include <cstddef>
#include <cstring>
#include <random>

namespace neo_dfttest {

struct DFTBatchScratchBuffers {
    float* real {nullptr};
    fft::Complex* coefficients {nullptr};
    fft::Complex* removed_mean {nullptr};
};

struct DFTThreadWorkspaceView {
    float* accumulation {nullptr};
    DFTBatchScratchBuffers batch;
    std::mt19937* dither_rng {nullptr};
    float* dither_buffer {nullptr};
};

inline DFTThreadWorkspaceView dft_thread_workspace(DFTThreadScratchSlot& slot) noexcept {
    return DFTThreadWorkspaceView{
        slot.ebuff.data(),
        DFTBatchScratchBuffers{slot.dftr.data(), slot.dftc.data(), slot.dftc2.data()},
        slot.rng.get(),
        slot.dither_buffer.data(),
    };
}

inline DFTThreadWorkspaceView dft_thread_workspace(const DFTKernelContext& context, unsigned int thread_id) noexcept {
    return dft_thread_workspace(context.scratch.slots[thread_id]);
}

inline void clear_accumulation(DFTThreadWorkspaceView workspace, std::size_t elements) noexcept {
    std::memset(workspace.accumulation, 0, elements * sizeof(float));
}

} // namespace neo_dfttest
