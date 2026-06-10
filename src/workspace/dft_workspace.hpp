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

struct DftWorkspaceLease {
    unsigned int slot_id {0};
    DFTThreadWorkspaceView host;

    [[nodiscard]] DFTThreadWorkspaceView host_view() const noexcept {
        return host;
    }
};

inline DFTThreadWorkspaceView dft_thread_workspace(DFTThreadScratchSlot& slot) noexcept {
    return DFTThreadWorkspaceView{
        slot.ebuff.data(),
        DFTBatchScratchBuffers{slot.dftr.data(), slot.dftc.data(), slot.dftc2.data()},
        slot.rng.get(),
        slot.dither_buffer.data(),
    };
}

inline DftWorkspaceLease dft_workspace_lease(DFTThreadScratch& scratch, const unsigned int slot_id) noexcept {
    return DftWorkspaceLease{
        slot_id,
        dft_thread_workspace(scratch.slots[slot_id])
    };
}

inline void clear_accumulation(DFTThreadWorkspaceView workspace, std::size_t elements) noexcept {
    std::memset(workspace.accumulation, 0, elements * sizeof(float));
}

} // namespace neo_dfttest
