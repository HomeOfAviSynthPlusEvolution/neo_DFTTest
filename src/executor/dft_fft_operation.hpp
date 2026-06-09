#pragma once

#include "dft_common.h"

namespace neo_dfttest {

class DftFftOperations {
public:
    explicit DftFftOperations(const DFTKernelContext& context) noexcept
        : context_(context) {}

    [[nodiscard]] fft::Completion submit_forward(
        DFTMutableRealBatchView input,
        DFTMutableComplexBatchView output
    ) const {
        return context_.fft.backend->submit_r2c(
            context_.fft.forward,
            fft::R2CBatch{input.fft_view(), output.fft_view(), input.count}
        );
    }

    [[nodiscard]] fft::Completion submit_inverse(
        DFTMutableComplexBatchView input,
        DFTMutableRealBatchView output
    ) const {
        return context_.fft.backend->submit_c2r(
            context_.fft.inverse,
            fft::C2RBatch{input.fft_view(), output.fft_view(), input.count}
        );
    }

    void submit_inverse_and_wait(
        DFTMutableComplexBatchView input,
        DFTMutableRealBatchView output
    ) const noexcept {
        submit_inverse(input, output).wait();
    }

private:
    const DFTKernelContext& context_;
};

} // namespace neo_dfttest
