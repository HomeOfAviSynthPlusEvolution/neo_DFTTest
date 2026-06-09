#pragma once

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <random>
#include <cstdlib>
#include <cstring>
#include <string.h>

#include "fft/fft_backend.hpp"
#include "memory/aligned_buffer.hpp"

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define VS_RESTRICT restrict
#elif defined(__cplusplus) || defined(_MSC_VER)
#define VS_RESTRICT __restrict
#else
#define VS_RESTRICT
#endif

#define EXTRA(a,b) (((a) % (b)) ? ((b) - ((a) % (b))) : 0)

struct DFTKernelContext;
struct DFTTestData;

struct DFTPlaneBytes {
    const unsigned char* data {nullptr};
    int stride_bytes {0};
};

struct DFTMutablePlaneBytes {
    unsigned char* data {nullptr};
    int stride_bytes {0};
};

struct DFTMutableFloatSpan {
    float* data {nullptr};
    int size {0};
};

struct DFTConstFloatSpan {
    const float* data {nullptr};
    int size {0};
};

struct DFTConstFloatPlane {
    const float* data {nullptr};
    int width {0};
    int height {0};
    int stride_elements {0};
};

struct DFTMutableBytePlane {
    uint8_t* data {nullptr};
    int width {0};
    int height {0};
    int stride_elements {0};
};

template<typename T>
struct DFTConstSampleBlock {
    const T* data {nullptr};
    int stride_elements {0};
    int size {0};
};

struct DFTDitherContext {
    float multiplier {1.0f};
    int peak {1};
    int mode {0};
    std::mt19937* rng {nullptr};
    DFTMutableFloatSpan buffer;
};

struct DFTFilterInput {
    DFTMutableFloatSpan coefficients;
    DFTConstFloatSpan sigmas;
    DFTConstFloatSpan pmins;
    DFTConstFloatSpan pmaxs;
    DFTConstFloatSpan sigmas2;
};

enum class DFTFilterKind {
    wiener,
    hard_threshold,
    multiplier,
    range_multiplier,
    range_wiener,
    wiener_power,
    wiener_sqrt,
};

struct DFTFilterPlan {
    DFTFilterKind kind {DFTFilterKind::wiener};
    bool custom_f0_beta {false};
};

using DFTCopyPadFunction = void (*)(int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) noexcept;
using DFTFilterCoefficientsFunction = void (*)(DFTFilterInput input);
using DFTProcessSpatialFunction = void (*)(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context);
using DFTProcessTemporalFunction = void (*)(unsigned int thread_id, int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, int temporal_position, const DFTKernelContext& context);

struct DFTFftState {
    neo_dfttest::fft::Backend* backend {nullptr};
    neo_dfttest::fft::Plan forward;
    neo_dfttest::fft::Plan inverse;
};

struct DFTClipFormat {
    int num_planes {0};
    int bytes_per_sample {0};
    int bits_per_sample {0};
    bool integer {true};
    int width {0};
    int height {0};
    int subsampling_w {0};
    int subsampling_h {0};
};

struct DFTBlockSettings {
    int spatial_size {16};
    int spatial_overlap {12};
    int temporal_size {3};
    int temporal_overlap {0};
    int spatial_window {0};
    int temporal_window {7};
    float spatial_beta {2.5f};
    float temporal_beta {2.5f};
    float f0_beta {1.0f};
    bool zero_mean {true};
    int dither_mode {0};
    int worker_threads {4};
};

struct DFTPlaneGeometry {
    std::array<int, 4> process {};
    std::array<int, 4> width {};
    std::array<int, 4> height {};
    std::array<int, 4> pad_width {};
    std::array<int, 4> pad_height {};
    std::array<int, 4> pad_stride {};
    std::array<int, 4> pad_block_size {};
    std::array<int, 4> e_stride {};
    std::array<int, 4> e_height {};
};

struct DFTSampleScale {
    float divisor {1.0f};
    float multiplier {1.0f};
    int peak {1};
};

struct DFTDerivedGeometry {
    int block_area {0};
    int block_volume {0};
    int complex_count {0};
    int coefficient_count {0};
    int transform_type {0};
    int spatial_center {0};
    int step {1};
    bool custom_f0_beta {false};
};

inline int dft_scratch_real_stride(const DFTDerivedGeometry& derived) noexcept {
    return ((derived.block_volume + 7) | 15) + 1;
}

inline int dft_scratch_complex_stride(const DFTDerivedGeometry& derived) noexcept {
    return ((derived.complex_count + 7) | 15) + 1;
}

constexpr int kMaxDftFftBatchSlots = 16;
constexpr int kDftFftPipelineSlots = 2;

struct DFTBlockBatch {
    std::array<int, kMaxDftFftBatchSlots> x_offsets {};
    int count {0};
};

inline int dft_fft_batch_capacity(const DFTBlockSettings& block, int backend_max_batch_size) noexcept {
    const int backend_limit = std::max(1, backend_max_batch_size);
    return std::min({std::max(1, block.worker_threads), kMaxDftFftBatchSlots, backend_limit});
}

inline int dft_fft_scratch_slots(const DFTBlockSettings& block, int backend_max_batch_size) noexcept {
    return dft_fft_batch_capacity(block, backend_max_batch_size) * kDftFftPipelineSlots;
}

inline float* dft_real_batch_data(float* base, int stride, int index) noexcept {
    return base + stride * index;
}

inline neo_dfttest::fft::Complex* dft_complex_batch_data(neo_dfttest::fft::Complex* base, int stride, int index) noexcept {
    return base + stride * index;
}

struct DFTCoefficientTables {
    neo_dfttest::AlignedBuffer<float> window;
    neo_dfttest::AlignedBuffer<float> sigmas;
    neo_dfttest::AlignedBuffer<float> sigmas2;
    neo_dfttest::AlignedBuffer<float> pmins;
    neo_dfttest::AlignedBuffer<float> pmaxs;
    neo_dfttest::AlignedBuffer<neo_dfttest::fft::Complex> window_dft;
};

struct DFTThreadScratchSlot {
    DFTThreadScratchSlot() = default;
    DFTThreadScratchSlot(const DFTThreadScratchSlot&) = delete;
    DFTThreadScratchSlot& operator=(const DFTThreadScratchSlot&) = delete;
    DFTThreadScratchSlot(DFTThreadScratchSlot&&) noexcept = default;
    DFTThreadScratchSlot& operator=(DFTThreadScratchSlot&&) noexcept = default;

    neo_dfttest::AlignedBuffer<float> dither_buffer;
    std::unique_ptr<std::mt19937> rng;
    neo_dfttest::AlignedBuffer<float> ebuff;
    neo_dfttest::AlignedBuffer<float> dftr;
    neo_dfttest::AlignedBuffer<neo_dfttest::fft::Complex> dftc;
    neo_dfttest::AlignedBuffer<neo_dfttest::fft::Complex> dftc2;
    int fft_scratch_slots {0};
};

struct DFTThreadScratch {
    std::vector<DFTThreadScratchSlot> slots;
};

struct DFTKernelDispatch {
    DFTCopyPadFunction copy_pad {nullptr};
    DFTFilterCoefficientsFunction filter_coefficients {nullptr};
    DFTFilterPlan filter_plan;
    DFTProcessSpatialFunction process_spatial {nullptr};
    DFTProcessTemporalFunction process_temporal {nullptr};
};

struct DFTKernelContext {
    const DFTFftState& fft;
    const DFTClipFormat& format;
    const DFTBlockSettings& block;
    const DFTPlaneGeometry& planes;
    const DFTSampleScale& sample;
    const DFTDerivedGeometry& derived;
    const DFTCoefficientTables& coefficients;
    DFTThreadScratch& scratch;
    DFTFilterCoefficientsFunction filter_coefficients {nullptr};
    DFTFilterPlan filter_plan;
};

struct DFTTestData {
    DFTFftState fft;
    DFTClipFormat format;
    DFTBlockSettings block;
    DFTPlaneGeometry planes;
    DFTSampleScale sample;
    DFTDerivedGeometry derived;
    DFTCoefficientTables coefficients;
    mutable DFTThreadScratch scratch;
    DFTKernelDispatch kernels;
};

inline DFTKernelContext make_kernel_context(const DFTTestData& state) noexcept {
    return DFTKernelContext{
        state.fft,
        state.format,
        state.block,
        state.planes,
        state.sample,
        state.derived,
        state.coefficients,
        state.scratch,
        state.kernels.filter_coefficients,
        state.kernels.filter_plan
    };
}

inline float* complex_float_data(neo_dfttest::fft::Complex* data) noexcept {
    return reinterpret_cast<float*>(data);
}

inline const float* complex_float_data(const neo_dfttest::fft::Complex* data) noexcept {
    return reinterpret_cast<const float*>(data);
}

inline DFTFilterInput make_filter_input(
    neo_dfttest::fft::Complex* coefficients,
    const DFTKernelContext& context
) noexcept {
    const float* pmins = context.filter_plan.custom_f0_beta
        ? &context.block.f0_beta
        : context.coefficients.pmins.data();
    const int pmins_size = context.filter_plan.custom_f0_beta ? 1 : context.derived.coefficient_count;

    return DFTFilterInput{
        DFTMutableFloatSpan{complex_float_data(coefficients), context.derived.coefficient_count},
        DFTConstFloatSpan{context.coefficients.sigmas.data(), context.derived.coefficient_count},
        DFTConstFloatSpan{pmins, pmins_size},
        DFTConstFloatSpan{context.coefficients.pmaxs.data(), context.derived.coefficient_count},
        DFTConstFloatSpan{context.coefficients.sigmas2.data(), context.derived.coefficient_count}
    };
}

struct NPInfo {
    int fn, b, y, x;
};

DFTKernelDispatch selectFunctions(const unsigned ftype, const unsigned opt, const DFTClipFormat& format, const DFTBlockSettings& block) noexcept;
void createWindow(DFTMutableFloatSpan window, const int tmode, const int smode, const DFTBlockSettings& block, const DFTDerivedGeometry& derived) noexcept;
std::vector<float> parseSigmaLocation(const std::vector<float>& s, const float sigma, const float pfact);
float interp(const float pf, DFTConstFloatSpan points) noexcept;
float getSVal(const int pos, const int len, DFTConstFloatSpan points, float & pf) noexcept;
void remove_mean_scalar(DFTMutableFloatSpan coefficients, DFTConstFloatSpan reference, DFTMutableFloatSpan removed) noexcept;
template<typename T>
void load_windowed_block_scalar(DFTConstSampleBlock<T> source, DFTConstFloatSpan window, DFTMutableFloatSpan destination, float divisor) noexcept;
void dither_u8_scalar(DFTConstFloatPlane source, DFTMutableBytePlane destination, DFTDitherContext context) noexcept;
