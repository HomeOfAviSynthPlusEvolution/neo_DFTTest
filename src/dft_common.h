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

template<typename T>
struct DFTConstSamplePlane {
    const T* data {nullptr};
    int stride_elements {0};
};

template<typename T>
struct DFTMutableSamplePlane {
    T* data {nullptr};
    int stride_elements {0};
};

template<typename T>
inline DFTConstSamplePlane<T> dft_const_sample_plane(DFTPlaneBytes plane) noexcept {
    return DFTConstSamplePlane<T>{
        reinterpret_cast<const T*>(plane.data),
        plane.stride_bytes / static_cast<int>(sizeof(T))
    };
}

template<typename T>
inline DFTMutableSamplePlane<T> dft_mutable_sample_plane(DFTMutablePlaneBytes plane) noexcept {
    return DFTMutableSamplePlane<T>{
        reinterpret_cast<T*>(plane.data),
        plane.stride_bytes / static_cast<int>(sizeof(T))
    };
}

struct DFTMutableFloatSpan {
    float* data {nullptr};
    int size {0};
};

struct DFTConstFloatSpan {
    const float* data {nullptr};
    int size {0};
};

struct DFTConstComplexSpan {
    const neo_dfttest::fft::Complex* data {nullptr};
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
constexpr int kMaxDftTemporalFrames = 15;

struct DFTBatchPolicy {
    int cpu_worker_threads {1};
    int fft_max_batch_size {1};
    int executor_preferred_batch_size {1};
};

struct DFTBlockJob {
    int plane {-1};
    int y {0};
    int x {0};
    int temporal_position {0};
};

struct DFTBlockBatch {
    std::array<DFTBlockJob, kMaxDftFftBatchSlots> jobs {};
    int count {0};
};

inline const DFTBlockJob& dft_block_job(const DFTBlockBatch& batch, const int index) noexcept {
    return batch.jobs[static_cast<std::size_t>(index)];
}

inline DFTBatchPolicy make_cpu_dft_batch_policy(const DFTBlockSettings& block, const neo_dfttest::fft::BackendCapabilities& backend) noexcept {
    const int worker_threads = std::max(1, block.worker_threads);
    return DFTBatchPolicy{
        worker_threads,
        std::max(1, backend.max_batch_size),
        worker_threads
    };
}

inline int dft_fft_batch_capacity(const DFTBatchPolicy& policy) noexcept {
    const int executor_limit = std::max(1, policy.executor_preferred_batch_size);
    const int backend_limit = std::max(1, policy.fft_max_batch_size);
    return std::min({executor_limit, kMaxDftFftBatchSlots, backend_limit});
}

inline int dft_fft_scratch_slots(const DFTBatchPolicy& policy) noexcept {
    return dft_fft_batch_capacity(policy) * kDftFftPipelineSlots;
}

inline float* dft_real_batch_data(float* base, int stride, int index) noexcept {
    return base + stride * index;
}

inline neo_dfttest::fft::Complex* dft_complex_batch_data(neo_dfttest::fft::Complex* base, int stride, int index) noexcept {
    return base + stride * index;
}

struct DFTMutableRealBatchView {
    float* data {nullptr};
    int stride_elements {0};
    int count {0};
    neo_dfttest::fft::MemoryDomain domain {neo_dfttest::fft::MemoryDomain::host};

    [[nodiscard]] float* block(const int index) const noexcept {
        return dft_real_batch_data(data, stride_elements, index);
    }

    [[nodiscard]] neo_dfttest::fft::RealBatchView fft_view() const noexcept {
        return neo_dfttest::fft::RealBatchView{data, stride_elements, domain};
    }
};

struct DFTMutableComplexBatchView {
    neo_dfttest::fft::Complex* data {nullptr};
    int stride_elements {0};
    int count {0};
    neo_dfttest::fft::MemoryDomain domain {neo_dfttest::fft::MemoryDomain::host};

    [[nodiscard]] neo_dfttest::fft::Complex* block(const int index) const noexcept {
        return dft_complex_batch_data(data, stride_elements, index);
    }

    [[nodiscard]] neo_dfttest::fft::ComplexBatchView fft_view() const noexcept {
        return neo_dfttest::fft::ComplexBatchView{data, stride_elements, domain};
    }
};

inline DFTMutableRealBatchView dft_host_real_batch_view(float* data, const int stride, const int count) noexcept {
    return DFTMutableRealBatchView{data, stride, count, neo_dfttest::fft::MemoryDomain::host};
}

inline DFTMutableComplexBatchView dft_host_complex_batch_view(
    neo_dfttest::fft::Complex* data,
    const int stride,
    const int count
) noexcept {
    return DFTMutableComplexBatchView{data, stride, count, neo_dfttest::fft::MemoryDomain::host};
}

struct DFTCoefficientTables {
    neo_dfttest::AlignedBuffer<float> window;
    neo_dfttest::AlignedBuffer<float> sigmas;
    neo_dfttest::AlignedBuffer<float> sigmas2;
    neo_dfttest::AlignedBuffer<float> pmins;
    neo_dfttest::AlignedBuffer<float> pmaxs;
    neo_dfttest::AlignedBuffer<neo_dfttest::fft::Complex> window_dft;
};

struct DFTCoefficientView {
    DFTConstFloatSpan window;
    DFTConstFloatSpan sigmas;
    DFTConstFloatSpan sigmas2;
    DFTConstFloatSpan pmins;
    DFTConstFloatSpan pmaxs;
    DFTConstComplexSpan window_dft;
};

inline DFTCoefficientView dft_host_coefficient_view(
    const DFTCoefficientTables& coefficients,
    const DFTDerivedGeometry& derived
) noexcept {
    return DFTCoefficientView{
        DFTConstFloatSpan{coefficients.window.data(), derived.block_volume},
        DFTConstFloatSpan{coefficients.sigmas.data(), derived.coefficient_count},
        DFTConstFloatSpan{coefficients.sigmas2.data(), derived.coefficient_count},
        DFTConstFloatSpan{coefficients.pmins.data(), derived.coefficient_count},
        DFTConstFloatSpan{coefficients.pmaxs.data(), derived.coefficient_count},
        DFTConstComplexSpan{coefficients.window_dft.data(), derived.coefficient_count}
    };
}

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
    DFTFilterPlan filter_plan;
};

struct DFTKernelContext {
    const DFTFftState& fft;
    const DFTClipFormat& format;
    const DFTBlockSettings& block;
    const DFTPlaneGeometry& planes;
    const DFTSampleScale& sample;
    const DFTDerivedGeometry& derived;
    const DFTBatchPolicy& batch_policy;
    DFTCoefficientView coefficients;
    DFTThreadScratch& scratch;
    DFTFilterPlan filter_plan;
};

struct DFTTestData {
    DFTFftState fft;
    DFTClipFormat format;
    DFTBlockSettings block;
    DFTPlaneGeometry planes;
    DFTSampleScale sample;
    DFTDerivedGeometry derived;
    DFTBatchPolicy batch_policy;
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
        state.batch_policy,
        dft_host_coefficient_view(state.coefficients, state.derived),
        state.scratch,
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
    const DFTFilterPlan filter_plan,
    const DFTBlockSettings& block,
    const DFTDerivedGeometry& derived,
    DFTCoefficientView coefficient_tables
) noexcept {
    const float* pmins = filter_plan.custom_f0_beta
        ? &block.f0_beta
        : coefficient_tables.pmins.data;
    const int pmins_size = filter_plan.custom_f0_beta ? 1 : derived.coefficient_count;

    return DFTFilterInput{
        DFTMutableFloatSpan{complex_float_data(coefficients), derived.coefficient_count},
        DFTConstFloatSpan{coefficient_tables.sigmas.data, derived.coefficient_count},
        DFTConstFloatSpan{pmins, pmins_size},
        DFTConstFloatSpan{coefficient_tables.pmaxs.data, derived.coefficient_count},
        DFTConstFloatSpan{coefficient_tables.sigmas2.data, derived.coefficient_count}
    };
}

inline DFTFilterInput make_filter_input(
    neo_dfttest::fft::Complex* coefficients,
    const DFTKernelContext& context
) noexcept {
    return make_filter_input(
        coefficients,
        context.filter_plan,
        context.block,
        context.derived,
        context.coefficients
    );
}

struct DFTFilterBatchOperation {
    DFTFilterPlan plan;
    DFTMutableComplexBatchView coefficients;
    const DFTBlockSettings* block {nullptr};
    const DFTDerivedGeometry* derived {nullptr};
    DFTCoefficientView coefficient_tables;

    [[nodiscard]] DFTFilterInput input(const int index) const noexcept {
        return make_filter_input(
            coefficients.block(index),
            plan,
            *block,
            *derived,
            coefficient_tables
        );
    }
};

inline DFTFilterBatchOperation dft_filter_batch_operation(
    DFTMutableComplexBatchView coefficients,
    const DFTKernelContext& context
) noexcept {
    return DFTFilterBatchOperation{
        context.filter_plan,
        coefficients,
        &context.block,
        &context.derived,
        context.coefficients
    };
}

struct DFTFilterStageOperation {
    DFTFilterBatchOperation batch;
    neo_dfttest::fft::MemoryDomain domain {neo_dfttest::fft::MemoryDomain::host};

    [[nodiscard]] DFTFilterInput input(const int index) const noexcept {
        return batch.input(index);
    }
};

inline DFTFilterStageOperation dft_filter_stage_operation(
    DFTMutableComplexBatchView coefficients,
    const DFTKernelContext& context
) noexcept {
    return DFTFilterStageOperation{
        dft_filter_batch_operation(coefficients, context),
        coefficients.domain
    };
}

struct DFTInverseAccumulationOperation {
    const float* inverse {nullptr};
    const float* window {nullptr};
    float* output {nullptr};
    int output_stride {0};
    const DFTKernelContext* context {nullptr};
    neo_dfttest::fft::MemoryDomain domain {neo_dfttest::fft::MemoryDomain::host};
};

inline DFTInverseAccumulationOperation dft_inverse_accumulation_operation(
    const float* inverse,
    const float* window,
    float* output,
    const DFTKernelContext& context,
    const int output_stride,
    const neo_dfttest::fft::MemoryDomain domain = neo_dfttest::fft::MemoryDomain::host
) noexcept {
    return DFTInverseAccumulationOperation{
        inverse,
        window,
        output,
        output_stride,
        &context,
        domain
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
