#ifndef DFT_COMMON_H
#define DFT_COMMON_H

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

#ifdef HAS_EXECUTION
  #include <execution>
#endif

#ifndef __cpp_lib_execution
  #undef ENABLE_PAR
#endif

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

struct DFTTestData;

using DFTCopyPadFunction = void (*)(int plane, const unsigned char *, int, unsigned char *, const DFTTestData *) noexcept;
using DFTFilterCoefficientsFunction = void (*)(float *, const float *, const int, const float *, const float *, const float *);
using DFTProcessSpatialFunction = void (*)(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const DFTTestData *);
using DFTProcessTemporalFunction = void (*)(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const int, const DFTTestData *);

struct DFTFftState {
    neo_dfttest::fft::Backend* backend {nullptr};
    neo_dfttest::fft::Plan forward {nullptr};
    neo_dfttest::fft::Plan inverse {nullptr};
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
    std::array<int, 4> e_batch_size {};
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
};

struct DFTThreadScratch {
    std::vector<DFTThreadScratchSlot> slots;
};

struct DFTKernelDispatch {
    DFTCopyPadFunction copy_pad {nullptr};
    DFTFilterCoefficientsFunction filter_coefficients {nullptr};
    DFTProcessSpatialFunction process_spatial {nullptr};
    DFTProcessTemporalFunction process_temporal {nullptr};
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

struct NPInfo {
    int fn, b, y, x;
};

DFTKernelDispatch selectFunctions(const unsigned ftype, const unsigned opt, const DFTTestData& d) noexcept;
void createWindow(float * VS_RESTRICT hw, const int tmode, const int smode, const DFTTestData * d) noexcept;
float * parseSigmaLocation(const std::vector<float> s, int & poscnt, const float sigma, const float pfact);
float interp(const float pf, const float * pv, const int cnt) noexcept;
float getSVal(const int pos, const int len, const float * pv, const int cnt, float & pf) noexcept;
void removeMean_c(float * VS_RESTRICT dftc, const float * dftgc, const int ccnt, float * VS_RESTRICT dftc2) noexcept;
template<typename T>
void proc0_c(const T * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept;
void dither_c(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode, std::mt19937 &rng, float *dither_buff) noexcept;

#endif
