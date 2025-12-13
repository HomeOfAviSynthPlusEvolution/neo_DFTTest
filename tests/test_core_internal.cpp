#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <chrono>

// Define necessary macros for the included files
#ifndef VS_TARGET_CPU_X86
#define VS_TARGET_CPU_X86
#endif

#include "core.h"
#include "vectorclass/vectormath_exp.h"

// Include source files in separate namespaces to access static functions
// and avoid naming conflicts.

namespace Ref {
    // Include C reference implementation
    // We need to be careful about relative paths.
    // Since we are writing to tests/test_core_internal.cpp, the src directory is ../src
    // However, the tool uses paths relative to workspace root.
    // The preprocessor #include looks relative to the file or in include paths.
    // We will assume the build system sets the include path correctly to workspace root
    // or we use relative paths from this file.
    #include "../src/core_C.cpp"
}

namespace SSE2 {
    #include "../src/core_SSE2.cpp"
}

namespace AVX2 {
    #include "../src/core_AVX2.cpp"
}
void dither_c(float const*, unsigned char*, int, int, int, int, float, int, int, MTRand&, float*) noexcept {};

// Helper class for memory alignment
template <typename T>
class AlignedBuffer {
public:
    AlignedBuffer(size_t size, size_t alignment = 64) : size_(size), alignment_(alignment) {
#ifdef _WIN32
        data_ = static_cast<T*>(_aligned_malloc(size * sizeof(T), alignment));
#else
        data_ = static_cast<T*>(aligned_alloc(alignment, size * sizeof(T)));
#endif
        if (!data_) throw std::bad_alloc();
        // Initialize with zeros
        std::memset(data_, 0, size * sizeof(T));
    }

    ~AlignedBuffer() {
#ifdef _WIN32
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    T* data() { return data_; }
    const T* data() const { return data_; }
    size_t size() const { return size_; }

private:
    T* data_;
    size_t size_;
    size_t alignment_;
};

// Helper to fill buffer with random data
template <typename T>
void fill_random(T* data, size_t size, double min_val, double max_val) {
    static std::mt19937 rng(static_cast<uint32_t>(std::chrono::system_clock::now().time_since_epoch().count()));
    if constexpr (std::is_integral_v<T>) {
        std::uniform_int_distribution<int> dist(static_cast<int>(min_val), static_cast<int>(max_val));
        std::generate(data, data + size, [&]() { return static_cast<T>(dist(rng)); });
    } else {
        std::uniform_real_distribution<float> dist(static_cast<float>(min_val), static_cast<float>(max_val));
        std::generate(data, data + size, [&]() { return static_cast<T>(dist(rng)); });
    }
}

// Helper to fill buffer with random data, duplicating values for complex number pairs
void fill_random_pairs(float* data, size_t count, double min_val, double max_val) {
    static std::mt19937 rng(static_cast<uint32_t>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<float> dist((float)min_val, (float)max_val);
    for (size_t i = 0; i < count; i += 2) {
        float val = dist(rng);
        data[i] = val;
        if (i + 1 < count) data[i + 1] = val;
    }
}

// --- proc0 Tests ---

template <typename T>
class Proc0Test : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup if needed
    }
};

using Proc0Types = ::testing::Types<uint8_t, uint16_t, float>;
TYPED_TEST_SUITE(Proc0Test, Proc0Types);

TYPED_TEST(Proc0Test, MatchesReference) {
    using T = TypeParam;
    constexpr int p0 = 128; // stride for s0
    constexpr int p1 = 64;  // stride for s1 and d, and width
    constexpr int h = 64;   // height
    constexpr float divisor = 0.5f;

    AlignedBuffer<T> s0(p0 * h);
    AlignedBuffer<float> s1(p1 * h);
    AlignedBuffer<float> d_ref(p1 * h);
    AlignedBuffer<float> d_sse2(p1 * h);
    AlignedBuffer<float> d_avx2(p1 * h);

    fill_random(s0.data(), s0.size(), 0, 100);
    fill_random(s1.data(), s1.size(), 0.0, 1.0);

    // Call implementations
    Ref::proc0(s0.data(), s1.data(), d_ref.data(), p0, p1, divisor);
    SSE2::proc0(s0.data(), s1.data(), d_sse2.data(), p0, p1, divisor);
    AVX2::proc0(s0.data(), s1.data(), d_avx2.data(), p0, p1, divisor);

    // Verify
    for (size_t i = 0; i < d_ref.size(); ++i) {
        ASSERT_NEAR(d_ref.data()[i], d_sse2.data()[i], 1e-5) << "SSE2 mismatch at index " << i;
        ASSERT_NEAR(d_ref.data()[i], d_avx2.data()[i], 1e-5) << "AVX2 mismatch at index " << i;
    }
}

// --- proc1 Tests ---

TEST(Proc1Test, MatchesReference) {
    constexpr int p0 = 64; // stride for s0, s1
    constexpr int p1 = 128; // stride for d
    constexpr int h = p0;

    AlignedBuffer<float> s0(p0 * h);
    AlignedBuffer<float> s1(p0 * h);
    AlignedBuffer<float> d_ref(p1 * h);
    AlignedBuffer<float> d_sse2(p1 * h);
    AlignedBuffer<float> d_avx2(p1 * h);

    fill_random(s0.data(), s0.size(), -10.0, 10.0);
    fill_random(s1.data(), s1.size(), 0.0, 1.0);
    fill_random(d_ref.data(), d_ref.size(), 0.0, 0.0); // Init output with zeros or distinct values

    // Copy initial state to other buffers
    std::copy(d_ref.data(), d_ref.data() + d_ref.size(), d_sse2.data());
    std::copy(d_ref.data(), d_ref.data() + d_ref.size(), d_avx2.data());

    Ref::proc1(s0.data(), s1.data(), d_ref.data(), p0, p1);
    SSE2::proc1(s0.data(), s1.data(), d_sse2.data(), p0, p1);
    AVX2::proc1(s0.data(), s1.data(), d_avx2.data(), p0, p1);

    for (size_t i = 0; i < d_ref.size(); ++i) {
        ASSERT_NEAR(d_ref.data()[i], d_sse2.data()[i], 1e-4) << "SSE2 mismatch at index " << i;
        ASSERT_NEAR(d_ref.data()[i], d_avx2.data()[i], 1e-4) << "AVX2 mismatch at index " << i;
    }
}

// --- removeMean Tests ---

TEST(RemoveMeanTest, MatchesReference) {
    constexpr int ccnt = 256;
    AlignedBuffer<float> dftc_ref(ccnt);
    AlignedBuffer<float> dftc_sse2(ccnt);
    AlignedBuffer<float> dftc_avx2(ccnt);

    AlignedBuffer<float> dftgc(ccnt);
    AlignedBuffer<float> dftc2_ref(ccnt);
    AlignedBuffer<float> dftc2_sse2(ccnt);
    AlignedBuffer<float> dftc2_avx2(ccnt);

    fill_random(dftc_ref.data(), ccnt, -100.0, 100.0);
    fill_random(dftgc.data(), ccnt, 1.0, 10.0); // Avoid division by zero

    std::copy(dftc_ref.data(), dftc_ref.data() + ccnt, dftc_sse2.data());
    std::copy(dftc_ref.data(), dftc_ref.data() + ccnt, dftc_avx2.data());

    Ref::removeMean(dftc_ref.data(), dftgc.data(), ccnt, dftc2_ref.data());
    SSE2::removeMean(dftc_sse2.data(), dftgc.data(), ccnt, dftc2_sse2.data());
    AVX2::removeMean(dftc_avx2.data(), dftgc.data(), ccnt, dftc2_avx2.data());

    for (int i = 0; i < ccnt; ++i) {
        ASSERT_NEAR(dftc_ref.data()[i], dftc_sse2.data()[i], 1e-4) << "dftc SSE2 mismatch at " << i;
        ASSERT_NEAR(dftc_ref.data()[i], dftc_avx2.data()[i], 1e-4) << "dftc AVX2 mismatch at " << i;

        ASSERT_NEAR(dftc2_ref.data()[i], dftc2_sse2.data()[i], 1e-4) << "dftc2 SSE2 mismatch at " << i;
        ASSERT_NEAR(dftc2_ref.data()[i], dftc2_avx2.data()[i], 1e-4) << "dftc2 AVX2 mismatch at " << i;
    }
}

// --- addMean Tests ---

TEST(AddMeanTest, MatchesReference) {
    constexpr int ccnt = 256;
    AlignedBuffer<float> dftc_ref(ccnt);
    AlignedBuffer<float> dftc_sse2(ccnt);
    AlignedBuffer<float> dftc_avx2(ccnt);

    AlignedBuffer<float> dftc2(ccnt);

    fill_random(dftc_ref.data(), ccnt, -100.0, 100.0);
    fill_random(dftc2.data(), ccnt, -50.0, 50.0);

    std::copy(dftc_ref.data(), dftc_ref.data() + ccnt, dftc_sse2.data());
    std::copy(dftc_ref.data(), dftc_ref.data() + ccnt, dftc_avx2.data());

    Ref::addMean(dftc_ref.data(), ccnt, dftc2.data());
    SSE2::addMean(dftc_sse2.data(), ccnt, dftc2.data());
    AVX2::addMean(dftc_avx2.data(), ccnt, dftc2.data());

    for (int i = 0; i < ccnt; ++i) {
        ASSERT_NEAR(dftc_ref.data()[i], dftc_sse2.data()[i], 1e-4) << "SSE2 mismatch at " << i;
        ASSERT_NEAR(dftc_ref.data()[i], dftc_avx2.data()[i], 1e-4) << "AVX2 mismatch at " << i;
    }
}

// --- filter Tests ---

template <typename T>
class FilterTest : public ::testing::Test {};

using FilterTypes = ::testing::Types<
    std::integral_constant<int, 0>,
    std::integral_constant<int, 1>,
    std::integral_constant<int, 2>,
    std::integral_constant<int, 3>,
    std::integral_constant<int, 4>,
    std::integral_constant<int, 5>,
    std::integral_constant<int, 6>
>;
TYPED_TEST_SUITE(FilterTest, FilterTypes);

TYPED_TEST(FilterTest, MatchesReference) {
    constexpr int type = TypeParam::value;
    constexpr int ccnt = 256;

    AlignedBuffer<float> dftc_ref(ccnt);
    AlignedBuffer<float> dftc_sse2(ccnt);
    AlignedBuffer<float> dftc_avx2(ccnt);

    AlignedBuffer<float> sigmas(ccnt);
    AlignedBuffer<float> sigmas2(ccnt);
    AlignedBuffer<float> pmin(ccnt);
    AlignedBuffer<float> pmax(ccnt);

    fill_random(dftc_ref.data(), ccnt, -10.0, 10.0);
    fill_random_pairs(sigmas.data(), ccnt, 0.1, 5.0);
    fill_random_pairs(sigmas2.data(), ccnt, 0.1, 5.0);
    fill_random_pairs(pmin.data(), ccnt, 0.0, 2.0);
    fill_random_pairs(pmax.data(), ccnt, 3.0, 10.0);

    // For type 5, pmin[0] is used as beta
    if constexpr (type == 5) {
         pmin.data()[0] = 2.0f; // beta
    }

    std::copy(dftc_ref.data(), dftc_ref.data() + ccnt, dftc_sse2.data());
    std::copy(dftc_ref.data(), dftc_ref.data() + ccnt, dftc_avx2.data());

    Ref::filter_c<type>(dftc_ref.data(), sigmas.data(), ccnt, pmin.data(), pmax.data(), sigmas2.data());
    SSE2::filter_sse2<type>(dftc_sse2.data(), sigmas.data(), ccnt, pmin.data(), pmax.data(), sigmas2.data());
    AVX2::filter_avx2<type>(dftc_avx2.data(), sigmas.data(), ccnt, pmin.data(), pmax.data(), sigmas2.data());

    for (int i = 0; i < ccnt; ++i) {
        ASSERT_NEAR(dftc_ref.data()[i], dftc_sse2.data()[i], 1e-4) << "SSE2 mismatch at " << i;
        ASSERT_NEAR(dftc_ref.data()[i], dftc_avx2.data()[i], 1e-4) << "AVX2 mismatch at " << i;
    }
}

// --- cast Tests ---

template <typename T>
class CastTest : public ::testing::Test {};

using CastTypes = ::testing::Types<uint8_t, uint16_t, float>;
TYPED_TEST_SUITE(CastTest, CastTypes);

TYPED_TEST(CastTest, MatchesReference) {
    using T = TypeParam;
    constexpr int width = 128;
    constexpr int height = 64;
    constexpr int dstStride = width;
    constexpr int ebpStride = width;
    constexpr float multiplier = 1.0f;
    int peak = 255;
    if constexpr (std::is_same_v<T, uint16_t>) {
        peak = 65535;
    } else if constexpr (std::is_same_v<T, float>) {
        peak = 0; // not used for float
    }

    AlignedBuffer<float> ebp(ebpStride * height);
    AlignedBuffer<T> dst_ref(dstStride * height);
    AlignedBuffer<T> dst_sse2(dstStride * height);
    AlignedBuffer<T> dst_avx2(dstStride * height);

    fill_random(ebp.data(), ebp.size(), 0.0, 255.0);

    Ref::cast(ebp.data(), dst_ref.data(), width, height, dstStride, ebpStride, multiplier, peak);
    SSE2::cast(ebp.data(), dst_sse2.data(), width, height, dstStride, ebpStride, multiplier, peak);
    AVX2::cast(ebp.data(), dst_avx2.data(), width, height, dstStride, ebpStride, multiplier, peak);

    for (size_t i = 0; i < dst_ref.size(); ++i) {
        if constexpr (std::is_floating_point_v<T>) {
             ASSERT_NEAR(dst_ref.data()[i], dst_sse2.data()[i], 1e-5) << "SSE2 mismatch at " << i;
             ASSERT_NEAR(dst_ref.data()[i], dst_avx2.data()[i], 1e-5) << "AVX2 mismatch at " << i;
        } else {
             ASSERT_EQ(dst_ref.data()[i], dst_sse2.data()[i]) << "SSE2 mismatch at " << i;
             ASSERT_EQ(dst_ref.data()[i], dst_avx2.data()[i]) << "AVX2 mismatch at " << i;
        }
    }
}
