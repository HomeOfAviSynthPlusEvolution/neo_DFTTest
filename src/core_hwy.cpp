// For highway, we need to be able to include this file multiple times.
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/core_hwy.cpp"

#include "hwy/foreach_target.h"
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h" // For Pow, Sqrt, etc.
#include "hwy/aligned_allocator.h"

#include "dft_common.h" // Contains DFTTestData struct and FFTW types
#include "core.h" // For removeMean_c and dither_c declarations
#include <iostream>
#include <type_traits> // Required for std::is_same_v

HWY_BEFORE_NAMESPACE();
namespace neo_dfttest {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

// =============================================================================
// HWY IMPLEMENTATIONS OF CORE FUNCTIONS
// =============================================================================

// Implements Proc0 using Highway
template<typename T>
void Proc0(const T* _s0, const float* _s1, float* d, const int p0, const int p1, const float divisor) {
    using D = hn::ScalableTag<float>;
    const D d_f; // Using d_f for float
    const size_t N = hn::Lanes(d_f);

    const auto divisor_v = hn::Set(d_f, divisor);
    const auto mul_255 = hn::Set(d_f, 255.0f);

    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v += N) {
            auto s1_vec = hn::Load(d_f, _s1 + v);
            decltype(s1_vec) s0_vec;

            if constexpr (std::is_same_v<T, uint8_t>) {
                const hn::Rebind<uint8_t,  decltype(d_f)> d_u8_N; 
                const hn::Rebind<uint16_t, decltype(d_f)> d_u16_N;
                const hn::Rebind<uint32_t, decltype(d_f)> d_u32_N;
                s0_vec = hn::ConvertTo(d_f, hn::PromoteTo(d_u32_N, hn::PromoteTo(d_u16_N, hn::LoadN(d_u8_N, reinterpret_cast<const uint8_t*>(_s0) + v, N))));
            } else if constexpr (std::is_same_v<T, uint16_t>) {
                const hn::Rebind<uint16_t, decltype(d_f)> d_u16_N;
                const hn::Rebind<uint32_t, decltype(d_f)> d_u32_N;
                s0_vec = hn::ConvertTo(d_f, hn::PromoteTo(d_u32_N, hn::LoadN(d_u16_N, reinterpret_cast<const uint16_t*>(_s0) + v, N)));
                s0_vec = hn::Mul(s0_vec, divisor_v);
            } else { // float
                s0_vec = hn::Load(d_f, reinterpret_cast<const float*>(_s0) + v);
                s0_vec = hn::Mul(s0_vec, mul_255);
            }
            
            auto res = hn::Mul(s0_vec, s1_vec);
            hn::StoreU(res, d_f, d + v);
        }

        _s0 += p0;
        _s1 += p1;
        d += p1;
    }
}

// Implements Proc1 using Highway
inline void Proc1(const float* _s0, const float* _s1, float* _d, const int p0, const int p1) {
    using D = hn::ScalableTag<float>;
    const D d_f;
    const size_t N = hn::Lanes(d_f);

    for (int u = 0; u < p0; u++) {
        for (size_t v = 0; v < p0; v += N) {
            const auto s0_vec = hn::Load(d_f, _s0 + v);
            const auto s1_vec = hn::Load(d_f, _s1 + v);
            const auto d_vec = hn::Load(d_f, _d + v);
            hn::StoreU(hn::MulAdd(s0_vec, s1_vec, d_vec), d_f, _d + v);
        }

        _s0 += p0;
        _s1 += p0;
        _d += p1;
    }
}

// Implements Proc1Partial using Highway
inline void Proc1Partial(const float* _s0, const float* _s1, float* _d, const int p0, const int p1) {
    using D = hn::ScalableTag<float>;
    const D d_f;
    const size_t N = hn::Lanes(d_f);

    for (int u = 0; u < p0; u++) {
        size_t v = 0;
        for (; v + N <= static_cast<size_t>(p0); v += N) {
            const auto s0_vec = hn::LoadU(d_f, _s0 + v);
            const auto s1_vec = hn::LoadU(d_f, _s1 + v);
            const auto d_vec = hn::LoadU(d_f, _d + v);
            hn::StoreU(hn::MulAdd(s0_vec, s1_vec, d_vec), d_f, _d + v);
        }

        if (v < static_cast<size_t>(p0)) {
            const size_t remaining = p0 - v;
            const auto s0_vec = hn::LoadU(d_f, _s0 + v);
            const auto s1_vec = hn::LoadU(d_f, _s1 + v);
            const auto d_vec = hn::LoadU(d_f, _d + v);
            auto res = hn::MulAdd(s0_vec, s1_vec, d_vec);
            hn::StoreN(res, d_f, _d + v, remaining); // StoreN handles partial writes safely
        }

        _s0 += p0;
        _s1 += p0;
        _d += p1;
    }
}

// Custom RcpNr using Newton-Raphson iteration
template<class D, class V>
HWY_INLINE V RcpNr(const D d, V v) {
    auto r = hn::ApproximateReciprocal(v);
    // Newton-Raphson iteration: r = r * (2 - v * r)
    auto two = hn::Set(d, 2.0f);
    return hn::Mul(r, hn::NegMulAdd(v, r, two));
}

// Implements the filter logic using Highway
template<int type>
void Filter(float * dftc, const float * _sigmas, const int ccnt, const float * _pmin, const float * _pmax, const float * _sigmas2) {
    using D = hn::ScalableTag<float>;
    const D d_f;
    const size_t N = hn::Lanes(d_f);

    const auto zero = hn::Zero(d_f);
    const auto epsilon = hn::Set(d_f, 1e-15f);

    // Specific to type 5
    decltype(zero) beta = (type == 5) ? hn::Set(d_f, _pmin[0]) : zero;

    for (int h = 0; h < ccnt; h += 2 * N) { // ccnt is count of complex floats, so 2*N floats
        auto r = hn::Undefined(d_f);
        auto i = hn::Undefined(d_f);
        hn::LoadInterleaved2(d_f, dftc + h, r, i);

        auto psd = hn::MulAdd(r, r, hn::Mul(i, i));

        auto sigmas = hn::Undefined(d_f);
        auto dummy = hn::Undefined(d_f);
        hn::LoadInterleaved2(d_f, _sigmas + h, sigmas, dummy);

        if constexpr (type == 0) {
            auto num = hn::Sub(psd, sigmas);
            auto den = hn::Add(psd, epsilon);
            auto rcp = RcpNr(d_f, den);
            auto mult = hn::Max(zero, hn::Mul(num, rcp));
            
            r = hn::Mul(r, mult);
            i = hn::Mul(i, mult);
        }
        else if constexpr (type == 1) {
            auto mask = hn::Lt(psd, sigmas);
            r = hn::IfThenElse(mask, zero, r);
            i = hn::IfThenElse(mask, zero, i);
        }
        else if constexpr (type == 2) {
            r = hn::Mul(r, sigmas);
            i = hn::Mul(i, sigmas);
        }
        else if constexpr (type == 3) {
            auto sigmas2_vec = hn::Undefined(d_f); hn::LoadInterleaved2(d_f, _sigmas2 + h, sigmas2_vec, dummy);
            auto pmin_vec = hn::Undefined(d_f); hn::LoadInterleaved2(d_f, _pmin + h, pmin_vec, dummy);
            auto pmax_vec = hn::Undefined(d_f); hn::LoadInterleaved2(d_f, _pmax + h, pmax_vec, dummy);

            auto mask = hn::And(hn::Ge(psd, pmin_vec), hn::Le(psd, pmax_vec));
            r = hn::Mul(r, hn::IfThenElse(mask, sigmas, sigmas2_vec));
            i = hn::Mul(i, hn::IfThenElse(mask, sigmas, sigmas2_vec));
        }
        else if constexpr (type == 4) {
            auto pmin_vec = hn::Undefined(d_f); hn::LoadInterleaved2(d_f, _pmin + h, pmin_vec, dummy);
            auto pmax_vec = hn::Undefined(d_f); hn::LoadInterleaved2(d_f, _pmax + h, pmax_vec, dummy);

            auto psd_mod = hn::Add(psd, epsilon); // original SSE2 adds epsilon directly to psd
            auto psd_pmax = hn::Mul(psd_mod, pmax_vec);
            auto term1 = hn::Add(psd_mod, pmin_vec);
            auto term2 = hn::Add(psd_mod, pmax_vec);
            auto den = hn::Mul(term1, term2);

            auto rcp = RcpNr(d_f, den);
            auto inner = hn::Mul(psd_pmax, rcp);
            auto mult = hn::Mul(sigmas, hn::Sqrt(inner));

            r = hn::Mul(r, mult);
            i = hn::Mul(i, mult);
        }
        else if constexpr (type == 5) {
            auto num = hn::Sub(psd, sigmas);
            auto den = hn::Add(psd, epsilon);
            auto rcp = RcpNr(d_f, den);
            auto base = hn::Max(zero, hn::Mul(num, rcp));
            
            //  auto mult = hn::Pow(d_f, base, beta); // Use hwy::contrib::math::Pow
            auto mult = hn::Exp(d_f, hn::Mul(hn::Log(d_f, base), beta));

            r = hn::Mul(r, mult);
            i = hn::Mul(i, mult);
        }
        else if constexpr (type == 6) {
            auto num = hn::Sub(psd, sigmas);
            auto den = hn::Add(psd, epsilon);
            auto rcp = RcpNr(d_f, den);
            auto val = hn::Max(zero, hn::Mul(num, rcp));
            auto mult = hn::Sqrt(val);

            r = hn::Mul(r, mult);
            i = hn::Mul(i, mult);
        }

        hn::StoreInterleaved2(r, i, d_f, dftc + h);
    }
}

// Implements Cast for dither == 0 using Highway
template<typename T>
void Cast(const float * ebp, T * dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride, const float multiplier, const int peak) {
    using D = hn::ScalableTag<float>;
    const D d_f;
    const size_t N = hn::Lanes(d_f);
    hn::Rebind<int32_t, decltype(d_f)> d_i32;
    hn::Rebind<int16_t, decltype(d_f)> d_i16;
    hn::Rebind<uint16_t, decltype(d_f)> d_u16;
    hn::Rebind<uint8_t, decltype(d_f)> d_u8;

    for (int y = 0; y < dstHeight; y++) {
        int x = 0;
        for (; x + N <= dstWidth; x += N) { // Iterate aligned blocks
            auto v = hn::LoadU(d_f, ebp + x);
            
            if constexpr (std::is_same_v<T, float>) {
                auto val = hn::Mul(v, hn::Set(d_f, 1.0f/255.0f));
                hn::StoreU(val, d_f, reinterpret_cast<float*>(dstp) + x);
            } else if constexpr (std::is_same_v<T, uint8_t>) {
                auto v_rounded = hn::Add(v, hn::Set(d_f, 0.5f));
                auto v_i32 = hn::ConvertTo(d_i32, v_rounded);

                auto zero = hn::Zero(d_i32);
                auto v255 = hn::Set(d_i32, 255);
                auto v_clamped_i32 = hn::Min(hn::Max(v_i32, zero), v255);

                auto v_i16 = hn::DemoteTo(d_i16, v_clamped_i32);
                auto v_u8 = hn::DemoteTo(d_u8, v_i16);
                hn::StoreN(v_u8, d_u8, reinterpret_cast<uint8_t*>(dstp) + x, N);
            } else if constexpr (std::is_same_v<T, uint16_t>) {
                auto v_scaled_rounded = hn::MulAdd(v, hn::Set(d_f, multiplier), hn::Set(d_f, 0.5f));
                auto v_i32 = hn::ConvertTo(d_i32, v_scaled_rounded);

                auto zero = hn::Zero(d_i32);
                auto v_peak = hn::Set(d_i32, peak);
                auto v_clamped_i32 = hn::Max(v_i32, zero);
                v_clamped_i32 = hn::Min(v_clamped_i32, v_peak);
                
                auto v_u16 = hn::DemoteTo(d_u16, v_clamped_i32);
                hn::StoreN(v_u16, d_u16, dstp + x, N);
            }
        }
        
        // Scalar fallback for remaining pixels
        for (; x < dstWidth; x++) {
             float val = ebp[x];
             if constexpr (std::is_same_v<T, float>) {
                 dstp[x] = val * (1.0f/255.0f);
             } else if constexpr (std::is_same_v<T, uint8_t>) {
                 int i = static_cast<int>(val + 0.5f);
                 dstp[x] = static_cast<uint8_t>(std::max(0, std::min(255, i)));
             } else { // uint16_t
                 int i = static_cast<int>(val * multiplier + 0.5f);
                 dstp[x] = static_cast<uint16_t>(std::max(0, std::min(peak, i)));
             }
        }
        ebp += ebpStride;
        dstp += dstStride;
    }
}

// Implements Dither using Highway
template<typename T>
void Dither(const float * ebp, T * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode, MTRand& rng, float *dither_buff) noexcept {
    if constexpr (std::is_same_v<T, uint8_t>) {
        dither_c(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak, dither_mode, rng, dither_buff);
    } else {
        Cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak);
    }
}

// Highway-optimized RemoveMean
void RemoveMean(float * dftc, const float * dftgc, const int ccnt, float * dftc2) {
    using D = hn::ScalableTag<float>;
    const D d_f;
    const size_t N = hn::Lanes(d_f);

    const auto gf_scalar = dftc[0] / dftgc[0];
    const auto gf = hn::Set(d_f, gf_scalar);

    for (size_t h = 0; h < ccnt; h += N) {
        const auto v_dftgc = hn::LoadU(d_f, dftgc + h);
        const auto v_dftc = hn::LoadU(d_f, dftc + h);
        
        const auto v_dftc2 = hn::Mul(gf, v_dftgc);
        hn::StoreU(v_dftc2, d_f, dftc2 + h);

        const auto v_new_dftc = hn::Sub(v_dftc, v_dftc2);
        hn::StoreU(v_new_dftc, d_f, dftc + h);
    }
}

// Highway-optimized AddMean
void AddMean(float * dftc, const int ccnt, const float * dftc2) {
    using D = hn::ScalableTag<float>;
    const D d_f;
    const size_t N = hn::Lanes(d_f);
    for (size_t h = 0; h < ccnt; h += N) {
        auto v1 = hn::LoadU(d_f, dftc + h);
        auto v2 = hn::LoadU(d_f, dftc2 + h);
        hn::StoreU(hn::Add(v1, v2), d_f, dftc + h);
    }
}

// Implements func_0 functionality using Highway
template<typename T>
void Func0(unsigned int thread_id, int plane, const unsigned char * src_ptr, unsigned char * dst_ptr, int dst_stride_bytes, const DFTTestData * d) {
    float * ebuff = d->ebuff[thread_id];
    const int width = d->padWidth[plane];
    const int height = d->padHeight[plane];
    const int eheight = d->eHeight[plane];
    const int srcStride = d->padStride[plane] / sizeof(T);
    const int ebpStride = d->eStride[plane];
    const int batch_size = d->eBatchSize[plane];

    memset(ebuff, 0, ebpStride * height * sizeof(float));
    
    for (int bk = 0; bk * batch_size < eheight; bk++) {
        auto block_start = bk * batch_size;
        auto block_end = std::min(block_start + batch_size, eheight);

        float * dftr = d->dftr[thread_id] + (((d->bvolume + 7) | 15) + 1) * bk;
        fftwf_complex * dftc = d->dftc[thread_id] + (((d->ccnt + 7) | 15) + 1) * bk;
        fftwf_complex * dftc2 = d->dftc2[thread_id] + (((d->ccnt + 7) | 15) + 1) * bk;

        const T * srcp = reinterpret_cast<const T *>(src_ptr) + srcStride * block_start;
        float * ebpSaved = ebuff + ebpStride * block_start;

        for (int y = block_start; y < block_end; y += d->inc) {
            for (int x = 0; x <= width - d->sbsize; x += d->inc) {
                Proc0(srcp + x, d->hw, dftr, srcStride, d->sbsize, d->divisor);

                d->fft->fftwf_execute_dft_r2c(d->ft, dftr, dftc);
                if (d->zmean)
                    RemoveMean(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(d->dftgc), d->ccnt2, reinterpret_cast<float *>(dftc2));

                d->filterCoeffs(reinterpret_cast<float *>(dftc), d->sigmas, d->ccnt2, d->uf0b ? &d->f0beta : d->pmins, d->pmaxs, d->sigmas2);

                if (d->zmean)
                    AddMean(reinterpret_cast<float *>(dftc), d->ccnt2, reinterpret_cast<const float *>(dftc2));
                d->fft->fftwf_execute_dft_c2r(d->fti, dftc, dftr);

                if (d->type & 1) { // spatial overlapping
                    using D_f = hn::ScalableTag<float>;
                    const size_t N_f = hn::Lanes(D_f()); // Get lane count for float
                    if (!(d->sbsize & (N_f - 1))) // Check alignment relative to Highway's float vector width
                         Proc1(dftr, d->hw, ebpSaved + x, d->sbsize, ebpStride);
                    else
                         Proc1Partial(dftr, d->hw, ebpSaved + x, d->sbsize, ebpStride);
                }
                else
                    ebpSaved[x + d->sbd1 * ebpStride + d->sbd1] = dftr[d->sbd1 * d->sbsize + d->sbd1] * d->hw[d->sbd1 * d->sbsize + d->sbd1];
            }

            srcp += srcStride * d->inc;
            ebpSaved += ebpStride * d->inc;
        }
    }

    int dstWidth = d->planeWidth[plane];
    int dstHeight = d->planeHeight[plane];
    int dstStride = dst_stride_bytes / sizeof(T);
    T * dstp = reinterpret_cast<T *>(dst_ptr);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    
    if (d->dither > 0)
        Dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak, d->dither, *d->rngs[thread_id], d->d_buffs[thread_id]);
    else
        Cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak);
}

// Implements func_1 functionality using Highway (temporal processing)
template<typename T>
void Func1(unsigned int thread_id, int plane, const unsigned char * src_ptr, unsigned char * dst_ptr, int dst_stride_bytes, const int pos, const DFTTestData * d) {
    float * ebuff = d->ebuff[thread_id];
    const int width = d->padWidth[plane];
    const int height = d->padHeight[plane];
    const int eheight = d->eHeight[plane];
    const int srcStride = d->padStride[plane] / sizeof(T);
    const int ebpStride = d->eStride[plane];
    const int batch_size = d->eBatchSize[plane];

    memset(ebuff, 0, ebpStride * height * sizeof(float));
    
  #ifdef ENABLE_PAR
    std::for_each_n(std::execution::par, reinterpret_cast<char*>(0), d->threads, [&](char&idx) {
        int bk = static_cast<int>(reinterpret_cast<intptr_t>(&idx));
  #else
    int bk = 0;
  #endif
        auto block_start = bk * batch_size;
        auto block_end = std::min(block_start + batch_size, eheight);

        float * dftr = d->dftr[thread_id] + (((d->bvolume + 7) | 15) + 1) * bk;
        fftwf_complex * dftc = d->dftc[thread_id] + (((d->ccnt + 7) | 15) + 1) * bk;
        fftwf_complex * dftc2 = d->dftc2[thread_id] + (((d->ccnt + 7) | 15) + 1) * bk;

        const T * srcp[15] = {}; // Max d->tbsize is 15 based on original code comments
        for (int i = 0; i < d->tbsize; i++)
            srcp[i] = reinterpret_cast<const T *>(src_ptr + d->padBlockSize[plane] * i) + srcStride * block_start;

        for (int y = block_start; y < block_end; y += d->inc) {
            for (int x = 0; x <= width - d->sbsize; x += d->inc) {
                for (int z = 0; z < d->tbsize; z++)
                    Proc0(srcp[z] + x, d->hw + d->barea * z, dftr + d->barea * z, srcStride, d->sbsize, d->divisor);

                d->fft->fftwf_execute_dft_r2c(d->ft, dftr, dftc);
                if (d->zmean)
                    RemoveMean(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(d->dftgc), d->ccnt2, reinterpret_cast<float *>(dftc2));

                d->filterCoeffs(reinterpret_cast<float *>(dftc), d->sigmas, d->ccnt2, d->uf0b ? &d->f0beta : d->pmins, d->pmaxs, d->sigmas2);

                if (d->zmean)
                    AddMean(reinterpret_cast<float *>(dftc), d->ccnt2, reinterpret_cast<const float *>(dftc2));
                d->fft->fftwf_execute_dft_c2r(d->fti, dftc, dftr);

                if (d->type & 1) { // spatial overlapping
                    using D_f = hn::ScalableTag<float>;
                    const size_t N_f = hn::Lanes(D_f());
                    if (!(d->sbsize & (N_f - 1)))
                        Proc1(dftr + pos * d->barea, d->hw + pos * d->barea, ebuff + y * ebpStride + x, d->sbsize, ebpStride);
                    else
                        Proc1Partial(dftr + pos * d->barea, d->hw + pos * d->barea, ebuff + y * ebpStride + x, d->sbsize, ebpStride);
                }
                else
                    ebuff[(y + d->sbd1) * ebpStride + x + d->sbd1] = dftr[pos * d->barea + d->sbd1 * d->sbsize + d->sbd1] * d->hw[pos * d->barea + d->sbd1 * d->sbsize + d->sbd1];
            }

            for (int q = 0; q < d->tbsize; q++)
                srcp[q] += srcStride * d->inc;
        }
  #ifdef ENABLE_PAR
    });
  #endif

    int dstWidth = d->planeWidth[plane];
    int dstHeight = d->planeHeight[plane];
    int dstStride = dst_stride_bytes / sizeof(T);
    T * dstp = reinterpret_cast<T *>(dst_ptr);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (d->dither > 0)
        Dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak, d->dither, *d->rngs[thread_id], d->d_buffs[thread_id]);
    else
        Cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak);
}

} // namespace HWY_NAMESPACE
} // namespace neo_dfttest

// Need to call HWY_AFTER_NAMESPACE() here
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "hwy/per_target.h"
namespace neo_dfttest {

HWY_EXPORT_T(Func0_u8, Func0<uint8_t>);
HWY_EXPORT_T(Func0_u16, Func0<uint16_t>);
HWY_EXPORT_T(Func0_f32, Func0<float>);

HWY_EXPORT_T(Func1_u8, Func1<uint8_t>);
HWY_EXPORT_T(Func1_u16, Func1<uint16_t>);
HWY_EXPORT_T(Func1_f32, Func1<float>);

HWY_EXPORT_T(Filter_Type0, Filter<0>);
HWY_EXPORT_T(Filter_Type1, Filter<1>);
HWY_EXPORT_T(Filter_Type2, Filter<2>);
HWY_EXPORT_T(Filter_Type3, Filter<3>);
HWY_EXPORT_T(Filter_Type4, Filter<4>);
HWY_EXPORT_T(Filter_Type5, Filter<5>);
HWY_EXPORT_T(Filter_Type6, Filter<6>);

// New exports for testing
HWY_EXPORT_T(Proc0_u8, Proc0<uint8_t>);
HWY_EXPORT_T(Proc0_u16, Proc0<uint16_t>);
HWY_EXPORT_T(Proc0_f32, Proc0<float>);
HWY_EXPORT(Proc1);
HWY_EXPORT(RemoveMean);
HWY_EXPORT(AddMean);
HWY_EXPORT_T(Cast_u8, Cast<uint8_t>);
HWY_EXPORT_T(Cast_u16, Cast<uint16_t>);
HWY_EXPORT_T(Cast_f32, Cast<float>);

// These functions will be exposed as non-templated entry points
// And will call the HWY_DYNAMIC_POINTER to get the right function.
using FilterFunc = void (*)(float *, const float *, const int, const float *, const float *, const float *);

FilterFunc GetHighwayFilter(int ftype, float f0beta) {
    int64_t chosen_target = hwy::DispatchedTarget();
    const char* target_name = hwy::TargetName(chosen_target);
    std::cout << "Instruction used: " << target_name << std::endl;
    if (ftype == 0) {
        if (std::abs(f0beta - 1.0f) < 0.00005f) return HWY_DYNAMIC_POINTER(Filter_Type0);
        else if (std::abs(f0beta - 0.5f) < 0.00005f) return HWY_DYNAMIC_POINTER(Filter_Type6);
        else return HWY_DYNAMIC_POINTER(Filter_Type5);
    } else if (ftype == 1) return HWY_DYNAMIC_POINTER(Filter_Type1);
    else if (ftype == 2) return HWY_DYNAMIC_POINTER(Filter_Type2);
    else if (ftype == 3) return HWY_DYNAMIC_POINTER(Filter_Type3);
    else return HWY_DYNAMIC_POINTER(Filter_Type4);
}

void GetHighwayFunc0(DFTTestData* d) {
    if (d->vi_bytesPerSample == 1) d->func_0 = HWY_DYNAMIC_POINTER(Func0_u8);
    else if (d->vi_bytesPerSample == 2) d->func_0 = HWY_DYNAMIC_POINTER(Func0_u16);
    else d->func_0 = HWY_DYNAMIC_POINTER(Func0_f32);
}

void GetHighwayFunc1(DFTTestData* d) {
    if (d->vi_bytesPerSample == 1) d->func_1 = HWY_DYNAMIC_POINTER(Func1_u8);
    else if (d->vi_bytesPerSample == 2) d->func_1 = HWY_DYNAMIC_POINTER(Func1_u16);
    else d->func_1 = HWY_DYNAMIC_POINTER(Func1_f32);
}

// Getters for internal testing
using Proc0Func_u8 = void (*)(const uint8_t*, const float*, float*, int, int, float);
Proc0Func_u8 GetHighwayProc0_u8() { return HWY_DYNAMIC_POINTER(Proc0_u8); }

using Proc0Func_u16 = void (*)(const uint16_t*, const float*, float*, int, int, float);
Proc0Func_u16 GetHighwayProc0_u16() { return HWY_DYNAMIC_POINTER(Proc0_u16); }

using Proc0Func_f32 = void (*)(const float*, const float*, float*, int, int, float);
Proc0Func_f32 GetHighwayProc0_f32() { return HWY_DYNAMIC_POINTER(Proc0_f32); }

using Proc1Func = void (*)(const float*, const float*, float*, int, int);
Proc1Func GetHighwayProc1() { return HWY_DYNAMIC_POINTER(Proc1); }

using RemoveMeanFunc = void (*)(float*, const float*, int, float*);
RemoveMeanFunc GetHighwayRemoveMean() { return HWY_DYNAMIC_POINTER(RemoveMean); }

using AddMeanFunc = void (*)(float*, int, const float*);
AddMeanFunc GetHighwayAddMean() { return HWY_DYNAMIC_POINTER(AddMean); }

using CastFunc_u8 = void (*)(const float *, uint8_t *, int, int, int, int, float, int);
CastFunc_u8 GetHighwayCast_u8() { return HWY_DYNAMIC_POINTER(Cast_u8); }

using CastFunc_u16 = void (*)(const float *, uint16_t *, int, int, int, int, float, int);
CastFunc_u16 GetHighwayCast_u16() { return HWY_DYNAMIC_POINTER(Cast_u16); }

using CastFunc_f32 = void (*)(const float *, float *, int, int, int, int, float, int);
CastFunc_f32 GetHighwayCast_f32() { return HWY_DYNAMIC_POINTER(Cast_f32); }

FilterFunc GetHighwayFilter_Type0() { return HWY_DYNAMIC_POINTER(Filter_Type0); }
FilterFunc GetHighwayFilter_Type1() { return HWY_DYNAMIC_POINTER(Filter_Type1); }
FilterFunc GetHighwayFilter_Type2() { return HWY_DYNAMIC_POINTER(Filter_Type2); }
FilterFunc GetHighwayFilter_Type3() { return HWY_DYNAMIC_POINTER(Filter_Type3); }
FilterFunc GetHighwayFilter_Type4() { return HWY_DYNAMIC_POINTER(Filter_Type4); }
FilterFunc GetHighwayFilter_Type5() { return HWY_DYNAMIC_POINTER(Filter_Type5); }
FilterFunc GetHighwayFilter_Type6() { return HWY_DYNAMIC_POINTER(Filter_Type6); }

} // namespace neo_dfttest

#endif // HWY_ONCE
