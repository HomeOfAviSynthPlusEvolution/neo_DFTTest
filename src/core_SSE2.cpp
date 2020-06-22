#ifdef VS_TARGET_CPU_X86
#include "core.h"

#include "vectorclass/vectormath_exp.h"

Vec4f rcpnr_ps(const Vec4f &a) {
  auto r = _mm_rcp_ps(a);
  return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, a), r));
}

template<typename T>
static inline void proc0(const T * s0, const float * s1, float * d, const int p0, const int p1, const float divisor) noexcept;

template<>
inline void proc0(const uint8_t * _s0, const float * _s1, float * d, const int p0, const int p1, const float divisor) noexcept {
    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v += 4) {
            const Vec4f s0 = to_float(Vec4i().load_4uc(_s0 + v));
            const Vec4f s1 = Vec4f().load(_s1 + v);
            (s0 * s1).store(d + v);
        }

        _s0 += p0;
        _s1 += p1;
        d += p1;
    }
}

template<>
inline void proc0(const uint16_t * _s0, const float * _s1, float * d, const int p0, const int p1, const float divisor) noexcept {
    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v += 4) {
            const Vec4f s0 = to_float(Vec4i().load_4us(_s0 + v));
            const Vec4f s1 = Vec4f().load(_s1 + v);
            (s0 * divisor * s1).store(d + v);
        }

        _s0 += p0;
        _s1 += p1;
        d += p1;
    }
}

template<>
inline void proc0(const float * _s0, const float * _s1, float * d, const int p0, const int p1, const float divisor) noexcept {
    for (int u = 0; u < p1; u++) {
        for (int v = 0; v < p1; v += 4) {
            const Vec4f s0 = Vec4f().load(_s0 + v);
            const Vec4f s1 = Vec4f().load(_s1 + v);
            (s0 * 255.0f * s1).store(d + v);
        }

        _s0 += p0;
        _s1 += p1;
        d += p1;
    }
}

static inline void proc1(const float * _s0, const float * _s1, float * _d, const int p0, const int p1) noexcept {
    for (int u = 0; u < p0; u++) {
        for (int v = 0; v < p0; v += 4) {
            const Vec4f s0 = Vec4f().load(_s0 + v);
            const Vec4f s1 = Vec4f().load(_s1 + v);
            const Vec4f d = Vec4f().load(_d + v);
            mul_add(s0, s1, d).store(_d + v);
        }

        _s0 += p0;
        _s1 += p0;
        _d += p1;
    }
}

static inline void proc1Partial(const float * _s0, const float * _s1, float * _d, const int p0, const int p1) noexcept {
    const int regularPart = p0 & -4;

    for (int u = 0; u < p0; u++) {
        int v;

        for (v = 0; v < regularPart; v += 4) {
            const Vec4f s0 = Vec4f().load(_s0 + v);
            const Vec4f s1 = Vec4f().load(_s1 + v);
            const Vec4f d = Vec4f().load(_d + v);
            mul_add(s0, s1, d).store(_d + v);
        }

        const Vec4f s0 = Vec4f().load(_s0 + v);
        const Vec4f s1 = Vec4f().load(_s1 + v);
        const Vec4f d = Vec4f().load(_d + v);
        mul_add(s0, s1, d).store_partial(p0 - v, _d + v);

        _s0 += p0;
        _s1 += p0;
        _d += p1;
    }
}

static inline void removeMean(float * _dftc, const float * _dftgc, const int ccnt, float * _dftc2) noexcept {
    const Vec4f gf = _dftc[0] / _dftgc[0];

    for (int h = 0; h < ccnt; h += 4) {
        const Vec4f dftgc = Vec4f().load_a(_dftgc + h);
        const Vec4f dftc = Vec4f().load_a(_dftc + h);
        const Vec4f dftc2 = gf * dftgc;
        dftc2.store_a(_dftc2 + h);
        (dftc - dftc2).store_a(_dftc + h);
    }
}

static inline void addMean(float * _dftc, const int ccnt, const float * _dftc2) noexcept {
    for (int h = 0; h < ccnt; h += 4) {
        const Vec4f dftc = Vec4f().load_a(_dftc + h);
        const Vec4f dftc2 = Vec4f().load_a(_dftc2 + h);
        (dftc + dftc2).store_a(_dftc + h);
    }
}

template<int type>
void filter_sse2(float * dftc, const float * sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept;

template<>
void filter_sse2<0>(float * dftc, const float * _sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 8) {
        Vec4f dftcLow = Vec4f().load_a(dftc + h);
        Vec4f dftcHigh = Vec4f().load_a(dftc + h + 4);
        Vec4f real = blend4f<0, 2, 4, 6>(dftcLow, dftcHigh);
        Vec4f imag = blend4f<1, 3, 5, 7>(dftcLow, dftcHigh);
        const Vec4f psd = mul_add(real, real, imag * imag);

        const Vec4f sigmasLow = Vec4f().load_a(_sigmas + h);
        const Vec4f sigmasHigh = Vec4f().load_a(_sigmas + h + 4);
        const Vec4f sigmas = blend4f<0, 2, 4, 6>(sigmasLow, sigmasHigh);

        const Vec4f mult = max((psd - sigmas) * rcpnr_ps(psd + 1e-15f), zero_4f());
        real *= mult;
        imag *= mult;

        dftcLow = blend4f<0, 4, 1, 5>(real, imag);
        dftcHigh = blend4f<2, 6, 3, 7>(real, imag);
        dftcLow.store_a(dftc + h);
        dftcHigh.store_a(dftc + h + 4);
    }
}

template<>
void filter_sse2<1>(float * dftc, const float * _sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 8) {
        Vec4f dftcLow = Vec4f().load_a(dftc + h);
        Vec4f dftcHigh = Vec4f().load_a(dftc + h + 4);
        Vec4f real = blend4f<0, 2, 4, 6>(dftcLow, dftcHigh);
        Vec4f imag = blend4f<1, 3, 5, 7>(dftcLow, dftcHigh);
        const Vec4f psd = mul_add(real, real, imag * imag);

        const Vec4f sigmasLow = Vec4f().load_a(_sigmas + h);
        const Vec4f sigmasHigh = Vec4f().load_a(_sigmas + h + 4);
        const Vec4f sigmas = blend4f<0, 2, 4, 6>(sigmasLow, sigmasHigh);

        const Vec4fb flag = (psd < sigmas);
        real = select(flag, zero_4f(), real);
        imag = select(flag, zero_4f(), imag);

        dftcLow = blend4f<0, 4, 1, 5>(real, imag);
        dftcHigh = blend4f<2, 6, 3, 7>(real, imag);
        dftcLow.store_a(dftc + h);
        dftcHigh.store_a(dftc + h + 4);
    }
}

template<>
void filter_sse2<2>(float * _dftc, const float * _sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 4) {
        const Vec4f dftc = Vec4f().load_a(_dftc + h);
        const Vec4f sigmas = Vec4f().load_a(_sigmas + h);
        (dftc * sigmas).store_a(_dftc + h);
    }
}

template<>
void filter_sse2<3>(float * dftc, const float * _sigmas, const int ccnt, const float * _pmin, const float * _pmax, const float * _sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 8) {
        Vec4f dftcLow = Vec4f().load_a(dftc + h);
        Vec4f dftcHigh = Vec4f().load_a(dftc + h + 4);
        Vec4f real = blend4f<0, 2, 4, 6>(dftcLow, dftcHigh);
        Vec4f imag = blend4f<1, 3, 5, 7>(dftcLow, dftcHigh);
        const Vec4f psd = mul_add(real, real, imag * imag);

        const Vec4f sigmasLow = Vec4f().load_a(_sigmas + h);
        const Vec4f sigmasHigh = Vec4f().load_a(_sigmas + h + 4);
        const Vec4f sigmas = blend4f<0, 2, 4, 6>(sigmasLow, sigmasHigh);

        const Vec4f sigmas2Low = Vec4f().load_a(_sigmas2 + h);
        const Vec4f sigmas2High = Vec4f().load_a(_sigmas2 + h + 4);
        const Vec4f sigmas2 = blend4f<0, 2, 4, 6>(sigmas2Low, sigmas2High);

        const Vec4f pminLow = Vec4f().load_a(_pmin + h);
        const Vec4f pminHigh = Vec4f().load_a(_pmin + h + 4);
        const Vec4f pmin = blend4f<0, 2, 4, 6>(pminLow, pminHigh);

        const Vec4f pmaxLow = Vec4f().load_a(_pmax + h);
        const Vec4f pmaxHigh = Vec4f().load_a(_pmax + h + 4);
        const Vec4f pmax = blend4f<0, 2, 4, 6>(pmaxLow, pmaxHigh);

        const Vec4fb flag = (psd >= pmin && psd <= pmax);
        real = select(flag, real * sigmas, real * sigmas2);
        imag = select(flag, imag * sigmas, imag * sigmas2);

        dftcLow = blend4f<0, 4, 1, 5>(real, imag);
        dftcHigh = blend4f<2, 6, 3, 7>(real, imag);
        dftcLow.store_a(dftc + h);
        dftcHigh.store_a(dftc + h + 4);
    }
}

template<>
void filter_sse2<4>(float * dftc, const float * _sigmas, const int ccnt, const float * _pmin, const float * _pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 8) {
        Vec4f dftcLow = Vec4f().load_a(dftc + h);
        Vec4f dftcHigh = Vec4f().load_a(dftc + h + 4);
        Vec4f real = blend4f<0, 2, 4, 6>(dftcLow, dftcHigh);
        Vec4f imag = blend4f<1, 3, 5, 7>(dftcLow, dftcHigh);
        const Vec4f psd = mul_add(real, real, imag * imag) + 1e-15f;

        const Vec4f sigmasLow = Vec4f().load_a(_sigmas + h);
        const Vec4f sigmasHigh = Vec4f().load_a(_sigmas + h + 4);
        const Vec4f sigmas = blend4f<0, 2, 4, 6>(sigmasLow, sigmasHigh);

        const Vec4f pminLow = Vec4f().load_a(_pmin + h);
        const Vec4f pminHigh = Vec4f().load_a(_pmin + h + 4);
        const Vec4f pmin = blend4f<0, 2, 4, 6>(pminLow, pminHigh);

        const Vec4f pmaxLow = Vec4f().load_a(_pmax + h);
        const Vec4f pmaxHigh = Vec4f().load_a(_pmax + h + 4);
        const Vec4f pmax = blend4f<0, 2, 4, 6>(pmaxLow, pmaxHigh);

        const Vec4f mult = sigmas * sqrt(psd * pmax / ((psd + pmin) * (psd + pmax)));
        real *= mult;
        imag *= mult;

        dftcLow = blend4f<0, 4, 1, 5>(real, imag);
        dftcHigh = blend4f<2, 6, 3, 7>(real, imag);
        dftcLow.store_a(dftc + h);
        dftcHigh.store_a(dftc + h + 4);
    }
}

template<>
void filter_sse2<5>(float * dftc, const float * _sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    const Vec4f beta = pmin[0];

    for (int h = 0; h < ccnt; h += 8) {
        Vec4f dftcLow = Vec4f().load_a(dftc + h);
        Vec4f dftcHigh = Vec4f().load_a(dftc + h + 4);
        Vec4f real = blend4f<0, 2, 4, 6>(dftcLow, dftcHigh);
        Vec4f imag = blend4f<1, 3, 5, 7>(dftcLow, dftcHigh);
        const Vec4f psd = mul_add(real, real, imag * imag);

        const Vec4f sigmasLow = Vec4f().load_a(_sigmas + h);
        const Vec4f sigmasHigh = Vec4f().load_a(_sigmas + h + 4);
        const Vec4f sigmas = blend4f<0, 2, 4, 6>(sigmasLow, sigmasHigh);

        const Vec4f mult = pow(max((psd - sigmas) * rcpnr_ps(psd + 1e-15f), zero_4f()), beta);
        real *= mult;
        imag *= mult;

        dftcLow = blend4f<0, 4, 1, 5>(real, imag);
        dftcHigh = blend4f<2, 6, 3, 7>(real, imag);
        dftcLow.store_a(dftc + h);
        dftcHigh.store_a(dftc + h + 4);
    }
}

template<>
void filter_sse2<6>(float * dftc, const float * _sigmas, const int ccnt, const float * pmin, const float * pmax, const float * sigmas2) noexcept {
    for (int h = 0; h < ccnt; h += 8) {
        Vec4f dftcLow = Vec4f().load_a(dftc + h);
        Vec4f dftcHigh = Vec4f().load_a(dftc + h + 4);
        Vec4f real = blend4f<0, 2, 4, 6>(dftcLow, dftcHigh);
        Vec4f imag = blend4f<1, 3, 5, 7>(dftcLow, dftcHigh);
        const Vec4f psd = mul_add(real, real, imag * imag);

        const Vec4f sigmasLow = Vec4f().load_a(_sigmas + h);
        const Vec4f sigmasHigh = Vec4f().load_a(_sigmas + h + 4);
        const Vec4f sigmas = blend4f<0, 2, 4, 6>(sigmasLow, sigmasHigh);

        const Vec4f mult = sqrt(max((psd - sigmas) * rcpnr_ps(psd + 1e-15f), zero_4f()));
        real *= mult;
        imag *= mult;

        dftcLow = blend4f<0, 4, 1, 5>(real, imag);
        dftcHigh = blend4f<2, 6, 3, 7>(real, imag);
        dftcLow.store_a(dftc + h);
        dftcHigh.store_a(dftc + h + 4);
    }
}

template<typename T>
static void cast(const float * ebp, T * dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride, const float multiplier, const int peak) noexcept;

template<>
void cast(const float * ebp, uint8_t * dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride, const float multiplier, const int peak) noexcept {
    for (int y = 0; y < dstHeight; y++) {
        for (int x = 0; x < dstWidth; x += 16) {
            const Vec4i srcp_4i_1 = truncate_to_int(Vec4f().load(ebp + x) + 0.5f);
            const Vec4i srcp_4i_2 = truncate_to_int(Vec4f().load(ebp + x + 4) + 0.5f);
            const Vec4i srcp_4i_3 = truncate_to_int(Vec4f().load(ebp + x + 8) + 0.5f);
            const Vec4i srcp_4i_4 = truncate_to_int(Vec4f().load(ebp + x + 12) + 0.5f);
            const Vec8s srcp_8s_1 = compress_saturated(srcp_4i_1, srcp_4i_2);
            const Vec8s srcp_8s_2 = compress_saturated(srcp_4i_3, srcp_4i_4);
            const Vec16uc srcp = compress_saturated_s2u(srcp_8s_1, srcp_8s_2);
            srcp.stream(dstp + x);
        }

        ebp += ebpStride;
        dstp += dstStride;
    }
}

template<>
void cast(const float * ebp, uint16_t * dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride, const float multiplier, const int peak) noexcept {
    for (int y = 0; y < dstHeight; y++) {
        for (int x = 0; x < dstWidth; x += 8) {
            const Vec4i srcp_4i_1 = truncate_to_int(mul_add(Vec4f().load(ebp + x), multiplier, 0.5f));
            const Vec4i srcp_4i_2 = truncate_to_int(mul_add(Vec4f().load(ebp + x + 4), multiplier, 0.5f));
            const Vec8us srcp = compress_saturated_s2u(srcp_4i_1, srcp_4i_2);
            min(srcp, peak).stream(dstp + x);
        }

        ebp += ebpStride;
        dstp += dstStride;
    }
}

template<>
void cast(const float * ebp, float * dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride, const float multiplier, const int peak) noexcept {
    for (int y = 0; y < dstHeight; y++) {
        for (int x = 0; x < dstWidth; x += 4) {
            const Vec4f srcp = Vec4f().load(ebp + x);
            (srcp * (1.0f / 255.0f)).stream(dstp + x);
        }

        ebp += ebpStride;
        dstp += dstStride;
    }
}

template<typename T>
static inline void dither(const float * ebp, T * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode) noexcept {
    cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak);
}

template<>
void dither<uint8_t>(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode) noexcept {
    dither_c(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, multiplier, peak, dither_mode);
}

template<typename T>
void func_0_sse2(unsigned int thread_id, int plane, const unsigned char * src_ptr, unsigned char * dst_ptr, int dst_stride_bytes, const DFTTestData * d) noexcept {
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

        const T * srcp = reinterpret_cast<const T *>(src_ptr) + srcStride * block_start;
        float * ebpSaved = ebuff + ebpStride * block_start;

        for (int y = block_start; y < block_end; y += d->inc) {
            for (int x = 0; x <= width - d->sbsize; x += d->inc) {
                proc0(srcp + x, d->hw, dftr, srcStride, d->sbsize, d->divisor);

                d->fft->fftwf_execute_dft_r2c(d->ft, dftr, dftc);
                if (d->zmean)
                    removeMean(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(d->dftgc), d->ccnt2, reinterpret_cast<float *>(dftc2));

                d->filterCoeffs(reinterpret_cast<float *>(dftc), d->sigmas, d->ccnt2, d->uf0b ? &d->f0beta : d->pmins, d->pmaxs, d->sigmas2);

                if (d->zmean)
                    addMean(reinterpret_cast<float *>(dftc), d->ccnt2, reinterpret_cast<const float *>(dftc2));
                d->fft->fftwf_execute_dft_c2r(d->fti, dftc, dftr);

                if (d->type & 1) { // spatial overlapping
                    if (!(d->sbsize & 3))
                        proc1(dftr, d->hw, ebpSaved + x, d->sbsize, ebpStride);
                    else
                        proc1Partial(dftr, d->hw, ebpSaved + x, d->sbsize, ebpStride);
                }
                else
                    ebpSaved[x + d->sbd1 * ebpStride + d->sbd1] = dftr[d->sbd1 * d->sbsize + d->sbd1] * d->hw[d->sbd1 * d->sbsize + d->sbd1];
            }

            srcp += srcStride * d->inc;
            ebpSaved += ebpStride * d->inc;
        }
  #ifdef ENABLE_PAR
    });
  #endif

    int dstWidth = d->vi_width;
    int dstHeight = d->vi_height;
    int dstStride = dst_stride_bytes / sizeof(T);
    if (plane > 0) {
        dstWidth >>= d->vi_subSamplingW;
        dstHeight >>= d->vi_subSamplingH;
    }
    T * dstp = reinterpret_cast<T *>(dst_ptr);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (d->dither > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak, d->dither);
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak);
}

template<typename T>
void func_1_sse2(unsigned int thread_id, int plane, const unsigned char * src_ptr, unsigned char * dst_ptr, int dst_stride_bytes, const int pos, const DFTTestData * d) noexcept {
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

        const T * srcp[15] = {};
        for (int i = 0; i < d->tbsize; i++)
            srcp[i] = reinterpret_cast<const T *>(src_ptr + d->padBlockSize[plane] * i) + srcStride * block_start;

        for (int y = block_start; y < block_end; y += d->inc) {
            for (int x = 0; x <= width - d->sbsize; x += d->inc) {
                for (int z = 0; z < d->tbsize; z++)
                    proc0(srcp[z] + x, d->hw + d->barea * z, dftr + d->barea * z, srcStride, d->sbsize, d->divisor);

                d->fft->fftwf_execute_dft_r2c(d->ft, dftr, dftc);
                if (d->zmean)
                    removeMean(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(d->dftgc), d->ccnt2, reinterpret_cast<float *>(dftc2));

                d->filterCoeffs(reinterpret_cast<float *>(dftc), d->sigmas, d->ccnt2, d->uf0b ? &d->f0beta : d->pmins, d->pmaxs, d->sigmas2);

                if (d->zmean)
                    addMean(reinterpret_cast<float *>(dftc), d->ccnt2, reinterpret_cast<const float *>(dftc2));
                d->fft->fftwf_execute_dft_c2r(d->fti, dftc, dftr);

                if (d->type & 1) { // spatial overlapping
                    if (!(d->sbsize & 3))
                        proc1(dftr + pos * d->barea, d->hw + pos * d->barea, ebuff + y * ebpStride + x, d->sbsize, ebpStride);
                    else
                        proc1Partial(dftr + pos * d->barea, d->hw + pos * d->barea, ebuff + y * ebpStride + x, d->sbsize, ebpStride);
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

    int dstWidth = d->vi_width;
    int dstHeight = d->vi_height;
    int dstStride = dst_stride_bytes / sizeof(T);
    if (plane > 0) {
        dstWidth >>= d->vi_subSamplingW;
        dstHeight >>= d->vi_subSamplingH;
    }
    T * dstp = reinterpret_cast<T *>(dst_ptr);
    const float * ebp = ebuff + ebpStride * ((height - dstHeight) / 2) + (width - dstWidth) / 2;
    if (d->dither > 0)
        dither(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak, d->dither);
    else
        cast(ebp, dstp, dstWidth, dstHeight, dstStride, ebpStride, d->multiplier, d->peak);
}

template void func_0_sse2<uint8_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template void func_0_sse2<uint16_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template void func_0_sse2<float>(unsigned int, int, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
template void func_1_sse2<uint8_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
template void func_1_sse2<uint16_t>(unsigned int, int, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
template void func_1_sse2<float>(unsigned int, int, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
#endif
