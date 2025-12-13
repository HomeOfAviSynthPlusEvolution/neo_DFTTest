#include "core.h"

template<typename T>
static void copyPad(int plane, const unsigned char * src_ptr, int src_stride_bytes, unsigned char * dst_ptr, const DFTTestData * d) noexcept {
    int srcWidth = d->planeWidth[plane];
    int srcHeight = d->planeHeight[plane];
    int dstWidth = d->padWidth[plane];
    int dstHeight = d->padHeight[plane];
    int dstStrideBytes = d->padStride[plane];
    int dstStride = d->padStride[plane] / sizeof(T);

    const int offy = (dstHeight - srcHeight) / 2;
    const int offx = (dstWidth - srcWidth) / 2;

    const unsigned char * scrp0 = src_ptr;
    unsigned char * dstp0 = dst_ptr + dstStrideBytes * offy + offx * sizeof(T);
    for (int h = 0; h < srcHeight; h++)
    {
        memcpy(dstp0, scrp0, srcWidth * sizeof(T));
        scrp0 += src_stride_bytes;
        dstp0 += dstStrideBytes;
    }
    
    T * dstp = reinterpret_cast<T *>(dst_ptr) + dstStride * offy;

    for (int y = offy; y < srcHeight + offy; y++) {
        int w = offx * 2;
        for (int x = 0; x < offx; x++, w--)
            dstp[x] = dstp[w];

        w = offx + srcWidth - 2;
        for (int x = offx + srcWidth; x < dstWidth; x++, w--)
            dstp[x] = dstp[w];

        dstp += dstStride;
    }

    int w = offy * 2;
    for (int y = 0; y < offy; y++, w--)
        memcpy(dst_ptr + dstStrideBytes * y, dst_ptr + dstStrideBytes * w, dstWidth * sizeof(T));

    w = offy + srcHeight - 2;
    for (int y = offy + srcHeight; y < dstHeight; y++, w--)
        memcpy(dst_ptr + dstStrideBytes * y, dst_ptr + dstStrideBytes * w, dstWidth * sizeof(T));
}

static double besselI0(double p) noexcept {
    p /= 2.;
    double n = 1., t = 1., d = 1.;
    int k = 1;
    double v;

    do {
        n *= p;
        d *= k;
        v = n / d;
        t += v * v;
    } while (++k < 15 && v > 1e-8);

    return t;
}

static double getWinValue(const double n, const double size, const int win, const double beta) noexcept {
    switch (win) {
    case 0: // hanning
        return 0.5 - 0.5 * std::cos(2. * M_PI * n / size);
    case 1: // hamming
        return 0.53836 - 0.46164 * std::cos(2. * M_PI * n / size);
    case 2: // blackman
        return 0.42 - 0.5 * std::cos(2. * M_PI * n / size) + 0.08 * std::cos(4. * M_PI * n / size);
    case 3: // 4 term blackman-harris
        return 0.35875 - 0.48829 * std::cos(2. * M_PI * n / size) + 0.14128 * std::cos(4. * M_PI * n / size) - 0.01168 * std::cos(6. * M_PI * n / size);
    case 4: // kaiser-bessel
    {
        const double v = 2. * n / size - 1.;
        return besselI0(M_PI * beta * std::sqrt(1. - v * v)) / besselI0(M_PI * beta);
    }
    case 5: // 7 term blackman-harris
        return 0.27105140069342415 -
               0.433297939234486060 * std::cos(2. * M_PI * n / size) +
               0.218122999543110620 * std::cos(4. * M_PI * n / size) -
               0.065925446388030898 * std::cos(6. * M_PI * n / size) +
               0.010811742098372268 * std::cos(8. * M_PI * n / size) -
               7.7658482522509342E-4 * std::cos(10. * M_PI * n / size) +
               1.3887217350903198E-5 * std::cos(12. * M_PI * n / size);
    case 6: // flat top
        return 0.2810639 - 0.5208972 * std::cos(2. * M_PI * n / size) + 0.1980399 * std::cos(4. * M_PI * n / size);
    case 7: // rectangular
        return 1.;
    case 8: // Bartlett
        return 2. / size * (size / 2. - std::abs(n - size / 2.));
    case 9: // Bartlett-Hann
        return 0.62 - 0.48 * (n / size - 0.5) - 0.38 * std::cos(2. * M_PI * n / size);
    case 10: // Nuttall
        return 0.355768 - 0.487396 * std::cos(2. * M_PI * n / size) + 0.144232 * std::cos(4. * M_PI * n / size) - 0.012604 * std::cos(6. * M_PI * n / size);
    case 11: // Blackman-Nuttall
        return 0.3635819 - 0.4891775 * std::cos(2. * M_PI * n / size) + 0.1365995 * std::cos(4. * M_PI * n / size) - 0.0106411 * std::cos(6. * M_PI * n / size);
    default:
        return 0.;
    }
}

static void normalizeForOverlapAdd(double * VS_RESTRICT hw, const int bsize, const int osize) noexcept {
    double * VS_RESTRICT nw = new double[bsize]();
    const int inc = bsize - osize;

    for (int q = 0; q < bsize; q++) {
        for (int h = q; h >= 0; h -= inc)
            nw[q] += hw[h] * hw[h];
        for (int h = q + inc; h < bsize; h += inc)
            nw[q] += hw[h] * hw[h];
    }

    for (int q = 0; q < bsize; q++)
        hw[q] /= std::sqrt(nw[q]);

    delete[] nw;
}

void createWindow(float * VS_RESTRICT hw, const int tmode, const int smode, const DFTTestData * d) noexcept {
    double * VS_RESTRICT tw = new double[d->tbsize];
    for (int j = 0; j < d->tbsize; j++)
        tw[j] = getWinValue(j + 0.5, d->tbsize, d->twin, d->tbeta);
    if (tmode == 1)
        normalizeForOverlapAdd(tw, d->tbsize, d->tosize);

    double * VS_RESTRICT sw = new double[d->sbsize];
    for (int j = 0; j < d->sbsize; j++)
        sw[j] = getWinValue(j + 0.5, d->sbsize, d->swin, d->sbeta);
    if (smode == 1)
        normalizeForOverlapAdd(sw, d->sbsize, d->sosize);

    const double nscale = 1. / std::sqrt(d->bvolume);
    for (int j = 0; j < d->tbsize; j++)
        for (int k = 0; k < d->sbsize; k++)
            for (int q = 0; q < d->sbsize; q++)
                hw[(j * d->sbsize + k) * d->sbsize + q] = static_cast<float>(tw[j] * sw[k] * sw[q] * nscale);

    delete[] tw;
    delete[] sw;
}

float * parseSigmaLocation(const std::vector<float> s, int & poscnt, const float sigma, const float pfact) {
    float * parray = nullptr;

    if (s.empty()) {
        parray = new float[4];
        parray[0] = 0.0f;
        parray[2] = 1.0f;
        parray[1] = parray[3] = std::pow(sigma, pfact);
        poscnt = 2;
    } else {
        bool found[2] = { false, false };
        poscnt = 0;

        for (int i = 0; i < s.size(); i += 2) {
            const float pos = s[i];

            if (pos < 0.0f || pos > 1.0f)
                throw strdup((std::string{ "sigma location - invalid pos (" } + std::to_string(pos) + ")").c_str());

            if (pos == 0.0f)
                found[0] = true;
            else if (pos == 1.0f)
                found[1] = true;

            poscnt++;
        }

        if (!found[0] || !found[1])
            throw "sigma location - one or more end points not provided";

        parray = new float[poscnt * 2];
        poscnt = 0;

        for (int i = 0; i < s.size(); i += 2) {
            parray[poscnt * 2 + 0] = s[i + 0];
            parray[poscnt * 2 + 1] = std::pow(s[i + 1], pfact);

            poscnt++;
        }

        for (int i = 1; i < poscnt; i++) {
            int j = i;
            const float t0 = parray[j * 2 + 0];
            const float t1 = parray[j * 2 + 1];

            while (j > 0 && parray[(j - 1) * 2] > t0) {
                parray[j * 2 + 0] = parray[(j - 1) * 2 + 0];
                parray[j * 2 + 1] = parray[(j - 1) * 2 + 1];
                j--;
            }

            parray[j * 2 + 0] = t0;
            parray[j * 2 + 1] = t1;
        }
    }

    return parray;
}

float interp(const float pf, const float * pv, const int cnt) noexcept {
    int lidx = 0;
    for (int i = cnt - 1; i >= 0; i--) {
        if (pv[i * 2] <= pf) {
            lidx = i;
            break;
        }
    }

    int hidx = cnt - 1;
    for (int i = 0; i < cnt; i++) {
        if (pv[i * 2] >= pf) {
            hidx = i;
            break;
        }
    }

    const float d0 = pf - pv[lidx * 2];
    const float d1 = pv[hidx * 2] - pf;

    if (hidx == lidx || d0 <= 0.0f)
        return pv[lidx * 2 + 1];
    if (d1 <= 0.0f)
        return pv[hidx * 2 + 1];

    const float tf = d0 / (d0 + d1);
    return pv[lidx * 2 + 1] * (1.0f - tf) + pv[hidx * 2 + 1] * tf;
}

float getSVal(const int pos, const int len, const float * pv, const int cnt, float & pf) noexcept {
    if (len == 1) {
        pf = 0.0f;
        return 1.0f;
    }

    const int ld2 = len / 2;
    if (pos > ld2)
        pf = (len - pos) / static_cast<float>(ld2);
    else
        pf = pos / static_cast<float>(ld2);

    return interp(pf, pv, cnt);
}

void selectFunctions(const unsigned ftype, const unsigned opt, DFTTestData * d) noexcept {
    if (ftype == 0) {
        if (std::abs(d->f0beta - 1.0f) < 0.00005f)
            d->filterCoeffs = filter_c<0>;
        else if (std::abs(d->f0beta - 0.5f) < 0.00005f)
            d->filterCoeffs = filter_c<6>;
        else
            d->filterCoeffs = filter_c<5>;
    } else if (ftype == 1) {
        d->filterCoeffs = filter_c<1>;
    } else if (ftype == 2) {
        d->filterCoeffs = filter_c<2>;
    } else if (ftype == 3) {
        d->filterCoeffs = filter_c<3>;
    } else {
        d->filterCoeffs = filter_c<4>;
    }

    if (d->vi_bytesPerSample == 1) {
        d->copyPad = copyPad<uint8_t>;
        d->func_0 = func_0_c<uint8_t>;
        d->func_1 = func_1_c<uint8_t>;
    } else if (d->vi_bytesPerSample == 2) {
        d->copyPad = copyPad<uint16_t>;
        d->func_0 = func_0_c<uint16_t>;
        d->func_1 = func_1_c<uint16_t>;
    } else {
        d->copyPad = copyPad<float>;
        d->func_0 = func_0_c<float>;
        d->func_1 = func_1_c<float>;
    }

    #ifdef VS_TARGET_CPU_X86
        const int iset = instrset_detect();

        if (opt == 8) {
            d->filterCoeffs = neo_dfttest::GetHighwayFilter(ftype, d->f0beta);
            neo_dfttest::GetHighwayFunc0(d);
            neo_dfttest::GetHighwayFunc1(d);
        } else if ((opt == 0 && iset >= 8) || opt == 3) {
            if (ftype == 0) {
                if (std::abs(d->f0beta - 1.0f) < 0.00005f)
                    d->filterCoeffs = filter_avx2<0>;
                else if (std::abs(d->f0beta - 0.5f) < 0.00005f)
                    d->filterCoeffs = filter_avx2<6>;
                else
                    d->filterCoeffs = filter_avx2<5>;
            } else if (ftype == 1) {
                d->filterCoeffs = filter_avx2<1>;
            } else if (ftype == 2) {
                d->filterCoeffs = filter_avx2<2>;
            } else if (ftype == 3) {
                d->filterCoeffs = filter_avx2<3>;
            } else {
                d->filterCoeffs = filter_avx2<4>;
            }

            if (d->vi_bytesPerSample == 1) {
                d->func_0 = func_0_avx2<uint8_t>;
                d->func_1 = func_1_avx2<uint8_t>;
            } else if (d->vi_bytesPerSample == 2) {
                d->func_0 = func_0_avx2<uint16_t>;
                d->func_1 = func_1_avx2<uint16_t>;
            } else {
                d->func_0 = func_0_avx2<float>;
                d->func_1 = func_1_avx2<float>;
            }
        } else if ((opt == 0 && iset >= 2) || opt == 2) {
            if (ftype == 0) {
                if (std::abs(d->f0beta - 1.0f) < 0.00005f)
                    d->filterCoeffs = filter_sse2<0>;
                else if (std::abs(d->f0beta - 0.5f) < 0.00005f)
                    d->filterCoeffs = filter_sse2<6>;
                else
                    d->filterCoeffs = filter_sse2<5>;
            } else if (ftype == 1) {
                d->filterCoeffs = filter_sse2<1>;
            } else if (ftype == 2) {
                d->filterCoeffs = filter_sse2<2>;
            } else if (ftype == 3) {
                d->filterCoeffs = filter_sse2<3>;
            } else {
                d->filterCoeffs = filter_sse2<4>;
            }

            if (d->vi_bytesPerSample == 1) {
                d->func_0 = func_0_sse2<uint8_t>;
                d->func_1 = func_1_sse2<uint8_t>;
            } else if (d->vi_bytesPerSample == 2) {
                d->func_0 = func_0_sse2<uint16_t>;
                d->func_1 = func_1_sse2<uint16_t>;
            } else {
                d->func_0 = func_0_sse2<float>;
                d->func_1 = func_1_sse2<float>;
            }
        }
    #endif
}
