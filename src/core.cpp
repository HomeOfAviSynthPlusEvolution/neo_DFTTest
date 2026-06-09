#include "core.h"

#include <stdexcept>

template<typename T>
static void copyPad(int plane, DFTPlaneBytes src, DFTMutablePlaneBytes dst, const DFTKernelContext& context) noexcept {
    int srcWidth = context.planes.width[plane];
    int srcHeight = context.planes.height[plane];
    int dstWidth = context.planes.pad_width[plane];
    int dstHeight = context.planes.pad_height[plane];
    int dstStrideBytes = dst.stride_bytes;
    int dstStride = dst.stride_bytes / sizeof(T);

    const int offy = (dstHeight - srcHeight) / 2;
    const int offx = (dstWidth - srcWidth) / 2;

    const unsigned char * scrp0 = src.data;
    unsigned char * dstp0 = dst.data + dstStrideBytes * offy + offx * sizeof(T);
    for (int h = 0; h < srcHeight; h++)
    {
        memcpy(dstp0, scrp0, srcWidth * sizeof(T));
        scrp0 += src.stride_bytes;
        dstp0 += dstStrideBytes;
    }
    
    T * dstp = reinterpret_cast<T *>(dst.data) + dstStride * offy;

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
        memcpy(dst.data + dstStrideBytes * y, dst.data + dstStrideBytes * w, dstWidth * sizeof(T));

    w = offy + srcHeight - 2;
    for (int y = offy + srcHeight; y < dstHeight; y++, w--)
        memcpy(dst.data + dstStrideBytes * y, dst.data + dstStrideBytes * w, dstWidth * sizeof(T));
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

static void normalizeForOverlapAdd(std::vector<double>& hw, const int osize) noexcept {
    const int bsize = static_cast<int>(hw.size());
    std::vector<double> nw(static_cast<std::size_t>(bsize), 0.0);
    const int inc = bsize - osize;

    for (int q = 0; q < bsize; q++) {
        for (int h = q; h >= 0; h -= inc)
            nw[q] += hw[h] * hw[h];
        for (int h = q + inc; h < bsize; h += inc)
            nw[q] += hw[h] * hw[h];
    }

    for (int q = 0; q < bsize; q++)
        hw[q] /= std::sqrt(nw[q]);
}

void createWindow(DFTMutableFloatSpan window, const int tmode, const int smode, const DFTBlockSettings& block, const DFTDerivedGeometry& derived) noexcept {
    std::vector<double> tw(static_cast<std::size_t>(block.temporal_size));
    for (int j = 0; j < block.temporal_size; j++)
        tw[j] = getWinValue(j + 0.5, block.temporal_size, block.temporal_window, block.temporal_beta);
    if (tmode == 1)
        normalizeForOverlapAdd(tw, block.temporal_overlap);

    std::vector<double> sw(static_cast<std::size_t>(block.spatial_size));
    for (int j = 0; j < block.spatial_size; j++)
        sw[j] = getWinValue(j + 0.5, block.spatial_size, block.spatial_window, block.spatial_beta);
    if (smode == 1)
        normalizeForOverlapAdd(sw, block.spatial_overlap);

    const double nscale = 1. / std::sqrt(derived.block_volume);
    for (int j = 0; j < block.temporal_size; j++)
        for (int k = 0; k < block.spatial_size; k++)
            for (int q = 0; q < block.spatial_size; q++)
                window.data[(j * block.spatial_size + k) * block.spatial_size + q] = static_cast<float>(tw[j] * sw[k] * sw[q] * nscale);
}

std::vector<float> parseSigmaLocation(const std::vector<float>& s, const float sigma, const float pfact) {
    if (s.empty()) {
        std::vector<float> parray(4);
        parray[0] = 0.0f;
        parray[2] = 1.0f;
        parray[1] = parray[3] = std::pow(sigma, pfact);
        return parray;
    } else {
        bool found[2] = { false, false };
        int poscnt = 0;

        for (int i = 0; i < s.size(); i += 2) {
            const float pos = s[i];

            if (pos < 0.0f || pos > 1.0f)
                throw std::runtime_error(std::string{ "sigma location - invalid pos (" } + std::to_string(pos) + ")");

            if (pos == 0.0f)
                found[0] = true;
            else if (pos == 1.0f)
                found[1] = true;

            poscnt++;
        }

        if (!found[0] || !found[1])
            throw std::runtime_error("sigma location - one or more end points not provided");

        std::vector<float> parray(static_cast<std::size_t>(poscnt) * 2);
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

        return parray;
    }
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

DFTKernelDispatch selectFunctions(const unsigned ftype, const unsigned opt, const DFTClipFormat& format, const DFTBlockSettings& block) noexcept {
    DFTKernelDispatch kernels {};

    if (ftype == 0) {
        if (std::abs(block.f0_beta - 1.0f) < 0.00005f)
            kernels.filter_coefficients = filter_c<0>;
        else if (std::abs(block.f0_beta - 0.5f) < 0.00005f)
            kernels.filter_coefficients = filter_c<6>;
        else
            kernels.filter_coefficients = filter_c<5>;
    } else if (ftype == 1) {
        kernels.filter_coefficients = filter_c<1>;
    } else if (ftype == 2) {
        kernels.filter_coefficients = filter_c<2>;
    } else if (ftype == 3) {
        kernels.filter_coefficients = filter_c<3>;
    } else {
        kernels.filter_coefficients = filter_c<4>;
    }

    if (format.bytes_per_sample == 1) {
        kernels.copy_pad = copyPad<uint8_t>;
        kernels.process_spatial = func_0_c<uint8_t>;
        kernels.process_temporal = func_1_c<uint8_t>;
    } else if (format.bytes_per_sample == 2) {
        kernels.copy_pad = copyPad<uint16_t>;
        kernels.process_spatial = func_0_c<uint16_t>;
        kernels.process_temporal = func_1_c<uint16_t>;
    } else {
        kernels.copy_pad = copyPad<float>;
        kernels.process_spatial = func_0_c<float>;
        kernels.process_temporal = func_1_c<float>;
    }

    if (opt == 0 || opt == 3 || opt == 8) {
        kernels.filter_coefficients = neo_dfttest::GetHighwayFilter(ftype, block.f0_beta);
        kernels.process_spatial = neo_dfttest::GetHighwayFunc0(format);
        kernels.process_temporal = neo_dfttest::GetHighwayFunc1(format);
    }

    return kernels;
}
