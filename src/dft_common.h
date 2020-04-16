#include <cstdint>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <execution>
#include <VSHelper.h>
#include "fftwlite.h"

#ifndef _WIN32
  #define _aligned_malloc(a,b) aligned_alloc(b,a)
  #define _aligned_free(a) free(a)
#endif

#define EXTRA(a,b) (((a) % (b)) ? ((b) - ((a) % (b))) : 0)

struct DFTTestData {
    FFTFunctionPointers* fft;
    int vi_numPlanes;
    int vi_bytesPerSample, vi_bitsPerSample;
    bool vi_integer;
    int vi_width, vi_height;
    int vi_subSamplingW, vi_subSamplingH;
    int sbsize {16}, sosize{12}, tbsize {3}, tosize {0}, swin {0}, twin {7};
    float sbeta {2.5f}, tbeta {2.5f}, f0beta {1.0f};
    bool zmean {true};
    int dither {0};
    int threads {6};
    int process[3];
    float divisor, multiplier;
    int peak, barea, bvolume, ccnt, type, sbd1, ccnt2, inc;
    bool uf0b;
    std::vector<unsigned char *> pad[3] {std::vector<unsigned char *>(), std::vector<unsigned char *>(), std::vector<unsigned char *>()};
    int padWidth[3], padHeight[3], padStride[3], padBlockSize[3], eStride[3], eHeight[3], eBatchSize[3];
    float * hw {nullptr}, * sigmas {nullptr}, * sigmas2 {nullptr}, * pmins {nullptr}, * pmaxs {nullptr};
    fftwf_complex * dftgc {nullptr};
    fftwf_plan ft {nullptr}, fti {nullptr};

    std::vector<float *> ebuff;
    std::vector<float *> dftr;
    std::vector<fftwf_complex *> dftc, dftc2;

    void (*copyPad)(int plane, const unsigned char *, int, unsigned char *, const DFTTestData *) noexcept;
    void (*filterCoeffs)(float *, const float *, const int, const float *, const float *, const float *) noexcept;
    void (*func_0)(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const DFTTestData *) noexcept;
    void (*func_1)(unsigned int thread_id, int plane, const unsigned char *, unsigned char *, int, const int, const DFTTestData *) noexcept;
};

struct NPInfo {
    int fn, b, y, x;
};

void selectFunctions(const unsigned ftype, const unsigned opt, DFTTestData * d) noexcept;
void createWindow(float * VS_RESTRICT hw, const int tmode, const int smode, const DFTTestData * d) noexcept;
float * parseSigmaLocation(const std::vector<float> s, int & poscnt, const float sigma, const float pfact);
float interp(const float pf, const float * pv, const int cnt) noexcept;
float getSVal(const int pos, const int len, const float * pv, const int cnt, float & pf) noexcept;
void removeMean_c(float * VS_RESTRICT dftc, const float * dftgc, const int ccnt, float * VS_RESTRICT dftc2) noexcept;
template<typename T>
void proc0_c(const T * s0, const float * s1, float * VS_RESTRICT d, const int p0, const int p1, const float divisor) noexcept;
void dither_c(const float * ebp, uint8_t * VS_RESTRICT dstp, const int dstWidth, const int dstHeight, const int dstStride, const int ebpStride,
                 const float multiplier, const int peak, const int dither_mode) noexcept;
