/*
 * Copyright 2020 Xinyue Lu
 *
 * DFTTest bridge.
 *
 */

#pragma once

#include <locale>
#include <numeric>
#include "dft_common.h"
#include "version.hpp"

struct DFTTest final : Filter {
  InDelegator* _in;
  DSVideoInfo out_vi;
  DFTTestData ep;
  FFTFunctionPointers fft;
  int fft_threads {2};
  IScriptEnvironment* avs_env{ nullptr };
  bool has_at_least_v12{ false };

  std::mutex thread_check_mutex;
  std::vector<int> thread_id_store;

  const char* VSName() const override { return "DFTTest"; }
  const char* AVSName() const override { return "neo_dfttest"; }
  const MtMode AVSMode() const override { return MT_NICE_FILTER; }
  const VSFilterMode VSMode() const override { return fmParallel; }
  const std::vector<Param> Params() const override {
    return std::vector<Param> {
      Param {"clip", Clip, false, true, true, false},
      Param {"ftype", Integer},
      Param {"sigma", Float},
      Param {"sigma2", Float},
      Param {"pmin", Float},
      Param {"pmax", Float},
      Param {"sbsize", Integer},
      Param {"smode", Integer},
      Param {"sosize", Integer},
      Param {"tbsize", Integer},
      Param {"tmode", Integer},
      Param {"tosize", Integer},
      Param {"swin", Integer},
      Param {"twin", Integer},
      Param {"sbeta", Float},
      Param {"tbeta", Float},
      Param {"zmean", Boolean},
      Param {"f0beta", Float},
      Param {"nlocation", Integer, true},
      Param {"alpha", Float},
      Param {"slocation", Float, true},
      Param {"ssx", Float, true},
      Param {"ssy", Float, true},
      Param {"sst", Float, true},
      Param {"ssystem", Integer},
      // Param {"sfile", String}, // From AVS version
      // Param {"sfile2", String}, // From AVS version
      Param {"dither", Integer}, // From AVS version
      Param {"planes", Integer, true, false, true},
      Param {"y", Integer, false, true, false},
      Param {"u", Integer, false, true, false},
      Param {"v", Integer, false, true, false},

      Param {"opt", Integer},
      Param {"threads", Integer},
      Param {"fft_threads", Integer}
    };
  }
  void Initialize(InDelegator* in, DSVideoInfo in_vi, FetchFrameFunctor* fetch_frame) override
  {
    // in_vi and fetch_frame are useless for source filter
    Filter::Initialize(in, in_vi, fetch_frame);

    this->avs_env = static_cast<IScriptEnvironment*>(in->GetEnv());
    this->has_at_least_v12 = in->IsAVS12();
    fft.load();

    ep.fft = &fft;
    ep.vi_bitsPerSample = in_vi.Format.BitsPerSample;
    ep.vi_bytesPerSample = in_vi.Format.BytesPerSample;
    ep.vi_integer = in_vi.Format.IsInteger;
    ep.vi_numPlanes = in_vi.Format.Planes;
    ep.vi_width = in_vi.Width;
    ep.vi_height = in_vi.Height;
    ep.vi_subSamplingH = in_vi.Format.SSH;
    ep.vi_subSamplingW = in_vi.Format.SSW;

    int ftype = 0;
    float sigma = 8.0f, sigma2 = 8.0f;
    float pmin = 0.0f, pmax = 500.0f;
    int smode = 1, tmode = 0;
    int opt = 0;

    in->Read("ftype", ftype);
    in->Read("sigma", sigma);
    in->Read("sigma2", sigma2);
    in->Read("pmin", pmin);
    in->Read("pmax", pmax);
    in->Read("smode", smode);
    in->Read("tmode", tmode);
    in->Read("sbsize", ep.sbsize);
    in->Read("sosize", ep.sosize);
    in->Read("tbsize", ep.tbsize);
    in->Read("tosize", ep.tosize);
    in->Read("swin", ep.swin);
    in->Read("twin", ep.twin);
    in->Read("sbeta", ep.sbeta);
    in->Read("tbeta", ep.tbeta);
    in->Read("zmean", ep.zmean);
    in->Read("f0beta", ep.f0beta);

    in->Read("opt", opt);
    in->Read("threads", ep.threads);

    in->Read("fft_threads", fft_threads);
    if (fft_threads < 1)
      fft_threads = 1;

    if (fft_threads > 1 && fft.has_threading()) {
      fft.fftwf_init_threads();
      fft.fftwf_plan_with_nthreads(fft_threads);
    }

    std::vector<int> nlocation;
    float alpha = ftype == 0 ? 5.0f : 7.0f;
    std::vector<float> slocation, ssx, ssy, sst;
    int ssystem = 0;
    std::string tmpstr;
    try { in->Read("nlocation", nlocation); }
    catch (const char *) {
      in->Read("nlocation", tmpstr);
      parse_array_string(tmpstr, nlocation);
    }

    in->Read("alpha", alpha);

    try { in->Read("slocation", slocation); }
    catch (const char *) {
      tmpstr.clear();
      in->Read("slocation", tmpstr);
      parse_array_string(tmpstr, slocation);
    }

    try { in->Read("ssx", ssx); }
    catch (const char *) {
      tmpstr.clear();
      in->Read("ssx", tmpstr);
      parse_array_string(tmpstr, ssx);
    }

    try { in->Read("ssy", ssy); }
    catch (const char *) {
      tmpstr.clear();
      in->Read("ssy", tmpstr);
      parse_array_string(tmpstr, ssy);
    }

    try { in->Read("sst", sst); }
    catch (const char *) {
      tmpstr.clear();
      in->Read("sst", tmpstr);
      parse_array_string(tmpstr, sst);
    }

    in->Read("ssystem", ssystem);
    in->Read("dither", ep.dither);

    try {
      ep.process[0] =
      ep.process[1] =
      ep.process[2] =
      ep.process[3] = 2;
      std::vector<int> user_planes {0, 1, 2};
      in->Read("planes", user_planes);
      for (auto &&p : user_planes)
      {
        if (p < in_vi.Format.Planes)
          ep.process[p] = 3;
        else
          throw "plane index out of range";
      }
    }
    catch (const char *) {
      ep.process[0] =
      ep.process[1] =
      ep.process[2] = 3;
      in->Read("y", ep.process[0]);
      in->Read("u", ep.process[1]);
      in->Read("v", ep.process[2]);
    }


    if (in_vi.Width <= 0 || in_vi.Height <= 0)
      throw "only constant format input supported";

    if ((in_vi.Format.IsInteger && in_vi.Format.BitsPerSample > 16) ||
        (in_vi.Format.IsFloat   && in_vi.Format.BitsPerSample != 32))
      throw "only 8-16 bit integer and 32 bit float input supported";

    if (ftype < 0 || ftype > 4)
      throw "ftype must be 0, 1, 2, 3, or 4";

    if (ep.sbsize < 1)
      throw "sbsize must be greater than or equal to 1";

    if (smode < 0 || smode > 1)
      throw "smode must be 0 or 1";

    if (smode == 0 && !(ep.sbsize & 1))
      throw "sbsize must be odd when using smode=0";

    if (smode == 0)
      ep.sosize = 0;

    if (ep.sosize < 0 || ep.sosize >= ep.sbsize)
      throw "sosize must be between 0 and sbsize-1 (inclusive)";

    if (ep.sosize > ep.sbsize / 2 && ep.sbsize % (ep.sbsize - ep.sosize) != 0)
      throw "spatial overlap greater than 50% requires that sbsize-sosize is a divisor of sbsize";

    if (ep.tbsize < 1 || ep.tbsize > 15)
      throw "tbsize must be between 1 and 15 (inclusive)";

    if (tmode != 0)
      throw "tmode must be 0. tmode=1 is not implemented";

    if (tmode == 0 && !(ep.tbsize & 1))
      throw "tbsize must be odd when using tmode=0";

    if (tmode == 0)
      ep.tosize = 0;

    if (ep.tosize < 0 || ep.tosize >= ep.tbsize)
      throw "tosize must be between 0 and tbsize-1 (inclusive)";

    if (ep.tosize > ep.tbsize / 2 && ep.tbsize % (ep.tbsize - ep.tosize) != 0)
      throw "temporal overlap greater than 50% requires that tbsize-tosize is a divisor of tbsize";

    if (ep.tbsize > in_vi.Frames)
      throw "tbsize must be less than or equal to the number of frames in the clip";

    if (ep.swin < 0 || ep.swin > 11)
      throw "swin must be between 0 and 11 (inclusive)";

    if (ep.twin < 0 || ep.twin > 11)
      throw "twin must be between 0 and 11 (inclusive)";

    if (nlocation.size() & 3)
      throw "the number of elements in nlocation must be a multiple of 4";

    if (alpha <= 0.0f)
      throw "alpha must be greater than 0.0";

    if (slocation.size() & 1)
      throw "the number of elements in slocation must be even";

    if (ssx.size() & 1)
      throw "the number of elements in ssx must be even";

    if (ssy.size() & 1)
      throw "the number of elements in ssy must be even";

    if (sst.size() & 1)
      throw "the number of elements in sst must be even";

    if (ssystem < 0 || ssystem > 1)
      throw "ssystem must be 0 or 1";

    if (opt < 0 || opt > 3)
      throw "opt must be 0, 1, 2, or 3";

    if (ep.threads <= 0)
      ep.threads = 4;

    if (ep.threads > 16)
      ep.threads = 16;

    #ifndef ENABLE_PAR
      ep.threads = 1;
    #endif

    selectFunctions(ftype, opt, &ep);

    if (in_vi.Format.IsInteger) {
      ep.multiplier = static_cast<float>(1 << (in_vi.Format.BitsPerSample - 8));
      ep.divisor = 1.0f / ep.multiplier;
      ep.peak = (1 << in_vi.Format.BitsPerSample) - 1;
    }

    if (ftype != 0)
      ep.f0beta = 1.0f;

    ep.barea = ep.sbsize * ep.sbsize;
    ep.bvolume = ep.barea * ep.tbsize;
    ep.ccnt = (ep.sbsize / 2 + 1) * ep.sbsize * ep.tbsize;
    ep.ccnt2 = ep.ccnt * 2;
    ep.type = tmode * 4 + (ep.tbsize > 1 ? 2 : 0) + smode;
    ep.sbd1 = ep.sbsize / 2;
    ep.uf0b = (std::abs(ep.f0beta - 1.0f) < 0.00005f) ? false : true;
    ep.inc = (ep.type & 1) ? ep.sbsize - ep.sosize : 1;

    for (int plane = 0; plane < ep.vi_numPlanes; plane++) {
      const int width = ep.vi_width >> (plane ? ep.vi_subSamplingW : 0);
      const int height = ep.vi_height >> (plane ? ep.vi_subSamplingH : 0);

      if (smode == 0) {
        const int ae = (ep.sbsize >> 1) << 1;
        ep.padWidth[plane] = width + ae;
        ep.padHeight[plane] = height + ae;
        ep.eHeight[plane] = height;
      } else {
        const int ae = std::max(ep.sbsize - ep.sosize, ep.sosize) * 2;
        ep.padWidth[plane] = width + EXTRA(width, ep.sbsize) + ae;
        ep.padHeight[plane] = height + EXTRA(height, ep.sbsize) + ae;
        ep.eHeight[plane] = (ep.padHeight[plane] - ep.sosize) / (ep.sbsize - ep.sosize) * (ep.sbsize - ep.sosize);
      }
      // round up to 64 for AVX512 alignment
      ep.padStride[plane] = ((ep.padWidth[plane] * in_vi.Format.BytesPerSample - 1) | (FRAME_ALIGN - 1)) + 1;
      ep.padBlockSize[plane] = ep.padStride[plane] * ep.padHeight[plane];
      ep.eStride[plane] = ((ep.padWidth[plane] * sizeof(float) - 1) | (FRAME_ALIGN - 1)) + 1;
      ep.eBatchSize[plane] = ((ep.eHeight[plane] - 1) / ep.threads / ep.inc + 1) * ep.inc;
    }

    ep.hw = (float*)_aligned_malloc((ep.bvolume + 7) * sizeof(float), FRAME_ALIGN);
    if (!ep.hw)
      throw "malloc failure (hw)";
    createWindow(ep.hw, tmode, smode, &ep);

    struct AlignedDeleter {
      void operator ()(void* p) const { _aligned_free(p); }
    };

    float * dftgr = (float*)_aligned_malloc((ep.bvolume + 7) * sizeof(float), FRAME_ALIGN);
    std::unique_ptr<float, AlignedDeleter> dftgr_smart(dftgr);
    ep.dftgc = (fftwf_complex*)_aligned_malloc((ep.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN);
    if (!dftgr || !ep.dftgc)
      throw "malloc failure (dftgr/dftgc)";

    {
      GlobalLockGuard fftw_lock(this->avs_env, "fftw", this->has_at_least_v12);
      if (ep.tbsize > 1) {
        ep.ft = fft.fftwf_plan_dft_r2c_3d(ep.tbsize, ep.sbsize, ep.sbsize, dftgr, ep.dftgc, FFTW_PATIENT | FFTW_DESTROY_INPUT);
        ep.fti = fft.fftwf_plan_dft_c2r_3d(ep.tbsize, ep.sbsize, ep.sbsize, ep.dftgc, dftgr, FFTW_PATIENT | FFTW_DESTROY_INPUT);
      } else {
        ep.ft = fft.fftwf_plan_dft_r2c_2d(ep.sbsize, ep.sbsize, dftgr, ep.dftgc, FFTW_PATIENT | FFTW_DESTROY_INPUT);
        ep.fti = fft.fftwf_plan_dft_c2r_2d(ep.sbsize, ep.sbsize, ep.dftgc, dftgr, FFTW_PATIENT | FFTW_DESTROY_INPUT);
      }
    }

    float wscale = 0.0f;

    const float * hwT = ep.hw;
    float * VS_RESTRICT dftgrT = dftgr;
    for (int s = 0; s < ep.tbsize; s++) {
      for (int i = 0; i < ep.sbsize; i++) {
        for (int k = 0; k < ep.sbsize; k++) {
          dftgrT[k] = 255.0f * hwT[k];
          wscale += hwT[k] * hwT[k];
        }
        hwT += ep.sbsize;
        dftgrT += ep.sbsize;
      }
    }
    fft.fftwf_execute_dft_r2c(ep.ft, dftgr, ep.dftgc);
    dftgr_smart.reset();

    wscale = 1.0f / wscale;
    const float wscalef = (ftype < 2) ? wscale : 1.0f;

    ep.sigmas = (float*)_aligned_malloc((ep.ccnt2 + 7) * sizeof(float), FRAME_ALIGN);
    ep.sigmas2 = (float*)_aligned_malloc((ep.ccnt2 + 7) * sizeof(float), FRAME_ALIGN);
    ep.pmins = (float*)_aligned_malloc((ep.ccnt2 + 7) * sizeof(float), FRAME_ALIGN);
    ep.pmaxs = (float*)_aligned_malloc((ep.ccnt2 + 7) * sizeof(float), FRAME_ALIGN);
    if (!ep.sigmas || !ep.sigmas2 || !ep.pmins || !ep.pmaxs)
      throw "malloc failure (sigmas/sigmas2/pmins/pmaxs)";

    if (!slocation.empty() || !ssx.empty() || !ssy.empty() || !sst.empty()) {
      int ndim = 3;
      if (ep.tbsize == 1)
        ndim -= 1;
      if (ep.sbsize == 1)
        ndim -= 2;

      const float ndiv = 1.0f / ndim;
      int tcnt = 0, sycnt = 0, sxcnt = 0;
      float * tdata, * sydata, * sxdata;

      if (!slocation.empty()) {
        tdata = parseSigmaLocation(slocation, tcnt, sigma, ssystem ? 1.0f : ndiv);
        sydata = parseSigmaLocation(slocation, sycnt, sigma, ssystem ? 1.0f : ndiv);
        sxdata = parseSigmaLocation(slocation, sxcnt, sigma, ssystem ? 1.0f : ndiv);
      } else {
        tdata = parseSigmaLocation(sst, tcnt, sigma, ndiv);
        sydata = parseSigmaLocation(ssy, sycnt, sigma, ndiv);
        sxdata = parseSigmaLocation(ssx, sxcnt, sigma, ndiv);
      }
      auto t_smart = std::unique_ptr<float[]>(tdata);
      auto sx_smart = std::unique_ptr<float[]>(sxdata);
      auto sy_smart = std::unique_ptr<float[]>(sydata);

      const int cpx = ep.sbsize / 2 + 1;
      float pft, pfy, pfx;

      for (int z = 0; z < ep.tbsize; z++) {
        const float tval = getSVal(z, ep.tbsize, tdata, tcnt, pft);

        for (int y = 0; y < ep.sbsize; y++) {
          const float syval = getSVal(y, ep.sbsize, sydata, sycnt, pfy);

          for (int x = 0; x < cpx; x++) {
            const float sxval = getSVal(x, ep.sbsize, sxdata, sxcnt, pfx);
            float val;

            if (ssystem) {
              const float dw = std::sqrt((pft * pft + pfy * pfy + pfx * pfx) / ndim);
              val = interp(dw, tdata, tcnt);
            } else {
              val = tval * syval * sxval;
            }

            const int pos = ((z * ep.sbsize + y) * cpx + x) * 2;
            ep.sigmas[pos] = ep.sigmas[pos + 1] = val / wscalef;
          }
        }
      }

      t_smart.reset();
      sy_smart.reset();
      sx_smart.reset();
    } else {
      for (int i = 0; i < ep.ccnt2; i++)
        ep.sigmas[i] = sigma / wscalef;
    }

    for (int i = 0; i < ep.ccnt2; i++) {
      ep.sigmas2[i] = sigma2 / wscalef;
      ep.pmins[i] = pmin / wscale;
      ep.pmaxs[i] = pmax / wscale;
    }

    if (!nlocation.empty() && ftype < 2) {
      memset(ep.sigmas, 0, ep.ccnt2 * sizeof(float));

      float * VS_RESTRICT hw2 = (float*)_aligned_malloc((ep.bvolume + 7) * sizeof(float), FRAME_ALIGN);
      if (!hw2)
        throw "malloc failure (hw2)";

      std::unique_ptr<float, AlignedDeleter> hw2_smart(hw2);
      createWindow(hw2, 0, 0, &ep);

      float * VS_RESTRICT dftr = (float*)_aligned_malloc((ep.bvolume + 7) * sizeof(float), FRAME_ALIGN);
      fftwf_complex * dftgc2 = (fftwf_complex*)_aligned_malloc((ep.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN);
      if (!dftr || !dftgc2)
        throw "malloc failure (dftr/dftgc2)";
      std::unique_ptr<float, AlignedDeleter> dftr_smart(dftr);
      std::unique_ptr<fftwf_complex, AlignedDeleter> dftgc2_smart(dftgc2);

      float wscale2 = 0.0f;
      int w = 0;
      for (int s = 0; s < ep.tbsize; s++) {
        for (int i = 0; i < ep.sbsize; i++) {
          for (int k = 0; k < ep.sbsize; k++, w++) {
            dftr[w] = 255.0f * hw2[w];
            wscale2 += hw2[w] * hw2[w];
          }
        }
      }
      wscale2 = 1.0f / wscale2;
      fft.fftwf_execute_dft_r2c(ep.ft, dftr, dftgc2);

      int nnpoints = 0;
      NPInfo * npts = new NPInfo[500];
      auto npts_smart = std::unique_ptr<NPInfo[]>(npts);

      for (int i = 0; i < nlocation.size(); i += 4) {
        const int fn = nlocation[i + 0];
        const int b = nlocation[i + 1];
        const int y = nlocation[i + 2];
        const int x = nlocation[i + 3];

        if (fn < 0 || fn > in_vi.Frames - ep.tbsize)
          throw strdup((std::string{ "invalid frame number in nlocation (" } + std::to_string(fn) + ")").c_str());

        if (b < 0 || b >= ep.vi_numPlanes)
          throw strdup((std::string{ "invalid plane number in nlocation (" } + std::to_string(b) + ")").c_str());

        const int height = ep.vi_height >> (b > 0 ? ep.vi_subSamplingH : 0);
        if (y < 0 || y > height - ep.sbsize)
          throw strdup((std::string{ "invalid y pos in nlocation (" } + std::to_string(y) + ")").c_str());

        const int width = ep.vi_width >> (b > 0 ? ep.vi_subSamplingW : 0);
        if (x < 0 || x > width - ep.sbsize)
          throw strdup((std::string{ "invalid x pos in nlocation (" } + std::to_string(x) + ")").c_str());

        if (nnpoints >= 500)
          throw "maximum number of entries in nlocation is 500";

        npts[nnpoints].fn = fn;
        npts[nnpoints].b = b;
        npts[nnpoints].y = y;
        npts[nnpoints].x = x;
        nnpoints++;
      }

      fftwf_complex * dftc = (fftwf_complex*)_aligned_malloc((ep.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN);
      fftwf_complex * dftc2 = (fftwf_complex*)_aligned_malloc((ep.ccnt + 7) * sizeof(fftwf_complex), FRAME_ALIGN);
      if (!dftc || !dftc2)
        throw "malloc failure (dftc/dftc2)";
      std::unique_ptr<fftwf_complex, AlignedDeleter> dftc_smart(dftc);
      std::unique_ptr<fftwf_complex, AlignedDeleter> dftc2_smart(dftc2);

      for (int ct = 0; ct < nnpoints; ct++) {
        for (int z = 0; z < ep.tbsize; z++) {
          DSFrame src = (*fetch_frame)(npts[ct].fn + z);
          int plane = npts[ct].b;
          auto src_ptr = src.SrcPointers[plane];
          const int stride_elements = src.StrideBytes[plane] / ep.vi_bytesPerSample;

          if (ep.vi_bytesPerSample == 1) {
            const uint8_t * srcp = src_ptr + stride_elements * npts[ct].y + npts[ct].x;
            proc0_c(srcp, hw2 + ep.barea * z, dftr + ep.barea * z, stride_elements, ep.sbsize, ep.divisor);
          } else if (ep.vi_bytesPerSample == 2) {
            const uint16_t * srcp = reinterpret_cast<const uint16_t *>(src_ptr) + stride_elements * npts[ct].y + npts[ct].x;
            proc0_c(srcp, hw2 + ep.barea * z, dftr + ep.barea * z, stride_elements, ep.sbsize, ep.divisor);
          } else {
            const float * srcp = reinterpret_cast<const float *>(src_ptr) + stride_elements * npts[ct].y + npts[ct].x;
            proc0_c(srcp, hw2 + ep.barea * z, dftr + ep.barea * z, stride_elements, ep.sbsize, ep.divisor);
          }
        }

        fft.fftwf_execute_dft_r2c(ep.ft, dftr, dftc);

        if (ep.zmean)
          removeMean_c(reinterpret_cast<float *>(dftc), reinterpret_cast<const float *>(dftgc2), ep.ccnt2, reinterpret_cast<float *>(dftc2));

        for (int h = 0; h < ep.ccnt2; h += 2) {
          const float psd = reinterpret_cast<float *>(dftc)[h] * reinterpret_cast<float *>(dftc)[h] + reinterpret_cast<float *>(dftc)[h + 1] * reinterpret_cast<float *>(dftc)[h + 1];
          ep.sigmas[h] += psd;
          ep.sigmas[h + 1] += psd;
        }

      }
      dftc_smart.reset();
      dftc2_smart.reset();

      hw2_smart.reset();
      dftr_smart.reset();
      dftgc2_smart.reset();
      npts_smart.reset();

      const float scale = 1.0f / nnpoints;
      for (int h = 0; h < ep.ccnt2; h++)
        ep.sigmas[h] *= scale * (wscale2 / wscale) * alpha;
    }
  }

  std::vector<int> RequestReferenceFrames(int n) const override
  {
    std::vector<int> req;
    if (ep.tbsize == 1) {
      req.push_back(n);
    } else {
      const int start = std::max(n - ep.tbsize / 2, 0);
      const int stop = std::min(n + ep.tbsize / 2, in_vi.Frames - 1);
      for (int i = start; i <= stop; i++)
        req.push_back(i);
    }
    return req;
  }

  DSFrame GetFrame(int n, std::unordered_map<int, DSFrame> in_frames) override
  {
    unsigned int thread_id;

    {
      std::lock_guard<std::mutex> lock(thread_check_mutex);
      // Find empty slot
      auto it = std::find(thread_id_store.begin(), thread_id_store.end(), 0);
      thread_id = static_cast<int>(std::distance(thread_id_store.begin(), it));
      if (it == thread_id_store.end()) {
        thread_id_store.push_back(1);

        while (ep.ebuff.size() <= thread_id)
          ep.ebuff.push_back((float *)_aligned_malloc(sizeof(float) * ep.eStride[0] * ep.padHeight[0], FRAME_ALIGN));
        while (ep.dftr.size() <= thread_id)
          ep.dftr.push_back((float *)_aligned_malloc(sizeof(float) * (((ep.bvolume + 7) | 15) + 1) * ep.threads, FRAME_ALIGN));
        while (ep.dftc.size() <= thread_id)
          ep.dftc.push_back((fftwf_complex*)_aligned_malloc(sizeof(fftwf_complex) * (((ep.ccnt + 7) | 15) + 1) * ep.threads, FRAME_ALIGN));
        while (ep.dftc2.size() <= thread_id)
          ep.dftc2.push_back((fftwf_complex*)_aligned_malloc(sizeof(fftwf_complex) * (((ep.ccnt + 7) | 15) + 1) * ep.threads, FRAME_ALIGN));
        while (ep.process[0] == 3 && ep.pad[0].size() <= thread_id)
          ep.pad[0].push_back((unsigned char *)_aligned_malloc(ep.padBlockSize[0] * ep.tbsize, FRAME_ALIGN));
        while (ep.process[1] == 3 && ep.pad[1].size() <= thread_id)
          ep.pad[1].push_back((unsigned char *)_aligned_malloc(ep.padBlockSize[1] * ep.tbsize, FRAME_ALIGN));
        while (ep.process[2] == 3 && ep.pad[2].size() <= thread_id)
          ep.pad[2].push_back((unsigned char *)_aligned_malloc(ep.padBlockSize[2] * ep.tbsize, FRAME_ALIGN));
        if (ep.dither > 0) {
          while (ep.rngs.size() <= thread_id)
            ep.rngs.push_back(std::make_unique<MTRand>());
          while (ep.d_buffs.size() <= thread_id)
            ep.d_buffs.push_back((float *)_aligned_malloc(sizeof(float) * 2 * ep.vi_width, FRAME_ALIGN));
        }
      }
      else
        thread_id_store[thread_id] = 1;
    }

    auto src0 = in_frames[n];
    auto dst = src0.Create(false);

    if (ep.tbsize == 1) {
      for (int p = 0; p < ep.vi_numPlanes; p++) {
        bool chroma = in_vi.Format.IsFamilyYUV && p > 0 && p < 3;
        auto height = in_vi.Height;
        auto width = in_vi.Width;
        auto src_stride = src0.StrideBytes[p];
        auto src_ptr = src0.SrcPointers[p];
        auto dst_stride = dst.StrideBytes[p];
        auto dst_ptr = dst.DstPointers[p];

        if (chroma) {
          height >>= in_vi.Format.SSH;
          width >>= in_vi.Format.SSW;
        }

        if (ep.process[p] == 3) {
          auto pad = ep.pad[p][thread_id];
          ep.copyPad(p, src_ptr, src_stride, pad, &ep);
          ep.func_0(thread_id, p, pad, dst_ptr, dst_stride, &ep);
        }
        else if (ep.process[p] == 2) {
          framecpy(dst_ptr, dst_stride, src_ptr, src_stride, width * in_vi.Format.BytesPerSample, height);
        }
      }

      {
        std::lock_guard<std::mutex> lock(thread_check_mutex);
        thread_id_store[thread_id] = 0;
      }
      
      return dst;
    }

    const int pos = ep.tbsize / 2;

    for (int p = 0; p < ep.vi_numPlanes; p++) {
      bool chroma = in_vi.Format.IsFamilyYUV && p > 0 && p < 3;
      auto height = in_vi.Height;
      auto width = in_vi.Width;
      auto src0_stride = src0.StrideBytes[p];
      auto src0_ptr = src0.SrcPointers[p];
      auto dst_stride = dst.StrideBytes[p];
      auto dst_ptr = dst.DstPointers[p];

      if (chroma) {
        height >>= in_vi.Format.SSH;
        width >>= in_vi.Format.SSW;
      }

      if (ep.process[p] == 3) {
        auto pad0 = ep.pad[p][thread_id];
        for (int i = 0; i < ep.tbsize; i++) {
          int fn = i + n - pos;
          int fn_real = std::min(std::max(fn, 0), in_vi.Frames - 1);

          auto src_stride = in_frames[fn_real].StrideBytes[p];
          auto src_ptr = in_frames[fn_real].SrcPointers[p];

          auto pad = pad0 + ep.padBlockSize[p] * i;
          ep.copyPad(p, src_ptr, src_stride, pad, &ep);
        }
        ep.func_1(thread_id, p, pad0, dst_ptr, dst_stride, pos, &ep);
      }
      else if (ep.process[p] == 2) {
        framecpy(dst_ptr, dst_stride, src0_ptr, src0_stride, width * in_vi.Format.BytesPerSample, height);
      }
    }

    thread_id_store[thread_id] = false;
    return dst;
  }

  void framecpy(unsigned char * dst_ptr, int dst_stride, const unsigned char * src_ptr, int src_stride, int width_byte, int height) {
    if (src_stride == dst_stride) {
      memcpy(dst_ptr, src_ptr, dst_stride * height);
      return;
    }
    for (int h = 0; h < height; h++)
    {
      memcpy(dst_ptr, src_ptr, width_byte);
      dst_ptr += dst_stride;
      src_ptr += src_stride;
    }
  }

  template<typename T>
  void parse_array_string(std::string str, std::vector<T>& sdata) {
    sdata.clear();
    std::stringstream ss(str);
    ss.imbue(std::locale::classic());
    float tmpf;
    while (!ss.eof()) {
      do {
        char ch = ss.peek();
        if (ch == ':' || ch == ' ' || ch == ',' )
          ss.ignore(1);
        else
          break;
      } while (!ss.eof());

      if (ss.eof())
        break;
      ss >> tmpf;
      if (ss.fail())
        throw strdup((std::string{ "Unable to parse string: " } + ss.str()).c_str());
      sdata.push_back(tmpf);
    }
  }

  ~DFTTest() {
    _aligned_free(ep.hw);
    _aligned_free(ep.dftgc);
    _aligned_free(ep.sigmas);
    _aligned_free(ep.sigmas2);
    _aligned_free(ep.pmins);
    _aligned_free(ep.pmaxs);

    for (auto &&buf : ep.ebuff)
      _aligned_free(buf);
    for (auto &&buf : ep.dftr)
      _aligned_free(buf);
    for (auto &&buf : ep.dftc)
      _aligned_free(buf);
    for (auto &&buf : ep.dftc2)
      _aligned_free(buf);
    for (auto &&buf : ep.pad[0])
      _aligned_free(buf);
    for (auto &&buf : ep.pad[1])
      _aligned_free(buf);
    for (auto &&buf : ep.pad[2])
      _aligned_free(buf);
    for (auto&& buf : ep.d_buffs)
      _aligned_free(buf);

    if (fft.library) {
        GlobalLockGuard lock(this->avs_env, "fftw", this->has_at_least_v12);

      if (ep.ft)
        fft.fftwf_destroy_plan(ep.ft);
      if (ep.fti)
        fft.fftwf_destroy_plan(ep.fti);
      fft.free();
    }
  }
};


namespace Plugin {
  const char* Identifier = "in.7086.neo_dfttest";
  const char* Namespace = "neo_dfttest";
  const char* Description = "Neo DFTTest Deband Filter " PLUGIN_VERSION " - 2D/3D frequency domain denoiser";
}

std::vector<register_vsfilter_proc> RegisterVSFilters()
{
  return std::vector<register_vsfilter_proc> { VSInterface::RegisterFilter<DFTTest> };
}

std::vector<register_avsfilter_proc> RegisterAVSFilters()
{
  return std::vector<register_avsfilter_proc> { AVSInterface::RegisterFilter<DFTTest> };
}
