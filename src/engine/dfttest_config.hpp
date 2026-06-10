#pragma once

#include "dft_common.h"

#include <dualsynth/param.hpp>
#include <dualsynth/video_filter.hpp>

#include <string>
#include <vector>

namespace neo_dfttest {

struct DfttestConfig {
  std::string fft_backend = "fftw";
  int fft_threads = 1;
  bool worker_threads_auto = true;
  bool fft_threads_auto = true;
  bool dither_seed_set = false;
  std::uint32_t dither_seed = 0;
  int ftype = 0;
  int smode = 1;
  int tmode = 0;
  int opt = 0;
  int ssystem = 0;
  float sigma = 8.0f;
  float sigma2 = 8.0f;
  float pmin = 0.0f;
  float pmax = 500.0f;
  float alpha = 5.0f;
  std::vector<int> nlocation;
  std::vector<float> slocation;
  std::vector<float> ssx;
  std::vector<float> ssy;
  std::vector<float> sst;

  static DfttestConfig read(const ds::ParamValues& values, DFTTestData& state);

  void configure_planes(const ds::ParamValues& values, DFTTestData& state) const;
  void validate(const ds::VideoInputInfo& input, const DFTTestData& state) const;
};

} // namespace neo_dfttest
