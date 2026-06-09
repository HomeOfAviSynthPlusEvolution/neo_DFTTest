#include "engine/dfttest_config.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace neo_dfttest {
namespace {

template <class T>
T unwrap(ds::Result<T> result) {
  if (!result.has_value()) {
    throw std::runtime_error(result.error().message);
  }
  return std::move(result.value());
}

bool has_param(const ds::ParamValues& params, const std::string& name) {
  return std::ranges::any_of(params.entries, [&](const ds::ParamEntry& entry) {
    return entry.name == name;
  });
}

int read_int(const ds::ParamValues& params, const std::string& name, int default_value) {
  return unwrap(params.get_int(name, default_value));
}

float read_float(const ds::ParamValues& params, const std::string& name, float default_value) {
  return static_cast<float>(unwrap(params.get_double(name, default_value)));
}

bool read_bool(const ds::ParamValues& params, const std::string& name, bool default_value) {
  return unwrap(params.get_bool(name, default_value));
}

std::vector<int> read_int_array(
  const ds::ParamValues& params,
  const std::string& name,
  std::vector<std::int64_t> default_value = {}
) {
  const auto values = unwrap(params.get_int_array(name, std::move(default_value)));
  std::vector<int> output;
  output.reserve(values.size());
  for (const std::int64_t value : values) {
    if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max()) {
      throw std::runtime_error("parameter '" + name + "' contains an integer outside int range");
    }
    output.push_back(static_cast<int>(value));
  }
  return output;
}

std::vector<float> read_float_array(
  const ds::ParamValues& params,
  const std::string& name,
  std::vector<double> default_value = {}
) {
  const auto values = unwrap(params.get_double_array(name, std::move(default_value)));
  std::vector<float> output;
  output.reserve(values.size());
  for (const double value : values) {
    output.push_back(static_cast<float>(value));
  }
  return output;
}

} // namespace

DfttestConfig DfttestConfig::read(const ds::ParamValues& values, DFTTestData& state) {
  DfttestConfig config{};

  config.ftype = read_int(values, "ftype", 0);
  config.sigma = read_float(values, "sigma", 8.0f);
  config.sigma2 = read_float(values, "sigma2", 8.0f);
  config.pmin = read_float(values, "pmin", 0.0f);
  config.pmax = read_float(values, "pmax", 500.0f);
  config.smode = read_int(values, "smode", 1);
  config.tmode = read_int(values, "tmode", 0);
  config.opt = read_int(values, "opt", 0);

  state.block.spatial_size = read_int(values, "sbsize", state.block.spatial_size);
  state.block.spatial_overlap = read_int(values, "sosize", state.block.spatial_overlap);
  state.block.temporal_size = read_int(values, "tbsize", state.block.temporal_size);
  state.block.temporal_overlap = read_int(values, "tosize", state.block.temporal_overlap);
  state.block.spatial_window = read_int(values, "swin", state.block.spatial_window);
  state.block.temporal_window = read_int(values, "twin", state.block.temporal_window);
  state.block.spatial_beta = read_float(values, "sbeta", state.block.spatial_beta);
  state.block.temporal_beta = read_float(values, "tbeta", state.block.temporal_beta);
  state.block.zero_mean = read_bool(values, "zmean", state.block.zero_mean);
  state.block.f0_beta = read_float(values, "f0beta", state.block.f0_beta);
  state.block.worker_threads = read_int(values, "threads", state.block.worker_threads);
  state.block.dither_mode = read_int(values, "dither", state.block.dither_mode);
  config.fft_threads = read_int(values, "fft_threads", config.fft_threads);
  if (config.fft_threads < 1) {
    config.fft_threads = 1;
  }

  if (config.smode == 0) {
    state.block.spatial_overlap = 0;
  }
  if (config.tmode == 0) {
    state.block.temporal_overlap = 0;
  }

  config.nlocation = read_int_array(values, "nlocation");
  config.alpha = read_float(values, "alpha", config.ftype == 0 ? 5.0f : 7.0f);
  config.slocation = read_float_array(values, "slocation");
  config.ssx = read_float_array(values, "ssx");
  config.ssy = read_float_array(values, "ssy");
  config.sst = read_float_array(values, "sst");
  config.ssystem = read_int(values, "ssystem", 0);

  if (state.block.worker_threads <= 0) {
    state.block.worker_threads = 4;
  }
  if (state.block.worker_threads > 16) {
    state.block.worker_threads = 16;
  }

  return config;
}

void DfttestConfig::configure_planes(const ds::ParamValues& values, DFTTestData& state) const {
  std::fill(std::begin(state.planes.process), std::end(state.planes.process), 2);

  if (has_param(values, "planes")) {
    const auto planes = read_int_array(values, "planes");
    for (const int plane : planes) {
      if (plane < 0 || plane >= state.format.num_planes) {
        throw std::runtime_error("plane index out of range");
      }
      state.planes.process[plane] = 3;
    }
    return;
  }

  state.planes.process[0] = 3;
  state.planes.process[1] = 3;
  state.planes.process[2] = 3;
  state.planes.process[3] = 2;
  state.planes.process[0] = read_int(values, "y", state.planes.process[0]);
  state.planes.process[1] = read_int(values, "u", state.planes.process[1]);
  state.planes.process[2] = read_int(values, "v", state.planes.process[2]);
  state.planes.process[3] = read_int(values, "a", state.planes.process[3]);
}

void DfttestConfig::validate(const ds::VideoInputInfo& input, const DFTTestData& state) const {
  if (input.width <= 0 || input.height <= 0) {
    throw std::runtime_error("only constant format input supported");
  }

  if (
    (state.format.integer && state.format.bits_per_sample > 16) ||
    (!state.format.integer && state.format.bits_per_sample != 32)
  ) {
    throw std::runtime_error("only 8-16 bit integer and 32 bit float input supported");
  }

  if (ftype < 0 || ftype > 4) {
    throw std::runtime_error("ftype must be 0, 1, 2, 3, or 4");
  }
  if (state.block.spatial_size < 1) {
    throw std::runtime_error("sbsize must be greater than or equal to 1");
  }
  if (smode < 0 || smode > 1) {
    throw std::runtime_error("smode must be 0 or 1");
  }
  if (smode == 0 && !(state.block.spatial_size & 1)) {
    throw std::runtime_error("sbsize must be odd when using smode=0");
  }
  if (state.block.spatial_overlap < 0 || state.block.spatial_overlap >= state.block.spatial_size) {
    throw std::runtime_error("sosize must be between 0 and sbsize-1 (inclusive)");
  }
  if (
    state.block.spatial_overlap > state.block.spatial_size / 2 &&
    state.block.spatial_size % (state.block.spatial_size - state.block.spatial_overlap) != 0
  ) {
    throw std::runtime_error(
      "spatial overlap greater than 50% requires that sbsize-sosize is a divisor of sbsize"
    );
  }
  if (state.block.temporal_size < 1 || state.block.temporal_size > 15) {
    throw std::runtime_error("tbsize must be between 1 and 15 (inclusive)");
  }
  if (tmode != 0) {
    throw std::runtime_error("tmode must be 0. tmode=1 is not implemented");
  }
  if (tmode == 0 && !(state.block.temporal_size & 1)) {
    throw std::runtime_error("tbsize must be odd when using tmode=0");
  }
  if (state.block.temporal_overlap < 0 || state.block.temporal_overlap >= state.block.temporal_size) {
    throw std::runtime_error("tosize must be between 0 and tbsize-1 (inclusive)");
  }
  if (
    state.block.temporal_overlap > state.block.temporal_size / 2 &&
    state.block.temporal_size % (state.block.temporal_size - state.block.temporal_overlap) != 0
  ) {
    throw std::runtime_error(
      "temporal overlap greater than 50% requires that tbsize-tosize is a divisor of tbsize"
    );
  }
  if (state.block.temporal_size > input.num_frames) {
    throw std::runtime_error("tbsize must be less than or equal to the number of frames in the clip");
  }
  if (state.block.spatial_window < 0 || state.block.spatial_window > 11) {
    throw std::runtime_error("swin must be between 0 and 11 (inclusive)");
  }
  if (state.block.temporal_window < 0 || state.block.temporal_window > 11) {
    throw std::runtime_error("twin must be between 0 and 11 (inclusive)");
  }
  if (nlocation.size() & 3U) {
    throw std::runtime_error("the number of elements in nlocation must be a multiple of 4");
  }
  if (alpha <= 0.0f) {
    throw std::runtime_error("alpha must be greater than 0.0");
  }
  if (slocation.size() & 1U) {
    throw std::runtime_error("the number of elements in slocation must be even");
  }
  if (ssx.size() & 1U) {
    throw std::runtime_error("the number of elements in ssx must be even");
  }
  if (ssy.size() & 1U) {
    throw std::runtime_error("the number of elements in ssy must be even");
  }
  if (sst.size() & 1U) {
    throw std::runtime_error("the number of elements in sst must be even");
  }
  if (ssystem < 0 || ssystem > 1) {
    throw std::runtime_error("ssystem must be 0 or 1");
  }
  if (opt != 8 && (opt < 0 || opt > 3)) {
    throw std::runtime_error("opt must be 0, 1, 2, 3, or 8");
  }
}

} // namespace neo_dfttest
