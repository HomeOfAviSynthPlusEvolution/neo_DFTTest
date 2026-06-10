#pragma once

#include "fft/fft_backend.hpp"

namespace neo_dfttest {

enum class DftTransferDirection {
  none,
  upload,
  download,
  bidirectional,
};

struct DftBufferMemoryPlan {
  fft::MemoryDomain domain {fft::MemoryDomain::host};
  DftTransferDirection transfer {DftTransferDirection::none};
};

struct DftMemoryPlan {
  DftBufferMemoryPlan source_frame;
  DftBufferMemoryPlan padded_frame;
  DftBufferMemoryPlan real_scratch;
  DftBufferMemoryPlan complex_scratch;
  DftBufferMemoryPlan coefficients;
  DftBufferMemoryPlan accumulation;
  DftBufferMemoryPlan output_frame;
  bool host_fallback {false};
};

inline DftBufferMemoryPlan dft_host_buffer_plan() noexcept {
  return DftBufferMemoryPlan{fft::MemoryDomain::host, DftTransferDirection::none};
}

inline DftMemoryPlan dft_host_memory_plan() noexcept {
  const DftBufferMemoryPlan host = dft_host_buffer_plan();
  return DftMemoryPlan{host, host, host, host, host, host, host, false};
}

inline DftMemoryPlan dft_vulkan_memory_plan(bool host_fallback) noexcept {
  const DftBufferMemoryPlan host {fft::MemoryDomain::host, DftTransferDirection::none};
  const DftBufferMemoryPlan upload {fft::MemoryDomain::device, DftTransferDirection::upload};
  const DftBufferMemoryPlan device {fft::MemoryDomain::device, DftTransferDirection::none};
  const DftBufferMemoryPlan download {fft::MemoryDomain::host, DftTransferDirection::download};
  return DftMemoryPlan{
    host,
    upload,
    device,
    device,
    device,
    device,
    download,
    host_fallback
  };
}

} // namespace neo_dfttest
