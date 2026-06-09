#pragma once

#include "dft_common.h"

#include <dualsynth/frame.hpp>

#include <cstddef>
#include <cstring>

namespace neo_dfttest::engine {

struct ReadPlane {
  const unsigned char* data {};
  int stride_bytes {};
  int width {};
  int height {};
  int row_bytes {};

  int stride_elements(int bytes_per_sample) const {
    return stride_bytes / bytes_per_sample;
  }

  template <class T>
  const T* typed_at(int x, int y, int bytes_per_sample) const {
    return reinterpret_cast<const T*>(data) + stride_elements(bytes_per_sample) * y + x;
  }
};

struct WritePlane {
  unsigned char* data {};
  int stride_bytes {};
  int width {};
  int height {};
  int row_bytes {};
};

inline ReadPlane read_plane(const ds::VideoFrameView& frame, int plane, const DFTTestData& state) {
  const ds::PlaneView& view = frame.plane(plane);
  const int width = state.planes.width[plane];
  const int height = state.planes.height[plane];
  return ReadPlane{
    static_cast<const unsigned char*>(view.data),
    static_cast<int>(view.stride_bytes),
    width,
    height,
    width * state.format.bytes_per_sample
  };
}

inline WritePlane write_plane(ds::MutableVideoFrameView& frame, int plane, const DFTTestData& state) {
  const ds::MutablePlaneView& view = frame.plane(plane);
  const int width = state.planes.width[plane];
  const int height = state.planes.height[plane];
  return WritePlane{
    static_cast<unsigned char*>(view.data),
    static_cast<int>(view.stride_bytes),
    width,
    height,
    width * state.format.bytes_per_sample
  };
}

inline void copy_plane_rows(const ReadPlane& src, const WritePlane& dst) {
  if (src.stride_bytes == dst.stride_bytes) {
    std::memcpy(
      dst.data,
      src.data,
      static_cast<std::size_t>(dst.stride_bytes) * static_cast<std::size_t>(dst.height)
    );
    return;
  }

  const unsigned char* src_row = src.data;
  unsigned char* dst_row = dst.data;
  for (int y = 0; y < dst.height; ++y) {
    std::memcpy(dst_row, src_row, static_cast<std::size_t>(dst.row_bytes));
    src_row += src.stride_bytes;
    dst_row += dst.stride_bytes;
  }
}

} // namespace neo_dfttest::engine
