#pragma once

#include <dualsynth/global_lock.hpp>
#include <dualsynth/video_filter.hpp>

#include <memory>

namespace neo_dfttest {

class DftTestEngine {
public:
  DftTestEngine();
  ~DftTestEngine();

  DftTestEngine(DftTestEngine&&) noexcept;
  DftTestEngine& operator=(DftTestEngine&&) noexcept;
  DftTestEngine(const DftTestEngine&) = delete;
  DftTestEngine& operator=(const DftTestEngine&) = delete;

  void initialize(
    const ds::VideoInputInfo& input,
    const ds::ParamValues* params,
    ds::HostGlobalLockCallbacks host_locks
  );
  void request_frames(ds::VideoRequestContext& context) const;
  void process_frame(int n, ds::VideoFrameProvider& provider, ds::MutableVideoFrameView dst);
  int cache_hints(int cachehints, int frame_range, int default_response);

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace neo_dfttest
