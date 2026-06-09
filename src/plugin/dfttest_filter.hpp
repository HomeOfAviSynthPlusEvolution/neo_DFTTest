#pragma once

#include <dualsynth/avisynth/video_bridge.hpp>
#include <dualsynth/format.hpp>
#include <dualsynth/param.hpp>
#include <dualsynth/video_bridge.hpp>
#include <dualsynth/video_filter.hpp>

#include "src/version.hpp"

#include <memory>

namespace neo_dfttest {

class DftTestEngine;

struct DFTTestCore {
  static constexpr const char* name = "DFTTest";
  static constexpr int input_count = 1;
  static constexpr ds::OutputOrigin output_origin = ds::OutputOrigin::fresh();

  struct State {
    State();
    explicit State(std::unique_ptr<DftTestEngine> engine);
    ~State();

    State(State&&) noexcept;
    State& operator=(State&&) noexcept;
    State(const State&) = delete;
    State& operator=(const State&) = delete;

    std::unique_ptr<DftTestEngine> engine;
  };

  static ds::Result<ds::VideoInitStateResult<State>> init(ds::VideoInitContext& context);
  static ds::Result<ds::VideoRequestResult> request(ds::VideoRequestContext& context);
  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& context);
  static int cache_hints(ds::VideoCacheHintsContext& context);
};

struct DFTTestBridge : ds::SingleInputVideoBridgeDefaults<DFTTestCore> {
  static constexpr const char* vs_name = "DFTTest";
  static constexpr const char* avs_name = "neo_dfttest";
  static constexpr const char* vs_signature = "";
  static constexpr const char* avs_signature = "";
  static constexpr const char* missing_input_error = "neo_dfttest: missing required video clip";
  static constexpr const char* vs_format_error =
    "neo_dfttest: only constant 8-16 bit integer and 32 bit float video is supported";
  static constexpr const char* avs_format_error =
    "neo_dfttest: only constant 8-16 bit integer and 32 bit float video is supported";
  static constexpr ds::avisynth::MtMode avs_mt_mode = ds::avisynth::MtMode::NiceFilter;

  static bool accepts_video_format(ds::VideoFormat format);
  static ds::FilterDescriptor descriptor();
};

namespace Plugin {
inline constexpr const char* Identifier = "in.7086.neo_dfttest";
inline constexpr const char* Namespace = "neo_dfttest";
inline constexpr const char* Description =
  "Neo DFTTest Deband Filter " PLUGIN_VERSION " - 2D/3D frequency domain denoiser";
} // namespace Plugin

} // namespace neo_dfttest
