#include "plugin/dfttest_filter.hpp"

#include "engine/dfttest_engine.hpp"

#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace neo_dfttest {

namespace {

ds::Error invalid_argument(std::string message) {
  return ds::Error{ds::ErrorCode::InvalidArgument, std::move(message)};
}

} // namespace

DFTTestCore::State::State() = default;

DFTTestCore::State::State(std::unique_ptr<DftTestEngine> input_engine)
  : engine(std::move(input_engine)) {}

DFTTestCore::State::~State() = default;

DFTTestCore::State::State(State&&) noexcept = default;

DFTTestCore::State& DFTTestCore::State::operator=(State&&) noexcept = default;

ds::Result<ds::VideoInitStateResult<DFTTestCore::State>> DFTTestCore::init(
  ds::VideoInitContext& context
) {
  try {
    const auto inputs = ds::collect_video_input_infos<DFTTestCore>(context.inputs);
    if (!inputs.has_value()) {
      return ds::Result<ds::VideoInitStateResult<State>>::failure(inputs.error());
    }

    auto engine = std::make_unique<DftTestEngine>();
    engine->initialize(inputs.value()[0], context.params, context.host_global_locks);
    const ds::VideoInputInfo& input = inputs.value()[0];
    return ds::Result<ds::VideoInitStateResult<State>>::success(
      ds::VideoInitStateResult<State>{
        ds::VideoOutputInfo{input.width, input.height, input.num_frames, input.format, input.fps},
        State{std::move(engine)}
      }
    );
  } catch (const char* error) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(invalid_argument(error));
  } catch (const std::exception& error) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(invalid_argument(error.what()));
  } catch (...) {
    return ds::Result<ds::VideoInitStateResult<State>>::failure(
      invalid_argument("neo_dfttest: unhandled initialization error")
    );
  }
}

ds::Result<ds::VideoRequestResult> DFTTestCore::request(ds::VideoRequestContext& context) {
  try {
    context.state<State>().engine->request_frames(context);
    return ds::Result<ds::VideoRequestResult>::success(ds::VideoRequestResult{});
  } catch (const std::exception& error) {
    return ds::Result<ds::VideoRequestResult>::failure(invalid_argument(error.what()));
  }
}

ds::Result<ds::VideoProcessResult> DFTTestCore::process(ds::VideoProcessContext& context) {
  try {
    context.state<State>().engine->process_frame(context.output_frame, context.frames, context.dst);
    return ds::Result<ds::VideoProcessResult>::success(ds::VideoProcessResult{});
  } catch (const char* error) {
    return ds::Result<ds::VideoProcessResult>::failure(invalid_argument(error));
  } catch (const std::exception& error) {
    return ds::Result<ds::VideoProcessResult>::failure(invalid_argument(error.what()));
  } catch (...) {
    return ds::Result<ds::VideoProcessResult>::failure(
      invalid_argument("neo_dfttest: unhandled processing error")
    );
  }
}

int DFTTestCore::cache_hints(ds::VideoCacheHintsContext& context) {
  return context.state<State>().engine->cache_hints(
    context.cachehints,
    context.frame_range,
    context.default_response
  );
}

bool DFTTestBridge::accepts_video_format(ds::VideoFormat format) {
  return format.plane_count >= 1 &&
    format.plane_count <= 4 &&
    (format.sample_format == ds::SampleFormat::UInt8 ||
     format.sample_format == ds::SampleFormat::UInt10 ||
     format.sample_format == ds::SampleFormat::UInt12 ||
     format.sample_format == ds::SampleFormat::UInt14 ||
     format.sample_format == ds::SampleFormat::UInt16 ||
     format.sample_format == ds::SampleFormat::Float32);
}

ds::FilterDescriptor DFTTestBridge::descriptor() {
  return ds::FilterDescriptor{
    "DFTTest",
    std::vector<ds::ParamSpec>{
      ds::ParamSpec{"clip", ds::ParamType::Clip, ds::ParamValue{}, true},
      ds::ParamSpec{"ftype", ds::ParamType::Integer},
      ds::ParamSpec{"sigma", ds::ParamType::Float},
      ds::ParamSpec{"sigma2", ds::ParamType::Float},
      ds::ParamSpec{"pmin", ds::ParamType::Float},
      ds::ParamSpec{"pmax", ds::ParamType::Float},
      ds::ParamSpec{"sbsize", ds::ParamType::Integer},
      ds::ParamSpec{"smode", ds::ParamType::Integer},
      ds::ParamSpec{"sosize", ds::ParamType::Integer},
      ds::ParamSpec{"tbsize", ds::ParamType::Integer},
      ds::ParamSpec{"tmode", ds::ParamType::Integer},
      ds::ParamSpec{"tosize", ds::ParamType::Integer},
      ds::ParamSpec{"swin", ds::ParamType::Integer},
      ds::ParamSpec{"twin", ds::ParamType::Integer},
      ds::ParamSpec{"sbeta", ds::ParamType::Float},
      ds::ParamSpec{"tbeta", ds::ParamType::Float},
      ds::ParamSpec{"zmean", ds::ParamType::Boolean},
      ds::ParamSpec{"f0beta", ds::ParamType::Float},
      ds::ParamSpec{"nlocation", ds::ParamType::Integer, ds::ParamValue{}, false, true},
      ds::ParamSpec{"alpha", ds::ParamType::Float},
      ds::ParamSpec{"slocation", ds::ParamType::Float, ds::ParamValue{}, false, true},
      ds::ParamSpec{"ssx", ds::ParamType::Float, ds::ParamValue{}, false, true},
      ds::ParamSpec{"ssy", ds::ParamType::Float, ds::ParamValue{}, false, true},
      ds::ParamSpec{"sst", ds::ParamType::Float, ds::ParamValue{}, false, true},
      ds::ParamSpec{"ssystem", ds::ParamType::Integer},
      ds::ParamSpec{"dither", ds::ParamType::Integer},
      ds::ParamSpec{"planes", ds::ParamType::Integer, ds::ParamValue{}, false, true, true, false},
      ds::ParamSpec{"y", ds::ParamType::Integer, ds::ParamValue{}, false, false, false, true},
      ds::ParamSpec{"u", ds::ParamType::Integer, ds::ParamValue{}, false, false, false, true},
      ds::ParamSpec{"v", ds::ParamType::Integer, ds::ParamValue{}, false, false, false, true},
      ds::ParamSpec{"a", ds::ParamType::Integer, ds::ParamValue{}, false, false, false, true},
      ds::ParamSpec{"opt", ds::ParamType::Integer},
      ds::ParamSpec{"threads", ds::ParamType::Integer},
      ds::ParamSpec{"fft_threads", ds::ParamType::Integer}
    }
  };
}

} // namespace neo_dfttest
