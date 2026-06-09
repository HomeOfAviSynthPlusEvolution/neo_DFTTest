#include <avisynth.h>
#include <vapoursynth/VapourSynth4.h>

#include <dualsynth/avisynth/video_bridge.hpp>
#include <dualsynth/vapoursynth/video_bridge.hpp>
#include <dualsynth/video_bridge.hpp>

#include "plugin/dfttest_filter.hpp"

#include <string>

#if defined(_WIN32)
#define NEO_DFTTEST_AVS_PLUGIN_EXPORT extern "C" __declspec(dllexport)
#elif defined(__clang__) || defined(__GNUC__)
#define NEO_DFTTEST_AVS_PLUGIN_EXPORT extern "C" __attribute__((visibility("default")))
#else
#define NEO_DFTTEST_AVS_PLUGIN_EXPORT extern "C"
#endif

const AVS_Linkage* AVS_linkage = nullptr;

namespace {

using DFTTestBridge = neo_dfttest::DFTTestBridge;

const char* vs_signature() {
  static const std::string signature = [] {
    const auto generated = ds::make_vapoursynth_signature(DFTTestBridge::descriptor());
    if (!generated.has_value()) {
      return std::string{"clip:vnode;"};
    }
    return generated.value();
  }();
  return signature.c_str();
}

const char* avs_signature() {
  static const std::string signature = [] {
    const auto generated = ds::make_avisynth_signature(DFTTestBridge::descriptor());
    if (!generated.has_value()) {
      return std::string{"c"};
    }
    return generated.value();
  }();
  return signature.c_str();
}

void VS_CC create_vapoursynth_dfttest(
  const VSMap* in,
  VSMap* out,
  void*,
  VSCore* core,
  const VSAPI* vsapi
) {
  ds::vapoursynth::create_video_filter_bridge<DFTTestBridge>(in, out, core, vsapi);
}

AVSValue __cdecl create_avisynth_dfttest(AVSValue args, void*, IScriptEnvironment* env) {
  return ds::avisynth::create_video_filter_bridge<DFTTestBridge>(args, env);
}

} // namespace

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
  vspapi->configPlugin(
    neo_dfttest::Plugin::Identifier,
    neo_dfttest::Plugin::Namespace,
    neo_dfttest::Plugin::Description,
    VS_MAKE_VERSION(9, 6),
    VAPOURSYNTH_API_VERSION,
    0,
    plugin
  );

  vspapi->registerFunction(
    DFTTestBridge::vs_name,
    vs_signature(),
    "clip:vnode;",
    create_vapoursynth_dfttest,
    nullptr,
    plugin
  );
}

NEO_DFTTEST_AVS_PLUGIN_EXPORT const char* __stdcall AvisynthPluginInit3(
  IScriptEnvironment* env,
  const AVS_Linkage* const vectors
) {
  AVS_linkage = vectors;
  env->AddFunction(
    DFTTestBridge::avs_name,
    avs_signature(),
    create_avisynth_dfttest,
    nullptr
  );
  ds::avisynth::set_video_filter_mt_mode<DFTTestBridge>(env);
  return neo_dfttest::Plugin::Description;
}
