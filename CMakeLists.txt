cmake_minimum_required (VERSION 3.24)
project (neo_dfttest CXX)
option(NEO_DFTTEST_ENABLE_VKFFT_VULKAN "Enable optional VkFFT Vulkan backend" OFF)

add_library(neo-dfttest SHARED
  src/plugin/plugin_entry.cpp
  src/plugin/dfttest_filter.cpp
  src/plugin/dfttest_filter.hpp
  src/version.rc
  src/core_dispatch.cpp
  src/core_window.cpp
  src/core.h
  src/core_scalar.cpp
  src/core_simd_hwy.cpp
  src/dft_common.h
  src/engine/dfttest_engine.cpp
  src/engine/dfttest_engine.hpp
  src/engine/dfttest_config.cpp
  src/engine/dfttest_config.hpp
  src/executor/cpu_dft_executor.cpp
  src/executor/dft_executor.hpp
  src/executor/dft_fft_operation.hpp
  src/executor/dft_memory_plan.hpp
  src/fft/fft_backend.hpp
  src/fft/dynamic_library.cpp
  src/fft/dynamic_library.hpp
  src/fft/fftw_abi.hpp
  src/fft/fftw_backend.cpp
  src/fft/fftw_runtime.cpp
  src/fft/fftw_runtime.hpp
  src/fft/pocketfft_backend.cpp
  src/fft/vkfft_vulkan_buffer.hpp
  src/fft/vkfft_vulkan_runtime.hpp
  src/engine/pixel_plane.hpp
  src/memory/aligned_buffer.hpp
  src/pipeline/dft_cpu_pipeline.hpp
  src/workspace/dft_workspace.hpp
)
set_property(TARGET neo-dfttest PROPERTY CXX_STANDARD 23)

add_executable(neo-dfttest_avs_frame_dump
  tests/tools/avs_frame_dump.cpp
)
set_property(TARGET neo-dfttest_avs_frame_dump PROPERTY CXX_STANDARD 17)
target_include_directories(neo-dfttest_avs_frame_dump PRIVATE include/avisynth)

if (NOT WIN32)
  target_link_libraries(neo-dfttest_avs_frame_dump dl)
endif()

target_compile_definitions(neo-dfttest PRIVATE _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)

include(FetchContent)

if(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)
  find_package(Vulkan REQUIRED)
  FetchContent_Declare(
    vkfft
    GIT_REPOSITORY https://github.com/DTolm/VkFFT.git
    GIT_TAG v1.3.4
  )
  FetchContent_MakeAvailable(vkfft)

  set(ENABLE_OPT 0 CACHE BOOL "" FORCE)
  FetchContent_Declare(
    glslang-main
    GIT_REPOSITORY https://github.com/KhronosGroup/glslang.git
    GIT_TAG 12.3.1
  )
  FetchContent_MakeAvailable(glslang-main)

  target_sources(neo-dfttest PRIVATE
    src/executor/vulkan_dft_executor.cpp
    src/fft/vkfft_vulkan_buffer.cpp
    src/fft/vkfft_vulkan_backend.cpp
    src/vulkan/vulkan_compute.cpp
    src/vulkan/vulkan_compute.hpp
  )
  set(NEO_DFTTEST_GENERATED_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/neo_dfttest")
  set(NEO_DFTTEST_LOAD_WINDOW_SPV "${NEO_DFTTEST_GENERATED_DIR}/load_window_comp_spv.inc")
  set(NEO_DFTTEST_FREQUENCY_OPS_SPV "${NEO_DFTTEST_GENERATED_DIR}/frequency_ops_comp_spv.inc")
  set(NEO_DFTTEST_ACCUMULATE_INVERSE_SPV "${NEO_DFTTEST_GENERATED_DIR}/accumulate_inverse_comp_spv.inc")
  set(NEO_DFTTEST_WRITE_OUTPUT_SPV "${NEO_DFTTEST_GENERATED_DIR}/write_output_comp_spv.inc")
  add_custom_command(
    OUTPUT "${NEO_DFTTEST_LOAD_WINDOW_SPV}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${NEO_DFTTEST_GENERATED_DIR}"
    COMMAND $<TARGET_FILE:glslang-standalone> --quiet -V -x -S comp -o "${NEO_DFTTEST_LOAD_WINDOW_SPV}" "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/load_window.comp"
    DEPENDS glslang-standalone "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/load_window.comp"
    VERBATIM
  )
  add_custom_command(
    OUTPUT "${NEO_DFTTEST_FREQUENCY_OPS_SPV}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${NEO_DFTTEST_GENERATED_DIR}"
    COMMAND $<TARGET_FILE:glslang-standalone> --quiet -V -x -S comp -o "${NEO_DFTTEST_FREQUENCY_OPS_SPV}" "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/frequency_ops.comp"
    DEPENDS glslang-standalone "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/frequency_ops.comp"
    VERBATIM
  )
  add_custom_command(
    OUTPUT "${NEO_DFTTEST_ACCUMULATE_INVERSE_SPV}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${NEO_DFTTEST_GENERATED_DIR}"
    COMMAND $<TARGET_FILE:glslang-standalone> --quiet -V -x -S comp -o "${NEO_DFTTEST_ACCUMULATE_INVERSE_SPV}" "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/accumulate_inverse.comp"
    DEPENDS glslang-standalone "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/accumulate_inverse.comp"
    VERBATIM
  )
  add_custom_command(
    OUTPUT "${NEO_DFTTEST_WRITE_OUTPUT_SPV}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${NEO_DFTTEST_GENERATED_DIR}"
    COMMAND $<TARGET_FILE:glslang-standalone> --quiet -V -x -S comp -o "${NEO_DFTTEST_WRITE_OUTPUT_SPV}" "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/write_output.comp"
    DEPENDS glslang-standalone "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders/write_output.comp"
    VERBATIM
  )
  add_custom_target(
    neo-dfttest-vulkan-shaders
    DEPENDS
      "${NEO_DFTTEST_LOAD_WINDOW_SPV}"
      "${NEO_DFTTEST_FREQUENCY_OPS_SPV}"
      "${NEO_DFTTEST_ACCUMULATE_INVERSE_SPV}"
      "${NEO_DFTTEST_WRITE_OUTPUT_SPV}"
  )
  add_dependencies(neo-dfttest neo-dfttest-vulkan-shaders)
  target_include_directories(neo-dfttest PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
    ${vkfft_SOURCE_DIR}/vkFFT
    ${glslang-main_SOURCE_DIR}/glslang/Include
  )
  target_compile_definitions(neo-dfttest PRIVATE
    NEO_DFTTEST_ENABLE_VKFFT_VULKAN=1
    VKFFT_BACKEND=0
    VK_API_VERSION=11
  )
  target_link_libraries(neo-dfttest Vulkan::Vulkan SPIRV glslang)
endif()

FetchContent_Declare(
  highway
  GIT_REPOSITORY https://github.com/google/highway.git
  GIT_TAG 1.3.0
)
set(HWY_ENABLE_CONTRIB OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(highway)

FetchContent_Declare(
  pocketfft
  GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
  GIT_TAG 5f27d5a8f51c5c25030cb22abf434decc9faf0ff
)
FetchContent_MakeAvailable(pocketfft)

set(DS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(DS_BUILD_ACCEPTANCE_PLUGIN OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  dualsynth2
  GIT_REPOSITORY https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2.git
  GIT_TAG master
)
FetchContent_MakeAvailable(dualsynth2)

find_package(Git REQUIRED)
execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --first-parent --tags --always
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_REPO_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REGEX REPLACE "^(r[0-9]+).*$" "\\1" VERSION "${GIT_REPO_VERSION}")
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-list --count HEAD
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_COMMIT_COUNT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/version.hpp.in"
  "${PROJECT_SOURCE_DIR}/src/version.hpp"
)
configure_file (
  "${PROJECT_SOURCE_DIR}/src/version.rc.in"
  "${PROJECT_SOURCE_DIR}/src/version.rc"
)

if(NOT MSVC)
  find_package(PkgConfig REQUIRED)

  pkg_check_modules(AVISYNTH avisynth)
  if(AVISYNTH_FOUND)
    include_directories(BEFORE include/avisynth)
    include_directories(${AVISYNTH_INCLUDE_DIRS})
  else()
    include_directories(include/avisynth)
  endif()

  pkg_check_modules(VAPOURSYNTH vapoursynth)
  if(VAPOURSYNTH_FOUND)
    include_directories(${VAPOURSYNTH_INCLUDE_DIRS})
  else()
    include_directories(include/vapoursynth)
  endif()
else()
  include_directories(include/avisynth)
  include_directories(include/vapoursynth)
endif()

include_directories(.)
include_directories(src)
include_directories(${highway_SOURCE_DIR})
include_directories(${pocketfft_SOURCE_DIR})

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if (CMAKE_GENERATOR_TOOLSET MATCHES "v[0-9]*_xp")
    target_compile_definitions(neo-dfttest PRIVATE WINVER=0x502 _WIN32_WINNT=0x502)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:threadSafeInit-")
  endif()
endif()

if (NOT WIN32)
  target_link_libraries(neo-dfttest dl)
endif()

target_link_libraries(neo-dfttest hwy)
target_link_libraries(neo-dfttest DualSynth::dualsynth)

add_custom_command(
  TARGET neo-dfttest POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:neo-dfttest> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:neo-dfttest>"
)
