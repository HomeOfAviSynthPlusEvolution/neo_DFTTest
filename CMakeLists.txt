cmake_minimum_required (VERSION 3.24)
project (neo_dfttest CXX)
add_library(neo-dfttest SHARED
  src/plugin/plugin_entry.cpp
  src/plugin/dfttest_filter.cpp
  src/plugin/dfttest_filter.hpp
  src/version.rc
  src/core.cpp
  src/core.h
  src/core_C.cpp
  src/core_hwy.cpp
  src/dft_common.h
  src/engine/dfttest_engine.cpp
  src/engine/dfttest_engine.hpp
  src/fftwlite.h
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

option(ENABLE_PAR "Enable C++17 Parallel Execution" OFF)
target_compile_definitions(neo-dfttest PRIVATE _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)

include(FetchContent)
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

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if (CMAKE_GENERATOR_TOOLSET MATCHES "v[0-9]*_xp")
    target_compile_definitions(neo-dfttest PRIVATE WINVER=0x502 _WIN32_WINNT=0x502)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:threadSafeInit-")
  endif()
endif()

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(execution HAS_EXECUTION)
if(HAS_EXECUTION)
  add_definitions(-DHAS_EXECUTION)
endif()
if(ENABLE_PAR AND HAS_EXECUTION)
  add_definitions(-DENABLE_PAR)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(neo-dfttest tbb)
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
