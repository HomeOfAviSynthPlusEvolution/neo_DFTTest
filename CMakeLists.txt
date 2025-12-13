cmake_minimum_required (VERSION 3.5)
project (neo_dfttest CXX)
add_library(neo-dfttest SHARED
  main.cpp
  src/version.rc
  src/vectorclass/instrset_detect.cpp
  src/core.cpp
  src/core.h
  src/core_AVX2.cpp
  src/core_C.cpp
  src/core_SSE2.cpp
  src/core_hwy.cpp
  src/dft_common.h
  src/fftwlite.h
  src/neo_dfttest.hpp
)
set_property(TARGET neo-dfttest PROPERTY CXX_STANDARD 17)
option(ENABLE_PAR "Enable C++17 Parallel Execution" ON)
add_compile_definitions(VS_TARGET_CPU_X86 _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)

include(FetchContent)
FetchContent_Declare(
  highway
  GIT_REPOSITORY https://github.com/google/highway.git
  GIT_TAG 1.3.0
)
FetchContent_MakeAvailable(highway)

find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} describe --first-parent --tags --always OUTPUT_VARIABLE GIT_REPO_VERSION)
string(REGEX REPLACE "(r[0-9]+).*" "\\1" VERSION ${GIT_REPO_VERSION})
execute_process(COMMAND ${GIT_EXECUTABLE} rev-list --count HEAD OUTPUT_VARIABLE GIT_COMMIT_COUNT OUTPUT_STRIP_TRAILING_WHITESPACE)

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
include_directories(include/dualsynth)
include_directories(${highway_SOURCE_DIR})

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_source_files_properties(src/core_SSE2.cpp PROPERTIES COMPILE_FLAGS "/arch:SSE2")
  set_source_files_properties(src/core_AVX2.cpp PROPERTIES COMPILE_FLAGS "/arch:AVX2")

  if (CMAKE_GENERATOR_TOOLSET MATCHES "v[0-9]*_xp")
    target_compile_definitions(neo-dfttest PRIVATE WINVER=0x502 _WIN32_WINNT=0x502)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:threadSafeInit-")
  endif()

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(src/core_SSE2.cpp PROPERTIES COMPILE_FLAGS "/arch:SSE2")
  set_source_files_properties(src/core_AVX2.cpp PROPERTIES COMPILE_FLAGS "/arch:AVX2")
  target_link_libraries(neo-dfttest libmmds)

else()
  set_source_files_properties(src/core_SSE2.cpp PROPERTIES COMPILE_FLAGS "-msse2")
  set_source_files_properties(src/core_AVX2.cpp PROPERTIES COMPILE_FLAGS "-mfma -mavx2")

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

target_link_libraries(neo-dfttest PRIVATE hwy)

add_custom_command(
  TARGET neo-dfttest POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:neo-dfttest> "../Release_${VERSION}/${_DIR}/$<TARGET_FILE_NAME:neo-dfttest>"
)
