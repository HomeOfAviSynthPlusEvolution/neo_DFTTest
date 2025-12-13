cmake_minimum_required (VERSION 3.14)
project (dfttest_tests CXX)

# Find and configure Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip # v1.14.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(core_internal_tests
  test_main.cpp
  test_core_internal.cpp
)

target_link_libraries(core_internal_tests
  GTest::gtest_main
)

# Relative to tests directory
target_include_directories(core_internal_tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/.. # for root level files like CMakeLists.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/avisynth
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/vapoursynth
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/dualsynth
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

set_property(TARGET core_internal_tests PROPERTY CXX_STANDARD 17)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_source_files_properties(test_core_internal.cpp PROPERTIES COMPILE_FLAGS "/arch:AVX2")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(test_core_internal.cpp PROPERTIES COMPILE_FLAGS "/arch:AVX2")
else()
  set_source_files_properties(test_core_internal.cpp PROPERTIES COMPILE_FLAGS "-mfma -mavx2")
endif()

add_test(NAME CoreInternalTests COMMAND core_internal_tests)

# Usage:
# mkdir -p build/msvc-x64
# cmake -G "Visual Studio 18 2026" -A x64 -B build/msvc-x64 -S . -DBUILD_GMOCK=OFF
# cmake --build build/msvc-x64 --config Release
# mkdir -p build/gcc-x64
# cmake -G "MSYS Makefiles" -B build/gcc-x64 -S . -DBUILD_GMOCK=OFF
# cmake --build build/gcc-x64 --config Release
