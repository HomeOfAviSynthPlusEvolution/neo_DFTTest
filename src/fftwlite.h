// Lite version of fftw header on base of fftw3.h
// some needed fftwf typedefs added for delayed loading
// (by Fizick)
//
#ifndef __FFTWLITE_H__
#define __FFTWLITE_H__

#if _WIN32
  #ifndef NOMINMAX
  #define NOMINMAX
  #endif

  #include <windows.h>
  typedef HMODULE lib_t;
  typedef FARPROC func_t;
#else
  #include <dlfcn.h>
  typedef void* lib_t;
  typedef void* func_t;
  #define __stdcall
#endif

typedef float fftwf_complex[2];
typedef struct fftwf_plan_s  *fftwf_plan;
typedef void (*fftwf_destroy_plan_proc) (fftwf_plan);
typedef void (*fftwf_execute_dft_r2c_proc) (fftwf_plan, float *realdata, fftwf_complex *fftsrc);
typedef void (*fftwf_execute_dft_c2r_proc) (fftwf_plan, fftwf_complex *fftsrc, float *realdata);

using R = float;
using C = fftwf_complex;
typedef fftwf_plan (*fftwf_plan_dft_r2c_2d_proc)(int n0, int n1, R *in, C *out, unsigned flags);
typedef fftwf_plan (*fftwf_plan_dft_c2r_2d_proc)(int n0, int n1, C *in, R *out, unsigned flags);
typedef fftwf_plan (*fftwf_plan_dft_r2c_3d_proc)(int n0, int n1, int n2, R *in, C *out, unsigned flags);
typedef fftwf_plan (*fftwf_plan_dft_c2r_3d_proc)(int n0, int n1, int n2, C *in, R *out, unsigned flags);

#define FFTW_MEASURE (0U)
#define FFTW_DESTROY_INPUT (1U << 0)
#define FFTW_UNALIGNED (1U << 1)
#define FFTW_CONSERVE_MEMORY (1U << 2)
#define FFTW_EXHAUSTIVE (1U << 3) /* NO_EXHAUSTIVE is default */
#define FFTW_PRESERVE_INPUT (1U << 4) /* cancels FFTW_DESTROY_INPUT */
#define FFTW_PATIENT (1U << 5) /* IMPATIENT is default */
#define FFTW_ESTIMATE (1U << 6)
#define FFTW_WISDOM_ONLY (1U << 21)

#define LOAD_FFT_FUNC(name) do {name = reinterpret_cast<name ## _proc>((void*)fftw3_address(#name)); if (name == NULL) throw "Library function is missing: " #name; } while(0)

struct FFTFunctionPointers {
  lib_t library;

  fftwf_destroy_plan_proc fftwf_destroy_plan;
  fftwf_execute_dft_r2c_proc fftwf_execute_dft_r2c;
  fftwf_execute_dft_c2r_proc fftwf_execute_dft_c2r;
  fftwf_plan_dft_r2c_2d_proc fftwf_plan_dft_r2c_2d;
  fftwf_plan_dft_c2r_2d_proc fftwf_plan_dft_c2r_2d;
  fftwf_plan_dft_r2c_3d_proc fftwf_plan_dft_r2c_3d;
  fftwf_plan_dft_c2r_3d_proc fftwf_plan_dft_c2r_3d;

  #if _WIN32
    void fftw3_open() {
      library = LoadLibraryW(L"libfftw3f-3");
      if (library == NULL)
        library = LoadLibraryW(L"fftw3");
      if (library == NULL)
        #ifdef _WIN32
          throw("libfftw3f-3.dll or fftw3.dll not found. Please put in PATH or use LoadDll() plugin");
        #else
          throw("libfftw3f_threads.so.3 not found. Please install libfftw3-single3 (deb) or fftw-devel (rpm) package");
        #endif
    }
    void fftw3_close() { FreeLibrary(library); }
    func_t fftw3_address(LPCSTR func) { return GetProcAddress(library, func); }
  #else
    void fftw3_open() { library = dlopen("libfftw3f_threads.so.3", RTLD_NOW); }
    void fftw3_close() { dlclose(library); }
    func_t fftw3_address(const char * func) { return dlsym(library, func); }
  #endif
  void load() {
    library = NULL;
    fftw3_open();
    if (library != NULL) {
      LOAD_FFT_FUNC(fftwf_destroy_plan);
      LOAD_FFT_FUNC(fftwf_execute_dft_r2c);
      LOAD_FFT_FUNC(fftwf_execute_dft_c2r);
      LOAD_FFT_FUNC(fftwf_plan_dft_r2c_2d);
      LOAD_FFT_FUNC(fftwf_plan_dft_c2r_2d);
      LOAD_FFT_FUNC(fftwf_plan_dft_r2c_3d);
      LOAD_FFT_FUNC(fftwf_plan_dft_c2r_3d);
    }
  }

  void free() {
    if (library != NULL) {
      fftw3_close();
    }
  }
};

#undef LOAD_FFT_FUNC

#endif
