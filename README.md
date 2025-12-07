# Neo DFTTest (forked from VapourSynth-DFTTest)

Neo DFTTest Copyright(C) 2020 Xinyue Lu, and previous developers

DFTTest is a 2D/3D Frequency Domain denoiser. It was originally written by tritical, and later modified by Ferenc Pintér aka pinterf for further improvement, high bit depth, and more. VapourSynth-DFTTest was ported to VapourSynth interface by HolyWu with further clean ups and efficient SSE2 and AVX2 SIMD routines. Kudos to them for creating and improving this fantastic tool.

This project backports VapourSynth-DFTTest to AviSynth+, with minor improvement on single precision floating point reciprocal using Newton Raphson method to produce bit identical result across C and SIMD routines.

Parameter names follow VapourSynth-DFTTest, some of which are different than the original AVS counterpart.

Requires libfftw3f-3.dll to be in the search path. Download from: http://www.fftw.org/install/windows.html .

## Usage

```python
# AviSynth+
LoadPlugin("neo-dfttest.dll")
neo_dfttest(clip, ftype=0, sigma=2.0, y=3, u=3, v=3, ...)
# VapourSynth
core.neo_dfttest.DFTTest(clip, ftype=0, sigma=2.0, planes=[0,1,2], ...)
```

Parameters:

[Check original dfttest usage documents.](https://github.com/pinterf/dfttest/blob/master/dfttest%20-%20README.txt)

[Check original VapourSynth-DFTTest usage documents.](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-DFTTest/blob/master/README.md)

- *sigma*, *sigma2*

    Default: 8.0.

- *sosize*

    Default: 12.

- *nlocation*, *alpha*

    Replace nfile and nstring in dfttest.

    In AviSynth+, *nlocation* is a string of group of integers, where delimiter can be a space, a comma, or a colon.

    In VapourSynth, *nlocation* is an array of integers.

- *slocation*, *ssx*, *ssy*, *sst*, *ssystem*

    Replace sstring by *slocation* and *ssystem*.

    In AviSynth+, *slocation*, *ssx*, *ssy*, *sst* are strings of group of floats, where delimiter can be a space, a comma, or a colon.

    '$' sign in sstring is now *ssystem = 1*.

    In VapourSynth, *slocation*, *ssx*, *ssy*, *sst* are arrays of floats.

- *dither*

    Only works in 8-bit and is the same as dfttest.

    Default: 0.

- *y*, *u*, *v*, *a* (AviSynth+ only)

    Whether a plane is to be filtered.

        1 - Do not touch, leaving garbage data
        2 - Copy from origin
        3 - Process

    Default: 3.

- *planes* (VapourSynth only)

    Planes to be filtered.

    Default: [0,1,2].

- *threads*

    Number of internal multi threads. When used in frame level multi threading, please adjust external threads and internal *threads* down to reduce overhead.

    For example, on a 6 cores 12 threads CPU, *threads=6* paired with *prefetch(2)* comes with great efficiency; *threads=6* paired with *prefetch(6)* is likely less efficient.

    Default: 4 for AviSynth+, 1 for VapourSynth. Max: 16.

- *fft_threads*

    Number of FFTW multi threads.

    Default: 2 for AviSynth+, 1 for VapourSynth.

- *opt*

    Sets which CPU optimizations to use.

        0 - Auto detect
        1 - Use C
        2 - Use SSE2 if supported
        3 - Use AVX2 if supported

## License

* GPLv2.
