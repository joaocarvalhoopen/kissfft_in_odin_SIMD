// Bindings name : kissfft_in_odin_SIMD
//
// Description : This are simple but very fast bindings in Odin to the f32 FFT and IFFT
//               of KISS_FFT with SIMD SSE optimization for complex input and output values.
//               The KISS_FFT is a small, fast, and simple FFT library written in C.
//               The trick that makes this bindings so fast is that they are not compiled
//               with GCC, like the KISS_FFT project, but with CLANG, with fast optimization
//               flags. This will do some auto vectorization.
//               The bindings are for the f32 type.
//               The KISS_FFT library is available at it's github page and has a good license.
//               Note that the allocation and dealocation of buffers are made in Odin, but the
//               FFT and IFFT and the creation of the plan are made in C.
//               I made 3 examples of it's usage, and it's very simple to use.
//               In principle this bindings can be used in multithreaded applications.
//               There is also the possibility to speed up the FFT and IFFT by using the CLANG
//               optimization flags ```-ffast-math````. See make file.
//               But you should do your own tests to see if it's worth it for your application.
//               The C files present in this lib are from KISSFFT.
//
//
// Author of the bindings : Joao Carvalho
// Date                   : 2024-07-20
//
// Original GitHub of KISS_FFT : 
//        https://github.com/mborgerding/kissfft
//
// To compile this lib do:
//
//       $ make kiss_fft_simd
//       $ make
//       $ time ./kiss_fft_simd.exe
//
//       or
//
//       $ make kiss_fft_simd
//       $ make opti_max
//       $ time ./kiss_fft_simd.exe
//
//
// License: The same of KISS_FFT library.
//          BSD-3-Clause
//
// Have fun.

package main

import "core:fmt"
import "core:strings"

import fft "./kiss_fft_simd_odin"

main:: proc () {
    fmt.printfln("Kiss FFT SIMD SSSE, run tests...\n")

    fft.test_kiss_fft_SIMD_FFT_IFFT( )

//    fft.test_kiss_fft( )

//    fft.test_kiss_fft_reuse_a_vector( )

//    fft.test_kiss_fft_FFT_IFFT( )

    fmt.printfln("..end Kiss FFT SIMD SSSE, tests done.\n")
}