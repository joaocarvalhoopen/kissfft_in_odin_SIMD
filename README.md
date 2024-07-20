# kissfft_in_odin_SIMD
KISSFFT SIMD very fast bindings for the Odin programming language. 

## Description
This are simple but very fast bindings in Odin to the f32 FFT and IFFT of KISS_FFT with SIMD SSE optimization for complex input and output values. The trick that makes this bindings so fast is the use of SIMD and that they are not compiled with GCC, like the KISS_FFT project, but with CLANG, with fast optimization flags. This will do some auto vectorization. The bindings are for the f32 type, __mm128 bits data type. The KISS_FFT library is available at it's github page and has a good license. Note that the allocation and dealocation of buffers are made in Odin, but the FFT and IFFT and the creation of the plan are made in C. I made 3 examples of it's usage, and it's very simple to use. In principle this bindings can be used in multithreaded applications. There is also the possibility to speed up the FFT and IFFT by using the CLANG optimization flags ```-ffast-math```. See Makefile. But you should do your own tests to see if it's worth it for your application. The C files present in this lib are from KISSFFT.

## Original GitHub of KISS_FFT 
[https://github.com/mborgerding/kissfft](https://github.com/mborgerding/kissfft)

## Speed up to the normal kiss_fft_in_odin
3.51x times faster. 

## To compile this lib do

``` bash
$ make kiss_fft_simd
$ make
$ time ./kiss_fft_simd.exe

or

$ make kiss_fft_simd
$ make opti_max
$ time ./kiss_fft_simd.exe
```

## License
The same of KISS_FFT library. BSD-3-Clause

## Have fun!
Best regards, <br>
Joao Carvalho
