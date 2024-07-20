all:
	odin build . -out:kiss_fft_simd.exe --debug

opti:
	odin build . -out:kiss_fft_simd.exe -o:speed

opti_max:	
	odin build . -out:kiss_fft_simd.exe -o:aggressive -microarch:native -no-bounds-check -disable-assert

# Target to build the kissfft optimized static library.
kiss_fft_simd: ./kiss_fft_simd_odin/lib_my_kiss_fft_simd.a

# Rule to create the static library.
./kiss_fft_simd_odin/lib_my_kiss_fft_simd.a: ./kiss_fft_simd_odin/my_kiss_fft_simd.o
	ar rcs ./kiss_fft_simd_odin/lib_my_kiss_fft_simd.a ./kiss_fft_simd_odin/my_kiss_fft_simd.o

# Rule to compile the object file.
# Note : The -ffast-math flag is used to enable fast math optimizations,
#        but I didn t test if the results are correct.
./kiss_fft_simd_odin/my_kiss_fft_simd.o: ./kiss_fft_simd_odin/kiss_fft.c
	clang -c ./kiss_fft_simd_odin/kiss_fft.c -o ./kiss_fft_simd_odin/my_kiss_fft_simd.o -O3 -march=native -funroll-loops -msse
# 	clang -c ./kiss_fft_simd_odin/kiss_fft.c -o ./kiss_fft_simd_odin/my_kiss_fft_simd.o -O3 -march=native -ffast-math -funroll-loops

clean:
	rm kiss_fft_simd.exe
	rm ./kiss_fft_simd_odin/my_kiss_fft_simd.o
	rm ./kiss_fft_simd_odin/lib_my_kiss_fft_simd.a

run:
	./kiss_fft_simd.exe



