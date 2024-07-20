// Bindings name : kiss_fft_simd_in_odin
//
// Description : This are simple but very fast bindings in Odin to the f32 FFT and IFFT
//               of KISS_FFT with SIMD.
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


package kiss_fft_simd_odin

import "base:runtime"

import "core:fmt"
import "core:strings"
import "core:c"
import "core:c/libc"
import "core:math"
import "core:math/cmplx"
import "core:mem"
import "core:os"

Dir_FFT  :: 0
Dir_IFFT :: 1

// External C object import.
when ODIN_OS == .Linux do foreign import foo { "./lib_my_kiss_fft_simd.a" }

foreign foo {

    // typical usage:      kiss_fft_cfg mycfg = kiss_fft_alloc( 1024, 0, NULL, NULL );
    kiss_fft_alloc :: proc "c" ( nfft : c.int, inverse_fft : c.int, mem : rawptr, lenmem : rawptr ) -> rawptr ---


    // usage: kiss_fft( mycfg, fin, fout );
    //
    // void KISS_FFT_API kiss_fft( kiss_fft_cfg cfg, const kiss_fft_cpx *fin, kiss_fft_cpx *fout );
    kiss_fft :: proc "c" ( cfg : rawptr, buf_in : [ ^ ]complex64, buf_out : [ ^ ]complex64 ) ---

    // void SSETools_pack_2_128_complex( float * target,
    //                                   float * source_fft_0,
    //                                   float * source_fft_1,
    //                                   float * source_fft_2,
    //                                   float * source_fft_3,
    //                                   unsigned long num_elem_fft )

    // num_elem_fft <--- Must be a multiple of 4.
    //                   TODO: Do better and partial fill the last SIMD vector.
    SSETools_pack_2_128_complex :: proc "c" ( target       : [ ^ ]f32,
                                              source_fft_0 : [ ^ ]f32,
                                              source_fft_1 : [ ^ ]f32,
                                              source_fft_2 : [ ^ ]f32,
                                              source_fft_3 : [ ^ ]f32,
                                              num_elem_fft : c.uint64_t ) ---

    // void SSETools_unpack_2_128_complex( float * target_fft_0,
    //                                     float * target_fft_1,
    //                                     float * target_fft_2,
    //                                     float * target_fft_3,
    //                                     float * source,
    //                                     unsigned long num_elem_fft )

    // num_elem_fft <--- Must be a multiple of 4.
    //                   TODO: Do better and partial fill the last SIMD vector.
    SSETools_unpack_2_128_complex :: proc "c" ( target_fft_0 : [ ^ ]f32,
                                                target_fft_1 : [ ^ ]f32,
                                                target_fft_2 : [ ^ ]f32,
                                                target_fft_3 : [ ^ ]f32,
                                                source       : [ ^ ]f32,
                                                num_elem_fft : c.uint64_t ) ---

}

apply_ifft_correction_factor :: proc ( vec     : [ ^ ]complex64,
                                       n_elems : c.int           ) {

    for i in 0 ..< n_elems {
        vec[ i ] = vec[ i ] / f32_to_c64( f32( n_elems ), 0.0 )
    }
}

f32_to_c64 :: #force_inline proc ( re, im : f32 ) -> complex64 {

    return complex64( complex( re, im ) )
}

print_vec :: proc ( var_name : string,
                    vec      : [ ^ ]complex64,
                    n_elems  : c.int           ) {

    for i in 0 ..< n_elems {
        fmt.printf( "%s[ %d ] = %f + %fi\n",
                    var_name,
                    i,
                    real( vec[ i ] ),
                    imag( vec[ i ] )  )
    }
}

print_vec_simd :: proc ( var_name : string,
                         vec      : [ ^ ]f32,
                         n_elems  : c.int     ) {

    for i in 0 ..< n_elems * 2 {
        fmt.printf( "%s_f32[ %d ] = %f\n",
                    var_name,
                    i,
                    vec[ i ] )
    }
}

my_make_aligned :: proc ( $T        : typeid/[]$E,
                          len       : int,
                          alignment : int,
                          allocator := context.allocator,
                          loc       := #caller_location   ) -> 
                        ( T, mem.Allocator_Error ) {
    runtime.make_slice_error_loc( loc, len )
    data, err := runtime.mem_alloc_bytes( size_of( E ) * len, alignment, allocator, loc )
    if data == nil && size_of( E ) != 0 {
        return nil, err
    }
    s := runtime.Raw_Slice{ raw_data( data ), len }
    return transmute( T )s, err
}

//
//  Test's
//

test_kiss_fft_SIMD_FFT_IFFT :: proc ( ) {

    fmt.printfln( "===>>> test SIMD FFT followed by IFFT.\n" )


    n_elems_fft : c.int = 16_000     // 1024    // 16

    len_rows :: 8 * 1_000    // 8 * 20_000     // 8 

    // Allocate the 2D Arrays buf_original and buf_destination
    rows_fft_buf_original : [ ][ ^ ]complex64 = make( [ ][ ^ ]complex64, len_rows )
    rows_fft_buf_out_dest : [ ][ ^ ]complex64 = make( [ ][ ^ ]complex64, len_rows )
    defer delete( rows_fft_buf_original )
    defer delete( rows_fft_buf_out_dest )
    
    // Multi-pointer slices used inside the SIMD functions,
    // this have to be aligned do 128 bits, 16 bytes, 4x f32.
    rows_fft_buf_in_simd  : [ ^ ]complex64
    rows_fft_buf_out_simd : [ ^ ]complex64
   
    //
    // FFT: Fast Fourier Transform plan creation.
    //

    cfg_fft := kiss_fft_alloc( n_elems_fft,
                               Dir_FFT,
                               nil,
                               nil )
    
    // Free the plan configuration.
    defer libc.free( cfg_fft )

    assert( cfg_fft != nil, "Failed to allocate kiss_fft_cfg_fft" )
    

    // Allocate input and output buffers

    // Allocate buf_original.
    for row in 0 ..< len_rows {
        rows_fft_buf_original[ row ] = make( [ ^ ]complex64, n_elems_fft )
    }

    defer { 
            for row in 0 ..< len_rows {
                free( rows_fft_buf_original[ row ] )
            }    
        }


    // Allocate buf_in_simd with aligned memory to 128 bits.
    n_elems_fft_simd := n_elems_fft * 4 

    aligment_128_in_bytes : int = 128 / 8     // 128 bits = 16 bytes 

    buf_in_tmp : [ ]complex64
    buf_in     : [ ^ ]complex64
    err_1      : runtime.Allocator_Error
    buf_in_tmp, err_1 = my_make_aligned( [  ]complex64,
                                         int( n_elems_fft_simd ),
                                         aligment_128_in_bytes )
    buf_in = raw_data( buf_in_tmp )

    rows_fft_buf_in_simd = buf_in

    defer free( rows_fft_buf_in_simd )



    // Allocate buf_out_simd with aligned memory to 128 bits.
    buf_out_tmp : [ ]complex64
    buf_out     : [ ^ ]complex64
    err_2       : runtime.Allocator_Error
    buf_out_tmp, err_2 = my_make_aligned( [  ]complex64,
                                          int( n_elems_fft_simd ),
                                          aligment_128_in_bytes )
    buf_out = raw_data( buf_out_tmp )

    rows_fft_buf_out_simd = buf_out

    defer free( rows_fft_buf_out_simd )


    // Allocate buf_out_dest.
    for row in 0 ..< len_rows {
        rows_fft_buf_out_dest[ row ] = make( [ ^ ]complex64, n_elems_fft )
    }
    defer { 
            for row in 0 ..< len_rows {
                free( rows_fft_buf_out_dest[ row ] )
            }
        }


    // Fill in buf_original and buf_in with some data
    for row in 0 ..< len_rows {
        for i in 0 ..< n_elems_fft {
            // re := math.cos( 2.0 * math.PI * ( f32( i ) / f32(  n_elems_fft / 3 ) ) )
            // rows_fft_buf_original[ row ][ i ] = f32_to_c64( re, re )

            re := f32( i )
            cm := f32( i ) + 0.5
            rows_fft_buf_original[ row ][ i ] = f32_to_c64( re, cm )
        }
    }

    /*
    // Print the rows of input buffers.
    for row in 0 ..< len_rows {
        // Print the input buffer
        row_str := fmt.aprintf( "buf_original _ rows[ %v ] ", row )
        print_vec( row_str, rows_fft_buf_original[ row ], n_elems_fft )
    }
    */

    assert( len_rows % 4 == 0, "len_rows must be a multiple of 4" )

    for row_div in 0 ..< len_rows / 4 {

        row := row_div * 4

        // Pack 4 rows of the original buffer into a SIMD buffer.
        target_simd := transmute( [ ^ ]f32 ) rows_fft_buf_in_simd

        source_fft_0 := transmute( [ ^ ]f32 ) rows_fft_buf_original[ row ]
        source_fft_1 := transmute( [ ^ ]f32 ) rows_fft_buf_original[ row + 1 ]
        source_fft_2 := transmute( [ ^ ]f32 ) rows_fft_buf_original[ row + 2 ]
        source_fft_3 := transmute( [ ^ ]f32 ) rows_fft_buf_original[ row + 3 ]

        num_elem_fft_u64 := c.uint64_t( n_elems_fft )

        SSETools_pack_2_128_complex( target_simd,
                                     source_fft_0,
                                     source_fft_1,
                                     source_fft_2,
                                     source_fft_3,
                                     num_elem_fft_u64 )


        /*
        // print the SIMD buffer.
        if row == 0 {
            fmt.printfln( "SIMD buffer rows[ %d ]", row )
            buf_in_simd_ptr := transmute( [ ^ ]f32 ) rows_fft_buf_in_simd
            print_vec_simd( "buf_in_simd", buf_in_simd_ptr, n_elems_fft * 4 )
        }
        */


        // Perform the FFT on a SIMD buffer x4 f32 ( x2 for real part and complex part).
        kiss_fft( cfg_fft, rows_fft_buf_in_simd, rows_fft_buf_out_simd );

        // Unpack the SIMD buffer into 4 rows of the destination buffer.
        target_fft_0 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row ]
        target_fft_1 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 1 ]
        target_fft_2 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 2 ]
        target_fft_3 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 3 ]

        source_simd := transmute( [ ^ ]f32 ) rows_fft_buf_out_simd

        SSETools_unpack_2_128_complex( target_fft_0,
                                       target_fft_1,
                                       target_fft_2,
                                       target_fft_3,
                                       source_simd,
                                       num_elem_fft_u64 )

    }

    /*
    // Print the rows of input buffers.
    for row in 0 ..< len_rows {
        // Print the input buffer
        row_str := fmt.aprintf( "buf_out_dest _ rows[ %v ] ", row )
        print_vec( row_str, rows_fft_buf_out_dest[ row ], n_elems_fft )
    }
    */
    

    //
    // IFFT: Inverse Fast Fourier Transform plan creation.
    //
    
    cfg_ifft := kiss_fft_alloc( n_elems_fft,
                                Dir_IFFT,
                                nil,
                                nil )
    
    // Free the plan configuration.
    defer libc.free( cfg_ifft )

    assert( cfg_ifft != nil, "Failed to allocate kiss_fft_cfg_ifft" )


    // Perform the IFFT
    for row in 0 ..< len_rows / 4 {

        row := row * 4

        // Pack 4 rows of the original buffer into a SIMD buffer.
        target_simd := transmute( [ ^ ]f32 ) rows_fft_buf_in_simd

        source_fft_0 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row ]
        source_fft_1 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 1 ]
        source_fft_2 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 2 ]
        source_fft_3 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 3 ]

        num_elem_fft_u64 := c.uint64_t( n_elems_fft )

        SSETools_pack_2_128_complex( target_simd,
                                     source_fft_0,
                                     source_fft_1,
                                     source_fft_2,
                                     source_fft_3,
                                     num_elem_fft_u64 )

        // Perform the IFFT on a SIMD buffer x4 f32 ( x2 for real part and complex part).
        kiss_fft( cfg_ifft, rows_fft_buf_in_simd, rows_fft_buf_out_simd );

        // Unpack the SIMD buffer into 4 rows of the destination buffer.
        target_fft_0 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row ]
        target_fft_1 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 1 ]
        target_fft_2 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 2 ]
        target_fft_3 := transmute( [ ^ ]f32 ) rows_fft_buf_out_dest[ row + 3 ]

        source_simd := transmute( [ ^ ]f32 ) rows_fft_buf_out_simd

        SSETools_unpack_2_128_complex( target_fft_0,
                                       target_fft_1,
                                       target_fft_2,
                                       target_fft_3,
                                       source_simd,
                                       num_elem_fft_u64 )
    }


    /*
    // Print the rows of input buffers.
    for row in 0 ..< len_rows {
        // Print the input buffer
        row_str := fmt.aprintf( "buf_out_dest _ rows[ %v ] ", row )
        print_vec( row_str, rows_fft_buf_out_dest[ row ], n_elems_fft )
    }
    */

    // Apply the correction factor to the FFT,
    // because the FFT made by KISS_FFT doesn't apply automatically the
    // correction factor.
    for row in 0 ..< len_rows {
        apply_ifft_correction_factor( rows_fft_buf_out_dest[ row ], n_elems_fft )
    }

    // Check if the original and the IFFT out_dest are the same.
    for row in 0 ..< len_rows {
        for i in 0 ..< n_elems_fft {
           /*
            fmt.printfln( "buf_original[ %d ] = %f + %fi  buf_out_dest[ %d ] = %f + %fi",
                          i,
                          real( rows_fft_buf_original[ row ][ i ] ),
                          imag( rows_fft_buf_original[ row ][ i ] ),
                          i,
                          real( rows_fft_buf_out_dest[ row ][ i ] ),
                          imag( rows_fft_buf_out_dest[ row ][ i ] ) )
           */

            // assert( cmplx.abs( rows_fft_buf_original[ row ][ i ] - rows_fft_buf_out_dest[ row ][ i ] ) < 0.0001, "IFFT failed" )
            val_abs := cmplx.abs( rows_fft_buf_original[ row ][ i ] - rows_fft_buf_out_dest[ row ][ i ] ) 
            if val_abs >= 0.01 { //  0.0001 
                fmt.printfln( "FTT of IFFT failed! val_abs : %v", val_abs )
                os.exit( 1 )
            }

        }
    }

    fmt.printfln( "...end test SIMD FFT followed by IFFT PASSED.\n" )
}

//
// Previous tests for KISS_FFT without SIMD.
// Not applyable to KISS_FFT SIMD.
// 

/*

test_kiss_fft :: proc ( ) {

    fmt.printfln( "===>>> test FFT.\n" )

    n_elems_fft : c.int = 1024
    cfg_fft := kiss_fft_alloc( n_elems_fft,
                               Dir_FFT,
                               nil,
                               nil )
    // free the configuration
    defer libc.free( cfg_fft )

    assert( cfg_fft != nil, "Failed to allocate kiss_fft_cfg_fft" )
    
    // Allocate input and output buffers
    buf_in  : [ ^ ]complex64 = make( [ ^ ]complex64, n_elems_fft )
    defer free( buf_in )

    buf_out : [ ^ ]complex64 = make( [ ^ ]complex64, n_elems_fft )
    defer free( buf_out )

    // Fill in buf_in with some data
    for i in 0 ..< n_elems_fft {
        re := math.cos( 2.0 * math.PI * ( f32( i ) / f32(  n_elems_fft / 10 ) ) )
        buf_in[ i ] = f32_to_c64( re, 0.0 )
    }

    // Print the input buffer
    print_vec( "buf_in", buf_in, n_elems_fft)

    // Perform the FFT
    kiss_fft( cfg_fft, buf_in, buf_out );

    // Print the buf_out ( output buffer ).
    print_vec( "buf_out", buf_out, n_elems_fft )

    fmt.printfln( "...end test FFT.\n" )
}

*/

/*

test_kiss_fft_reuse_a_vector :: proc ( ) {

    fmt.printfln( "===>>> test FFT reuse a vector.\n" )

    n_elems_fft : c.int = 16
    cfg_fft := kiss_fft_alloc( n_elems_fft,
                               Dir_FFT,
                               nil,
                               nil )
    // free the configuration
    defer libc.free( cfg_fft )

    assert( cfg_fft != nil, "Failed to allocate kiss_fft_cfg_fft" )
    
    // Allocate input and output buffers
    vec_in  : [ ]complex64 = make( [ ]complex64, n_elems_fft )
    defer delete( vec_in )

    vec_out : [  ]complex64 = make( [ ]complex64, n_elems_fft )
    defer delete( vec_out )

    // Get the raw data from the slices, and make the convertion of types.
    buf_in  : [ ^ ]complex64 = raw_data( vec_in )

    buf_out : [ ^ ]complex64 = raw_data( vec_out )

    // Fill in buf_in with some data
    for i in 0 ..< n_elems_fft {
        re := math.cos( 2.0 * math.PI * ( f32( i ) / f32(  n_elems_fft / 10 ) ) )
        vec_in[ i ] = f32_to_c64( re, 0.0 )
    }

    // Print the input buffer
    print_vec( "buf_in", buf_in, n_elems_fft)

    // Perform the FFT
    kiss_fft( cfg_fft, buf_in, buf_out );

    // Print the buf_out ( output buffer ).
    print_vec( "buf_out", buf_out, n_elems_fft )

    fmt.printfln( "...end test FFT reuse a vector.\n" )
}

*/

/*

test_kiss_fft_FFT_IFFT :: proc () {

    fmt.printfln( "===>>> test FFT followed by IFFT.\n" )

    n_elems_fft : c.int = 16 // 1024

    //
    // FFT: Fast Fourier Transform
    //

    cfg_fft := kiss_fft_alloc( n_elems_fft,
                               Dir_FFT,
                               nil,
                               nil )
    // Free the configuration.
    defer libc.free( cfg_fft )

    assert( cfg_fft != nil, "Failed to allocate kiss_fft_cfg_fft" )
    
    // Allocate input and output buffers
    buf_original  : [ ]complex64 = make( [ ]complex64, n_elems_fft )
    defer delete( buf_original )

    buf_in  : [ ^ ]complex64 = make( [ ^ ]complex64, n_elems_fft )
    defer free( buf_in )
    
    buf_out : [ ^ ]complex64 = make( [ ^ ]complex64, n_elems_fft )
    defer free( buf_out )

    // Fill in buf_in with some data
    for i in 0 ..< n_elems_fft {
        re := math.cos( 2.0 * math.PI * ( f32( i ) / f32(  n_elems_fft / 3 ) ) )
        // buf_in[ i ]       = f32_to_c64( re, 0.0 )
        // buf_original[ i ] = f32_to_c64( re, 0.0 )
        buf_in[ i ]       = f32_to_c64( re, re )
        buf_original[ i ] = f32_to_c64( re, re )
    }

    // Print the input buffer
    print_vec( "buf_in", buf_in, n_elems_fft )

    // Perform the FFT
    kiss_fft( cfg_fft, buf_in, buf_out );

    // Print the buf_out ( output buffer ).
    print_vec( "buf_out", buf_out, n_elems_fft )


    //
    // IFFT: Inverse Fast Fourier Transform
    //
    
    cfg_ifft := kiss_fft_alloc( n_elems_fft,
                                Dir_IFFT,
                                nil,
                                nil )
    // Free the configuration.
    defer libc.free( cfg_ifft )

    assert( cfg_ifft != nil, "Failed to allocate kiss_fft_cfg_ifft" )

    // Perform the IFFT
    kiss_fft( cfg_ifft, buf_out, buf_in );

    // Print the buf_in ( output buffer ).
    print_vec( "buf_in_that_is_ouput", buf_in, n_elems_fft )

    // Apply the correction factor.
    apply_ifft_correction_factor( buf_in, n_elems_fft )

    // Check if the original and the IFFT are the same.
    for i in 0 ..< n_elems_fft {
        fmt.printfln( "buf_original[ %d ] = %f + %fi  buf_in[ %d ] = %f + %fi\n",
                      i,
                      real( buf_original[ i ] ),
                      imag( buf_original[ i ] ),
                      i,
                      real( buf_in[ i ] ),
                      imag( buf_in[ i ] ) )

        assert( cmplx.abs( buf_original[ i ] - buf_in[ i ] ) < 0.0001, "IFFT failed" )
    }

    fmt.printfln( "...end test FFT followed by IFFT PASSED.\n" )
}

*/

