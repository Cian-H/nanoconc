using Test
using PropCheck

include("../anchors.jl")

import .Anchors: TEST_DIR, SRC_DIR

if !@isdefined TestUtils
    include(joinpath(TEST_DIR, "testutils.jl"))
end
if !@isdefined miemfp
    include(joinpath(SRC_DIR, "miemfp.jl"))
end
if !@isdefined FFIWraps
    include(joinpath(TEST_DIR, "ffi_wraps.jl"))
end

# function julia_vs_c(args::Tuple{Float64, Float64, Float64, UInt32, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}})
#     x, cxref_re, cxref_im, nang, cxs1_re, cxs1_im, cxs2_re, cxs2_im = args
function julia_vs_c(x, cxref_re, cxref_im, nang, cxs1_re, cxs1_im, cxs2_re, cxs2_im)
    cxref, cxs1, cxs2 = ComplexF64(cxref_re, cxref_im), ComplexF64.(cxs1_re, cxs1_im), ComplexF64.(cxs2_re, cxs2_im)
    x_c, cxref_c, nang_c, cxs1_c, cxs2_c = Float32(x), ComplexF32(cxref), UInt32(nang), ComplexF32.(cxs1), ComplexF32.(cxs2)
    return isapprox(
        miemfp.bhmie(x, cxref, nang),
        FFIWraps.bhmie_c(x_c, cxref_c, nang_c, cxs1_c, cxs2_c),
        rtol=0.1,
    )
end

# function julia_vs_fortran(args::Tuple{Float64, Float64, Float64, UInt32, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}})
#     x, cxref_re, cxref_im, nang, cxs1_re, cxs1_im, cxs2_re, cxs2_im = args
function julia_vs_fortran(x, cxref_re, cxref_im, nang, cxs1_re, cxs1_im, cxs2_re, cxs2_im)
    cxref, cxs1, cxs2 = ComplexF64(cxref_re, cxref_im), ComplexF64.(cxs1_re, cxs1_im), ComplexF64.(cxs2_re, cxs2_im)
    x_f, cxref_f, nang_f, cxs1_f, cxs2_f = Float32(x), ComplexF32(cxref), Int32(nang), ComplexF32.(cxs1), ComplexF32.(cxs2)
    b = miemfp.bhmie(x, cxref, nang)
    f = FFIWraps.bhmie_fortran(x_f, cxref_f, nang_f, cxs1_f, cxs2_f)
    # open("bhmie_julia_vs_fortran.txt", "a") do io
    #     println(io, "julia: ", b)
    #     println(io, "fortran: ", f)
    # end
    # return is_approx(b, f)
    return isapprox(
        miemfp.bhmie(x, cxref, nang),
        FFIWraps.bhmie_fortran(x_f, cxref_f, nang_f, cxs1_f, cxs2_f),
        rtol=0.1,
    )
end

# function julia_vs_fortran77(args::Tuple{Float64, Float64, Float64, UInt32, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}})
#     x, cxref_re, cxref_im, nang, cxs1_re, cxs1_im, cxs2_re, cxs2_im = args
function julia_vs_fortran77(x, cxref_re, cxref_im, nang, cxs1_re, cxs1_im, cxs2_re, cxs2_im)
    cxref, cxs1, cxs2 = ComplexF64(cxref_re, cxref_im), ComplexF64.(cxs1_re, cxs1_im), ComplexF64.(cxs2_re, cxs2_im)
    x_f, cxref_f, nang_f, cxs1_f, cxs2_f = Float32(x), ComplexF32(cxref), Int32(nang), ComplexF32.(cxs1), ComplexF32.(cxs2)
    return isapprox(
        miemfp.bhmie(x, cxref, nang),
        FFIWraps.bhmie_fortran77(x_f, cxref_f, nang_f, cxs1_f, cxs2_f),
        rtol=0.1,
    )
end

f64_gen = PropCheck.itype(Float64)
UInt32_gen = PropCheck.itype(UInt32)
f64_vec_gen = PropCheck.vector(isample(1:100), f64_gen)
bhmie_gen = PropCheck.interleave(
    f64_gen,
    f64_gen,
    f64_gen,
    UInt32_gen,
    f64_vec_gen,
    f64_vec_gen,
    f64_vec_gen,
    f64_vec_gen,
)

@testset "bhmie" begin
    @testset "miemfp.bhmie" begin
        c_check = PropCheck.check(julia_vs_c, bhmie_gen)
        c_result = c_check == true
        if !c_result
            println("Fail vs C, PropCheck:")
            display(c_check)
        end
        @test c_result
        f_check = PropCheck.check(julia_vs_fortran, bhmie_gen)
        f_result = f_check == true
        if !f_result
            println("Fail vs Fortran, PropCheck:")
            display(f_check)
        end
        @test f_result
        f77_check = PropCheck.check(julia_vs_fortran77, bhmie_gen)
        f77_result = f77_check == true
        if !f77_result
            println("Fail vs Fortran77, PropCheck:")
            display(f77_check)
        end
        @test f77_result
    end
end