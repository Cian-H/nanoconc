using Test
using Random
using PropCheck
using PyCall

if !@isdefined TestUtils
    include("testutils.jl")
end

TestUtils.init_pyenv()
TestUtils.singleton_include("../anchors.jl", :Anchors, @__MODULE__)

import .Anchors: TEST_DIR, SRC_DIR, ROOT_DIR

TestUtils.singleton_include(joinpath(SRC_DIR, "miemfp.jl"), :miemfp, @__MODULE__)
TestUtils.singleton_include(joinpath(TEST_DIR, "ffi_wraps.jl"), :FFIWraps, @__MODULE__)

@pyinclude(joinpath(TEST_DIR, "miemfp_tests.py"))


miemfp.bhmie(
    x::Float64,
    cxref::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
) = miemfp.bhmie(x, cxref, nang)

function miemfp.bhmie(
    x::Float64,
    cxref::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
    event::PyObject,
)
    event.clear()
    result = miemfp.bhmie(x, cxref, nang, s1, s2)
    event.set()
    return result
end

FFIWraps.bhmie_fortran(
    x::Float64,
    refrel::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
) = FFIWraps.bhmie_fortran(
    Float32(x),
    ComplexF32(refrel),
    Int32(nang),
    ComplexF32.(s1),
    ComplexF32.(s2),
)

function FFIWraps.bhmie_fortran(
    x::Float64,
    refrel::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
    event::PyObject,
)
    event.clear()
    result = FFIWraps.bhmie_fortran(x, refrel, nang, s1, s2)
    event.set()
    return result
end

FFIWraps.bhmie_fortran77(
    x::Float64,
    refrel::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
) = FFIWraps.bhmie_fortran77(
    Float32(x),
    ComplexF32(refrel),
    Int32(nang),
    ComplexF32.(s1),
    ComplexF32.(s2),
)

function FFIWraps.bhmie_fortran77(
    x::Float64,
    refrel::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
    event::PyObject,
)
    event.clear()
    result = FFIWraps.bhmie_fortran77(x, refrel, nang, s1, s2)
    event.set()
    return result
end

FFIWraps.bhmie_c(
    x::Float64,
    refrel::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
) = FFIWraps.bhmie_c(
    Float32(x),
    ComplexF32(refrel),
    UInt32(nang),
    ComplexF32.(s1),
    ComplexF32.(s2),
)

function FFIWraps.bhmie_c(
    x::Float64,
    refrel::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
    event::PyObject,
)
    event.clear()
    result = FFIWraps.bhmie_c(x, refrel, nang, s1, s2)
    event.set()
    return result
end

@testset "bhmie" begin
    @testset "Properties vs ffi" begin
        event1, event2 = py"asyncio.Event"(), py"asyncio.Event"()
        event1.set(), event2.set()
        result, output = py"compare_bhmie_functions"(miemfp.bhmie, FFIWraps.bhmie_fortran, event1, event2)
        TestUtils.@namedtest "Fortran" result
        result, output = py"compare_bhmie_functions"(miemfp.bhmie, FFIWraps.bhmie_fortran77, event1, event2)
        TestUtils.@namedtest "FORTRAN77" result
        result, output = py"compare_bhmie_functions"(miemfp.bhmie, FFIWraps.bhmie_c, event1, event2)
        TestUtils.@namedtest "C" result
    end
    @testset "Type stability" begin
        expected_return_type = Tuple{Float64,Float64,Float64,Array{ComplexF64,1},Array{ComplexF64,1}}
        x = 1.0
        cxref = 1.0 + 1.0im
        nang = 2
        s1 = [1.0 + 1.0im for _ in 1:nang]
        s2 = [1.0 + 1.0im for _ in 1:nang]
        out = @inferred(miemfp.bhmie(x, cxref, nang, s1, s2))
        @test typeof(out) === expected_return_type # Check that inferred type is correct
        @test out == miemfp.bhmie(x, cxref, nang, s1, s2) # Check that inferred result matches actual result
    end
end