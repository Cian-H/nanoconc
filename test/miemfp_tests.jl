using Test
using Random
using PropCheck
using Debugger
using PyCall


include("../anchors.jl")

import .Anchors: TEST_DIR, SRC_DIR, ROOT_DIR

if !@isdefined TestUtils
    include(joinpath(TEST_DIR, "testutils.jl"))
end
if !@isdefined miemfp
    include(joinpath(SRC_DIR, "miemfp.jl"))
end
if !@isdefined FFIWraps
    include(joinpath(TEST_DIR, "ffi_wraps.jl"))
end


# Set up the Python environment
run(`$ROOT_DIR/setup_venv.sh`)
ENV["PYTHON"] = joinpath(ROOT_DIR, ".venv/bin/python")

@pyinclude(joinpath(TEST_DIR, "miemfp_tests.py"))


miemfp.bhmie(
    x::Float64,
    cxref::ComplexF64,
    nang::Int64,
    s1::Vector{ComplexF64},
    s2::Vector{ComplexF64},
) = miemfp.bhmie(x, cxref, UInt32(nang))

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

@testset "miemfp" begin
    @testset "miemfp.bhmie" begin
        event1, event2 = py"asyncio.Event"(), py"asyncio.Event"()
        event1.set(), event2.set()
        result, output = py"compare_bhmie_functions"(miemfp.bhmie, FFIWraps.bhmie_fortran, event1, event2)
        @test result
        result, output = py"compare_bhmie_functions"(miemfp.bhmie, FFIWraps.bhmie_fortran77, event1, event2)
        @test result
        result, output = py"compare_bhmie_functions"(miemfp.bhmie, FFIWraps.bhmie_c, event1, event2)
        @test result
    end
end