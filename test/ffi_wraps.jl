module FFIWraps

include("../anchors.jl")

if !@isdefined TEST_DIR
    include("../anchors.jl")
    import .Anchors: TEST_DIR
end

BHMIELIBS_DIR = joinpath(TEST_DIR, ".cache/bhmie-libs/")

function __init__()
    build_script = joinpath(TEST_DIR, "build_ffi.sh")

    mkpath(BHMIELIBS_DIR)

    run(`$build_script $BHMIELIBS_DIR`)
end

function bhmie_c(x::Float32, cxref::ComplexF32, nang::UInt32, cxs1::Vector{ComplexF32}, cxs2::Vector{ComplexF32})
    # Pre-allocate memory for the output variables
    qext = Ref{Float32}(0.0)
    qsca = Ref{Float32}(0.0)
    qback = Ref{Float32}(0.0)
    gsca = Ref{Float32}(0.0)

    # Ensure cxs1 and cxs2 have proper sizes, as expected by the C function
    # For example, if they need to be of size `nang`, you should verify or resize them accordingly

    # Call the C function
    ccall((:bhmie, "$BHMIELIBS_DIR/bhmie-c/bhmie.so"), Cvoid,
        (Float32, ComplexF32, UInt32, Vector{ComplexF32}, Vector{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
        x, cxref, nang, cxs1, cxs2, qext, qsca, qback, gsca)

    # Return the output variables
    return qext[], qsca[], qback[], gsca[]
end

function bhmie_fortran(x::Float32, refrel::ComplexF32, nang::Int32, s1::Vector{ComplexF32}, s2::Vector{ComplexF32})
    # Pre-allocate output variables
    qext = Ref{Float32}(0.0)
    qsca = Ref{Float32}(0.0)
    qback = Ref{Float32}(0.0)
    gsca = Ref{Float32}(0.0)

    # Call the Fortran subroutine
    ccall((:bhmie_, "$BHMIELIBS_DIR/bhmie-f/bhmie.so"), Cvoid,
          (Ref{Float32}, Ref{ComplexF32}, Ref{Int32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
          x, refrel, nang, s1, s2, qext, qsca, qback, gsca)

    # Return the modified values
    return qext[], qsca[], qback[], gsca[]
end

function bhmie_fortran77(x::Float32, refrel::ComplexF32, nang::Int32, s1::Vector{ComplexF32}, s2::Vector{ComplexF32})
    # Pre-allocate output variables
    qext = Ref{Float32}(0.0)
    qsca = Ref{Float32}(0.0)
    qback = Ref{Float32}(0.0)
    gsca = Ref{Float32}(0.0)

    # Call the Fortran subroutine
    ccall((:bhmie_, "$BHMIELIBS_DIR/bhmie-f/bhmie.so"), Cvoid,
          (Ref{Float32}, Ref{ComplexF32}, Ref{Int32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
          x, refrel, nang, s1, s2, qext, qsca, qback, gsca)

    # Return the modified values
    return qext[], qsca[], qback[], gsca[]
end

end