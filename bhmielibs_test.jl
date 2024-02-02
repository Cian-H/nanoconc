run(`scripts/build_other_impls.sh`)

include("src/miemfp.jl")
# using miemfp

function bhmie_c(x::Float32, cxref::ComplexF32, nang::UInt32, cxs1::Vector{ComplexF32}, cxs2::Vector{ComplexF32})
    # Pre-allocate memory for the output variables
    qext = Ref{Float32}(0.0)
    qsca = Ref{Float32}(0.0)
    qback = Ref{Float32}(0.0)
    gsca = Ref{Float32}(0.0)

    # Ensure cxs1 and cxs2 have proper sizes, as expected by the C function
    # For example, if they need to be of size `nang`, you should verify or resize them accordingly

    # Call the C function
    ccall((:bhmie, ".bhmielibs/bhmie-c/bhmie.so"), Cvoid,
        (Float32, ComplexF32, UInt32, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
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
    ccall((:bhmie_, ".bhmielibs/bhmie-f/bhmie.so"), Cvoid,
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
    ccall((:bhmie_, ".bhmielibs/bhmie-f/bhmie.so"), Cvoid,
          (Ref{Float32}, Ref{ComplexF32}, Ref{Int32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
          x, refrel, nang, s1, s2, qext, qsca, qback, gsca)

    # Return the modified values
    return qext[], qsca[], qback[], gsca[]
end

# Create test data
x = Float32(1.0)  # Example value for x
cxref = ComplexF32(1.5, 0.5)  # Example complex refractive index
nang = UInt32(2)  # Example number of angles

# Example arrays for scattering amplitudes, initialized with dummy complex values
cxs1 = [ComplexF32(0.1, 0.2) for _ in 1:nang]
cxs2 = [ComplexF32(0.3, 0.4) for _ in 1:nang]

# Test C wrapper
qext, qsca, qback, gsca = bhmie_c(x, cxref, nang, cxs1, cxs2)
println("bhmie_c output: qext = $qext, qsca = $qsca, qback = $qback, gsca = $gsca")

# Test Fortran wrapper
qext, qsca, qback, gsca = bhmie_fortran(x, cxref, Int32(nang), cxs1, cxs2)
println("bhmie_fortran output: qext = $qext, qsca = $qsca, qback = $qback, gsca = $gsca")

# Test Fortran77 wrapper
qext, qsca, qback, gsca = bhmie_fortran77(x, cxref, Int32(nang), cxs1, cxs2)
println("bhmie_fortran77 output: qext = $qext, qsca = $qsca, qback = $qback, gsca = $gsca")

# Test Julia wrapper
qext, qsca, qback, gsca = miemfp.bhmie(Float64(x), ComplexF64(cxref), nang)
println("bhmie_julia output: qext = $qext, qsca = $qsca, qback = $qback, gsca = $gsca")