run(`scripts/build_other_impls.sh`)

# using Debugger
using BenchmarkTools

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

# Fixed testing values
nang = UInt32(2)  # Example number of angles

c_result = @benchmark bhmie_c(x, cxref, nang, cxs1, cxs2) setup=(
    x = rand(Float32);
    cxref = rand(ComplexF32);
    nang = UInt32($nang);
    cxs1 = rand(ComplexF32, $nang);
    cxs2 = rand(ComplexF32, $nang);
)

f_result = @benchmark bhmie_fortran(x, cxref, nang, cxs1, cxs2) setup=(
    x = rand(Float32);
    cxref = rand(ComplexF32);
    nang = Int32($nang);
    cxs1 = rand(ComplexF32, $nang);
    cxs2 = rand(ComplexF32, $nang);
)

f77_result = @benchmark bhmie_fortran77(x, cxref, nang, cxs1, cxs2) setup=(
    x = rand(Float32);
    cxref = rand(ComplexF32);
    nang = Int32($nang);
    cxs1 = rand(ComplexF32, $nang);
    cxs2 = rand(ComplexF32, $nang);
)

j_result = @benchmark miemfp.bhmie(Float64(x), ComplexF64(cxref), nang) setup=(
    x = rand(Float32);
    cxref = rand(ComplexF32);
    nang = UInt32($nang);
)

println("\nC Implementation")
display(c_result)
println("\nFortran Implementation")
display(f_result)
println("\nFortran 77 Implementation")
display(f77_result)
println("\nJulia Implementation")
display(j_result)
