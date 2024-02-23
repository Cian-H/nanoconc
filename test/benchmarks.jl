module Benchmarks

include("../anchors.jl")
include("ffi_wraps.jl")

import .Anchors.SRC_DIR
import .FFIWraps: bhmie_c, bhmie_fortran, bhmie_fortran77
using BenchmarkTools
using InteractiveUtils #! DEBUG

include("$SRC_DIR/miemfp.jl")

function bench_vs_ffi()
    # Fixed testing values
    nang = 2  # Example number of angles

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
    
    j_result = @benchmark miemfp.bhmie(Float64(x), ComplexF64(cxref), Int64(nang)) setup=(
        x = rand(Float32);
        cxref = rand(ComplexF32);
        nang = $nang;
    )

    println("\nC Implementation")
    display(c_result)
    println("\nFortran Implementation")
    display(f_result)
    println("\nFortran 77 Implementation")
    display(f77_result)
    println("\nJulia Implementation")
    display(j_result)

    return c_result, f_result, f77_result, j_result
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    result = Benchmarks.bench_vs_ffi()

    include("../anchors.jl")
    import .Anchors.ROOT_DIR
    using Pkg

    current_package_version = Pkg.TOML.parsefile("$ROOT_DIR/Project.toml")["version"]

    function display_to_file(io, x)
        show(IOContext(io, :limit => false, :color => true), "text/plain", x)
    end

    open("$ROOT_DIR/benchmarks/$current_package_version.ansi", "w") do io
        println(io, "C Implementation")
        display_to_file(io, result[1])
        println(io, "\n\nFortran Implementation")
        display_to_file(io, result[2])
        println(io, "\n\nFortran 77 Implementation")
        display_to_file(io, result[3])
        println(io, "\n\nJulia Implementation")
        display_to_file(io, result[4])
    end
end