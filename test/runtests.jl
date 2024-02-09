using Test

include("testutils.jl")
include("../src/nanoconc.jl")

# include("nanoconc_tests.jl")
include("miemfp_tests.jl")
# include("quantumcalc_tests.jl")
include("benchmarks.jl")

Benchmarks.bench_vs_ffi()