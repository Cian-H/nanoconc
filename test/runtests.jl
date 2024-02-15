using Test

include("testutils.jl")
include("../src/nanoconc.jl")

# Set up the Python environment
TestUtils.init_pyenv()

# include("nanoconc_tests.jl")
@testset "miemfp" include("miemfp_tests.jl")
# include("quantumcalc_tests.jl")
TestUtils.singleton_include("benchmarks.jl", :Benchmarks, @__MODULE__)

Benchmarks.bench_vs_ffi()