using Test

if !@isdefined TestUtils
    include("testutils.jl")
end
if !@isdefined nanoconc
    include("../src/nanoconc.jl")
end

@testset "nanoconc" begin
    @testset "qpredict" begin
        TestUtils.test_from_serialized(
            nanoconc.qpredict,
            "test/data/Main.nanoconc.qpredict.ser"
        )
    end

    @testset "nanoconc" begin
        TestUtils.test_from_serialized(
            nanoconc.abspredict,
            "test/data/Main.nanoconc.abspredict.ser"
        )
    end
end