using Test

if !@isdefined TestUtils
    include("testutils.jl")
end
if !@isdefined nanoconc
    include("../src/nanoconc.jl")
end

@testset "miemfp" begin
    @testset "miemfp.bhmie" begin
        TestUtils.test_from_serialized(
            nanoconc.miemfp.bhmie,
            "test/data/Main.nanoconc.miemfp.bhmie.ser"
        )
    end

    @testset "miemfp.mfp" begin
        TestUtils.test_from_serialized(
            nanoconc.miemfp.mfp,
            "test/data/Main.nanoconc.miemfp.mfp.ser"
        )
    end

    @testset "miemfp.qbare" begin
        TestUtils.test_from_serialized(
            nanoconc.miemfp.qbare,
            "test/data/Main.nanoconc.miemfp.qbare.ser"
        )
    end
end