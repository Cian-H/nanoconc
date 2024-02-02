using Test

if !@isdefined TestUtils
    include("testutils.jl")
end
if !@isdefined nanoconc
    include("../src/nanoconc.jl")
end

@testset "quantumcalc" begin
    @testset "quantumcalc.attencoeff" begin
        TestUtils.test_from_serialized(
            nanoconc.quantumcalc.attencoeff,
            "test/data/Main.nanoconc.quantumcalc.attencoeff.ser"
        )
    end

    @testset "quantumcalc.attentoabs" begin
        TestUtils.test_from_serialized(
            nanoconc.quantumcalc.attentoabs,
            "test/data/Main.nanoconc.quantumcalc.attentoabs.ser"
        )
    end

    @testset "quantumcalc.numtomol" begin
        TestUtils.test_from_serialized(
            nanoconc.quantumcalc.numtomol,
            "test/data/Main.nanoconc.quantumcalc.numtomol.ser"
        )
    end

    @testset "quantumcalc.predictabs" begin
        TestUtils.test_from_serialized(
            nanoconc.quantumcalc.predictabs,
            "test/data/Main.nanoconc.quantumcalc.predictabs.ser"
        )
    end
end