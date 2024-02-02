println("\nImporting modules...")

include("src/nanoconc.jl")
using CSV
using DataFrames
using Profile, ProfileVega

println("Loading parameters...")

refmed = 1.333 # refractive index of the medium
wavel1, wavel2 = 250., 900. # wavelengths to calculate between
numval = 0x00001000 # number of values to calculate in spectrum
particledata = [2.5 1.;
                 5. 1.;
                 10. 1.;
                 20. 1.;
                 30. 1.;
                 40. 1.;
                 50. 1.;
                 60. 1.;
                 70. 1.;
                 80. 1.;
                 90. 1.;
                 100. 1.]
materialdata = nanoconc.loadmaterial("Gold",disp=false)
ppml = 10000. # particles per ml nanoconc.quantumcalc.numtomol(ppml) -> conc
d0 = 1. # path length
params = (refmed,wavel1,wavel2,numval,particledata,materialdata,ppml,d0)
spectrum = nanoconc.abspredict(params...)

println("Calculating absorbance spectrum...")
print("Calculation time: ")

@time spectrum = nanoconc.abspredict(params...)
CSV.write("abs.csv", DataFrame(spectrum, :auto), writeheader=false)

# @profview nanoconc.abspredict(params...)
