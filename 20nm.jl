println("\nImporting modules...")

include("nanoconc.jl")
using PyPlot

println("Loading parameters...")

particledata = [20. 1.]

mat = nanoconc.loadmaterial("Gold",disp=false)
predict(x) = nanoconc.abspredict(1.333,450.,650.,0x00001000,x,mat,1000.,1.)

println("Calculating spectra...")
print("Calculation time: ")

@time spectrum = predict(particledata)

println("Finding lambda max...")
#println("\nDebug data:\nfindmax = $(findmax(spectrum[:,2]))\nmaxind = $(findmax(spectrum[:,2])[2])")
lmax = spectrum[findmax(spectrum[:,2])[2],1]
println("lambda max = $(lmax)")

println("Plotting spectra...")

plot(spectrum[:,1],spectrum[:,2])
println("Displaying spectra...")
show()
println("Script complete!")
