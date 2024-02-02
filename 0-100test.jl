println("\nImporting modules...")

include("nanoconc.jl")
using PyPlot
using JLD2

println("Loading parameters...")

particledata = ([2.5 1.],
                [5. 1.],
                [10. 1.],
                [20. 1.],
                [30. 1.],
                [40. 1.],
                [50. 1.],
                [60. 1.],
                [70. 1.],
                [80. 1.],
                [90. 1.],
                [100. 1.])

mat = nanoconc.loadmaterial("Gold",disp=false)
predict(x) = nanoconc.qpredict(1.333,250.,900.,0x00001000,x,mat)
println("Calculating spectra...")
print("Calculation time: ")

@time spectra = predict.(particledata)

println("Saving output...")

println("Plotting spectra...")

for spectrum in spectra
    plot(spectrum[:,1],spectrum[:,2])
end
println("Displaying spectra...")
println("Script complete!")
show()
println()
