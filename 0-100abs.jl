include("nanoconc.jl")
using PyPlot

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
predict(x) = nanoconc.abspredict(1.333,250.,900.,0x00001000,x,mat,1000.,1.)
spectra = predict.(particledata)

for spectrum in spectra
    plot(spectrum[:,1],spectrum[:,2])
end
show()
