"""
A module that uses the functions contained in miemfp and quantumcalc to find the
concentration of a nanoparticle colloid given enough information
"""
@fastmath module nanoconc

export fetchinfo, dispinfo, addmaterial, loadmaterial, listmaterials, delmaterial, editmaterial, qpredict, abspredict

include("miemfp.jl")
include("quantumcalc.jl")

using CSV
using DataFrames
using HDF5

const matfile = "data/materials.h5" # this is the location of the material datafile

"""
Fetches info string for a saved material
"""
function fetchinfo(material::String)::String
    return string("\nMaterial:\t\t", material,
        "\n- Valence:\t\t", h5readattr(matfile, material)["z"],
        "\n- Atomic Mass:\t\t", h5readattr(matfile, material)["am"],
        "amu\n- Density:\t\t", h5readattr(matfile, material)["rho"],
        "g/cm^3\n- Resistivity:\t\t", h5readattr(matfile, material)["res"],
        "g/cm^3\nDerived values:\n- Plasma frequency:\t", h5readattr(matfile, material)["omp"],
        "e14Hz\n- Collision frequency:\t", h5readattr(matfile, material)["om0"],
        "e14Hz\n- Fermi Velocity:\t", h5readattr(matfile, material)["fv"], "cm/s\n",
        h5readattr(matfile, material)["info"],
        "\nDescription:\n\"", h5readattr(matfile, material)["description"],
        "\"\n")
end

"""
Displays info for saved material
"""
function dispinfo(material::String)
    print(fetchinfo(material))
end

"""
stores material data for future use in the materials HDF5 data file with the given
material name and description. Requires valence (z), atomic mass (am), density (rho)(g/cm^3)
resistivity (res)(ohm*m) and a complex refractive index spectrum from a file at filepath.
Refractive index data is taken from an appropriately formatted csv file (formatted with
columns entitled "w","n" and "k" for wavelength in nm, n and k respectively)
"""
function addmaterial(z::Float64, am::Float64, rho::Float64, res::Float64,
    filepath::String, material::String, description::String;
    disp::Bool=true)
    try
        flag = h5open(matfile, "r") do file
            has(file, material)
        end
    catch
        flag = false
    end

    if !flag
        df::DataFrames.DataFrame = CSV.read(filepath)
        data::Array{Float64,2} = (convert(Array{Float64,2}, hcat(df[:w], df[:n], df[:k])))
        h5write(matfile, material, data)
        omp, om0, fv = quantumcalc.mieparams(z, am, rho, res)
        println(omp)
        println(om0)
        println(fv)
        h5writeattr(matfile, material, Dict(
            "z" => z,
            "am" => am,
            "rho" => rho,
            "res" => res,
            "omp" => omp,
            "om0" => om0,
            "fv" => fv,
            "description" => description,
            "info" => string("Wavelength Range:\t", minimum(data[:, 1]),
                "nm-", maximum(data[:, 1]),
                "nm\nNumber of datapoints:\t",
                size(data)[1])))
        if disp == true
            println("\nAdded material:")
            dispinfo(material)
        end
    elseif disp == true
        println("\nMaterial ", material, " Already exists!\n\nExisting material:")
        dispinfo(material)
    end
end

"""
An alternative version of the addmaterial function for materials with known mie
    parameters. Requires input of plasma freq (omp) in Hz/1e14, collision freq (om0) in Hz/1e14 and
    Fermi velocity (fv) in cm/s
"""
function addmaterial(omp::Float64, om0::Float64, fv::Float64,
    filepath::String, material::String, description::String;
    disp::Bool=true)
    try
        flag = h5open(matfile, "r") do file
            has(file, material)
        end
    catch
        flag = false
    end

    if !flag
        df::DataFrames.DataFrame = CSV.read(filepath)
        data::Array{Float64,2} = (convert(Array{Float64,2}, hcat(df[:w], df[:n], df[:k])))
        h5write(matfile, material, data)
        println(omp)
        println(om0)
        println(fv)
        h5writeattr(matfile, material, Dict(
            "omp" => omp,
            "om0" => om0,
            "fv" => fv,
            "description" => description,
            "info" => string("Wavelength Range:\t", minimum(data[:, 1]),
                "nm-", maximum(data[:, 1]),
                "nm\nNumber of datapoints:\t",
                size(data)[1])))
        if disp == true
            println("\nAdded material:")
            dispinfo(material)
        end
    elseif disp == true
        println("\nMaterial ", material, " Already exists!\n\nExisting material:")
        dispinfo(material)
    end
end

"""
This function loads the stored miemfp data for a saved material. omp, om0, fv and
refractive index array are returned as a static tuple that can be readily splatted
into the last few args of the miemfp.qbare and miemfp.qcoat functions
"""
function loadmaterial(material::String; disp::Bool=true
)::Tuple{Float64,Float64,Float64,Array{Float64,2}}
    local output::Tuple{Float64,Float64,Float64,Array{Float64,2}}
    h5open(matfile, "r") do file
        material_data = file[material]
        output = (
            read(material_data["omp"]),
            read(material_data["om0"]),
            read(material_data["fv"]),
            read(material_data)
        )
    end
    if disp == true
        println("\nLoading material:")
        dispinfo(material)
        println()
    end
    return output
end

"""
This function lists all materials for which refractive index data is stored
"""
function listmaterials()
    h5open(matfile, "r") do file
        flag::Bool = false
        for i in file
            flag = true
            break
        end
        if flag == true
            println("\nSaved materials:")
            for material in file
                println(name(material))
            end
        else
            println("\nNo materials saved!")
        end
    end
end

"""
This function deletes a material from the materials datafile
"""
function delmaterial(material::String; disp::Bool=true)
    h5open(matfile, "r+") do file
        if exists(file, material)
            if disp == true
                println("\nDeleting material:")
                dispinfo(material)
                println()
            end
            o_delete(file, material)
        else
            println("\nMaterial \"", material, "\" not found!")
        end
    end
end

"""
This function allows editing of specific parameters for a stored material
"""
function editmaterial(material::String, param::String, value; disp::Bool=true)
    h5open(matfile, "r+") do file
        if !exists(file, material)
            if disp
                println("\nMaterial \"", material, "\" not found!")
            end
            return
        end
    end
    matattr = h5readattr(matfile, material)
    matdata = h5read(matfile, material)
    if param == "refind"
        val0 = matdata
        value::typeof(matdata)
        matdata = value
    else
        val0 = matattr[param]
        value::typeof(val0)
        matattr[param] = value
    end
    h5open(matfile, "r+") do file
        o_delete(file, material)
    end
    h5write(matfile, material, matdata)
    h5writeattr(matfile, material, matattr)
    if disp
        println(material, " Parameter ", param, " changed from ", val0, " to ", value)
    end
end

"""
A function that predicts the expected average per particle extinction efficiency
spectrum for a colloid within a wavelength range (wavel1 - wavel2) with a set
number of points (numval) per spectrum. Also requires refractive index of the
medium, a particledata array containing the diameters of the particles in solution
(in nm) and their relative amounts, and the saved material data array
"""
function qpredict(refmed::Float64, wavel1::Float64, wavel2::Float64,
    numval::UInt32, particledata::Array{Float64,2},
    materialdata::Tuple{Float64,Float64,Float64,Array{Float64,2}}
)::Array{Float64,2}
    # normalise relative amounts of particles to a total of 1
    particledata[:, 2] = particledata[:, 2] ./ sum(particledata[:, 2])
    # convert particle diameters in particledata to radii
    particledata[:, 1] = particledata[:, 1] ./ 2
    # define a wrapper function for later vectorised calculations
    spec(x) = miemfp.qbare(wavel1, wavel2, numval, 0x00000002, refmed, x, materialdata...)
    # note: the UInt32 above for "scangles" must be at least 2 for accurate Qext. Higher "scangles"
    ##	does not increase Qext precision but will be important if modding to also account for Qsca
    # perform vectorised spectrum prediction on particledata array
    spectra = spec.(particledata[:, 1])
    # for each spectrum, multiply all values by their relative amounts from the particledata array
    @inbounds @simd for i in eachindex(spectra)
        spectra[i][:, 2] = particledata[i, 2] .* spectra[i][:, 2]
    end
    # prep an array for final predicted spectrum
    data::Array{Float64,2} = Array{Float64,2}(undef, numval, 2)
    data[:, 1] = spectra[1][:, 1]
    data[:, 2] = zeros(numval)
    # for each weighted spectrum, add it to the total spectrum in the "data" array
    @inbounds @simd for x in spectra
        data[:, 2] = data[:, 2] .+ x[:, 2]
    end
    return data
end

struct q_to_sigma
    geometric_cross_section::Vector{Float64}

    function q_to_sigma(diameter::Vector{Float64})::q_to_sigma
        new(π .* ((diameter ./ 2).^2))
    end 
end

function (p::q_to_sigma)(q::Matrix{Float64})::Array{Float64, 3}
    # Calculate the extinction cross-section
    hcat([q .* x for x in p.geometric_cross_section])
end

function sigmapredict(refmed::Float64, wavel1::Float64, wavel2::Float64,
    numval::UInt32, particledata::Array{Float64,2},
    materialdata::Tuple{Float64,Float64,Float64,Array{Float64,2}}
)::Array{Float64,2}
    # normalise relative amounts of particles to a total of 1
    particledata[:, 2] = particledata[:, 2] ./ sum(particledata[:, 2])
    # convert particle diameters in particledata to radii
    particledata[:, 1] = particledata[:, 1] ./ 2
    # define a wrapper function for later vectorised calculations
    spec(x) = miemfp.qbare(wavel1, wavel2, numval, 0x00000002, refmed, x, materialdata...)
    # note: the UInt32 above for "scangles" must be at least 2 for accurate Qext. Higher "scangles"
    ##	does not increase Qext precision but will be important if modding to also account for Qsca
    # perform vectorised spectrum prediction on particledata array
    spectra = spec.(particledata[:, 1])
    println(size(spectra))
    # for each spectrum, multiply all values by their relative amounts from the particledata array
    @inbounds @simd for i in eachindex(spectra)
        spectra[i][:, 2] = particledata[i, 2] .* spectra[i][:, 2]
    end
    # prep an array for final predicted spectrum
    data::Array{Float64,2} = Array{Float64,2}(undef, numval, 2)
    data[:, 1] = spectra[1][:, 1]
    data[:, 2] = zeros(numval)
    # for each weighted spectrum, add it to the total spectrum in the "data" array
    @inbounds @simd for x in spectra
        data[:, 2] = data[:, 2] .+ x[:, 2]
    end
    return data
end

function abspredict(refmed::Float64, wavel1::Float64, wavel2::Float64,
    numval::UInt32, particledata::Array{Float64,2},
    materialdata::Tuple{Float64,Float64,Float64,Array{Float64,2}},
    ppml::Float64, d0::Float64)::Array{Float64,2}
    data = qpredict(refmed, wavel1, wavel2, numval, particledata, materialdata)
    predict(qext) = quantumcalc.predictabs(
        qext,
        sum([particledata[i, 1] * particledata[i, 2] for i in 1:size(particledata)[1]]) / sum(particledata[:, 2]),
        ppml,
        d0
    )
    data[:, 2] = predict.(data[:, 2])
    return data
end

end
