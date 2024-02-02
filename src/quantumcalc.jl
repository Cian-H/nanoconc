# NOTE: Future functionality to add
# calculation of plasma frequency/collision frequency for non-metals
#       can be done based on electron's "effective mass" and charge
#       carrier density instead of e- density?

"""
A module containing functions that apply various quantum physics formulae
for use with "miemfp.jl" in spectrum prediction/concentration calculation
"""
@fastmath module quantumcalc

using Memoize

# Planck's constant in JS according to the CODATA 2014 definition
const h::Float64 = 6.626070040e-34
# Reduced Planck's Constant
const hbar::Float64 = h / (2 * pi)
# Fine structure constant based on the CODATA 2014 definition
const alpha::Float64 = 7.2973525664e-3
# Speed of light in vacuum by the CODATA 2014 definition
const c::Float64 = 299752458
# Elementary charge in C according to the CODATA 2014 definition
const ec::Float64 = 1.6021766208e-19
# Boltzmann Constant in J/K according to CODATA 2014 definition
const kb::Float64 = 1.38064852e-23
# Electron mass in kg according to the 2014 CODATA definition of 1amu and ref: DOI:10.1038/nature13026
const me::Float64 = 9.109383555654034e-31
# Bohr radius in m
const a0::Float64 = hbar / (me * c * alpha)
# Avogadro's constant in mol^-1 according to the CODATA 2014 definition
const A::Float64 = 6.022140857e23
# Precomputed constants for optimization of functions
const A_E6::Float64 = A * 1e6
const PI_SQUARED::Float64 = pi^2
const THREE_PI_SQUARED::Float64 = 3 * PI_SQUARED
const TWO_PI::Float64 = 2 * pi
const FOUR_PI::Float64 = 2 * TWO_PI
const TWO_PI_TIMES_11_44E15::Float64 = TWO_PI * 11.44e15
const EC_SQUARED::Float64 = ec^2
const A_PI_LOG10_EULER_OVER_1000::Float64 = (A * pi * log10(ℯ)) / 1000
const NUM_TO_MOL_CONV::Float64 = 1000 / A


"""
Calculates the free electron density (ne) in m^-3 of an elemental material given
its valence (z), atomic mass in amu (am) and mass density (rho) in g/cm^2
Source: ISBN 0-03-083993-9 (page 4)
"""
function metalfed(z::Float64, rho::Float64, am::Float64)::Float64
    (A_E6 * ((z * rho) / am))
end

"""
Calculates the Fermi wave vector (kf) of electrons in a material given
electron density (ne)
Source: ISBN 0-03-083993-9 (page 36)
"""
function fermivec(ne::Float64)::Float64
    cbrt(ne * THREE_PI_SQUARED)
end

"""
Calculates the Fermi velocity (fv) in m/s of electrons in a material given Fermi
wave vector (kf)
Source: ISBN 0-03-083993-9 (page 36)
"""
function fermivel(kf::Float64; mstar::Float64=me)::Float64
    (hbar / mstar) * kf
end

"""
Calculates the electron sphere radius (rs) given its free electron density (ne)
Source: ISBN 0-03-083993-9 (page 4)
"""
function spherad(ne::Float64)::Float64
    cbrt(3 / (FOUR_PI * ne))
end

"""
Calculates the plasma frequency (omp) of a material in Hz give its electron
sphere radius (rs)
Source: ISBN 0-03-083993-9 (page 758)
"""
function plasmafreq(rs::Float64)::Float64
    TWO_PI_TIMES_11_44E15 * ((rs / a0)^(-1.5))
end

"""
Calculates the collision frequency (om0) in Hz for a material given its free
electron density (ne) and its resistivity (res)
Source: ISBN 0-03-083993-9 (page 8)
"""
function collfreq(ne::Float64, res::Float64; mstar::Float64=me)::Float64
    r((EC_SQUARED) * ne * res) / mstar
end

"""
Calculates the plasma freq (omp) in Hz/1e14, collision freq (om0) in Hz/1e14 and
Fermi velocity (fv) in cm/s of an elemental material given its valence (z),
atomic mass (am) in amu, mass density (rho) in g/cm^2 and its resisitivity (res)
in Ohm*m as a tuple of (omp,om0,fv) for use in Mie model algorithms contained in
the "miemfp" module
"""
function mieparams(z::Float64, am::Float64, rho::Float64, res::Float64)
    ne::Float64 = metalfed(z, rho, am)::Tuple{Float64,Float64,Float64}
    @sync begin
        @async omp = plasmafreq(spherad(ne)) * 1e-14
        @async om0 = collfreq(ne, res) * 1e-14
        @async fv = fermivel(fermivec(ne)) * 1e2
    end

    (omp, om0, fv)
end

"""
Converts Qext to the molar attenuation coefficient (κ) based on the average radius
of the particles present. Formula sourced from ISSN 1061-933X
formula: κ = π * R^2 * Qext * Na * log(e) / 1000
"""
function attencoeff(qext::Float64, r::Float64)::Float64
    (r^2) * qext * A_PI_LOG10_EULER_OVER_1000
end

"""
Converts units per millilitre to moles per litre
"""
@memoize function numtomol(n::Float64)::Float64
    n * NUM_TO_MOL_CONV
end

"""
Finds absorbance based on molar attenutaion coefficient (κ), concentration (molar)
(c) and path length (in cm) (d0) based on the Beer/Lambert law
"""
function attentoabs(kappa::Float64, c::Float64, d0::Float64)::Float64
    kappa * c * d0
end

"""
Predicts the abs value as displayed on a spectrum from the extinction efficiency
(qext), average particle radius (avgr), particles per ml (ppml) and path length (l)
"""
function predictabs(qext::Float64, avgr::Float64, ppml::Float64, l::Float64)::Float64
    log10(attentoabs(attencoeff(qext, avgr), numtomol(ppml), l))
end

# Below is an old, failed version of the predictabs code above, kept in case useful in later debugging
#
# @doc """
# This function takes an absorbance (Abs), extinction coefficient (qext), average
# particle radius in nm (r) and path length (d0) and returns a count of particles
# per ml
# """ ->
# function partperml(Abs::Float64,qext::Float64,r::Float64,d0::Float64)
#     return (log(10) * Abs) / (pi * (r ^ 2) *qext * d0) ::Float64
# end
#
# @doc """
# This function takes number of particles per ml, extinction coefficient (qext),
# average particle radius in nm (r) and path length (d0) and returns a predicted
# absorbance value (from DOI: 10.1016/j.saa.2017.10.047)
# """ ->
# function abspredict(N::Float64,qext::Float64,r::Float64,d0::Float64)
#     # note: log10 here is not part of formula but it is required to convert to AU (i think)
#     return log10((pi * (r^2) *qext * d0 * N) / log(10)) ::Float64
# end

############## FUNCTIONS BELOW THIS POINT ARE PROTOTYPES #######################

# NOTE: if successful change mieparams to miemetal
function XXmieionic(z::Float64, rho::Float64, mm::Float64, shc::Float64, res::Float64)::Tuple{Float64,Float64,Float64}
    # mm is molar mass in g/mol
    ne::Float64 = metalfed(z, rho, mm)
    mstar::Float64 = (shc * hbar^2 * (3 * PI_SQUARED * ne)^(2 / 3)) / (PI_SQUARED * kb^2 * A * z)

    return (
        plasmafreq(spherad(ne)) * 1e-14,
        collfreq(ne, res, mstar=mstar) * 1e-14,
        fermivel(fermivec(ne), mstar=mstar) * 1e2
    )
end

function mieionic(z::Float64, rho::Float64, mm::Float64, shc::Float64, res::Float64)::Tuple{Float64,Float64,Float64}
    # mm is molar mass in g/mol
    ne::Float64 = metalfed(z, rho, mm)
    rs::Float64 = spherad(ne)
    mstar::Float64 = shc * me / 0.169 * z * ((rs / a0)^2) * 1e-4

    return (
        plasmafreq(rs) * 1e-14,
        collfreq(ne, res, mstar=mstar) * 1e-14,
        fermivel(fermivec(ne), mstar=mstar) * 1e2
    )
end

end
