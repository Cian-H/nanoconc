"""
A module for calculating the expected extinction coefficients and backscattering
parameters for coated/uncoated particles of a given size/material
"""
@fastmath module miemfp

using Base.Threads
using Memoize
using Interpolations
using StaticArrays

const PI_5_996E3::Float64 = pi * 5.996e3
const TWO_PI::Float64 = 2 * π
const HALF_PI::Float64 = π / 2
const BHCOAT_DEL::Float64 = 1e-8

# Avoiding recomputing amu on every loop saves a lot of unnecessary allocations
const PRECOMPUTED_AMU = @SVector [
    1.0, 0.5403023058681398, -0.4161468365471424, -0.9899924966004454, -0.6536436208636119, 0.28366218546322625,
    0.960170286650366, 0.7539022543433046, -0.14550003380861354, -0.9111302618846769, -0.8390715290764524,
    0.004425697988050785, 0.8438539587324921, 0.9074467814501962, 0.1367372182078336, -0.7596879128588213,
    -0.9576594803233847, -0.27516333805159693, 0.6603167082440802, 0.9887046181866692, 0.40808206181339196,
    -0.5477292602242684, -0.9999608263946371, -0.5328330203333975, 0.424179007336997, 0.9912028118634736,
    0.6469193223286404, -0.2921388087338362, -0.9626058663135666, -0.7480575296890004, 0.15425144988758405,
    0.9147423578045313, 0.8342233605065102, -0.013276747223059479, -0.8485702747846052, -0.9036922050915067,
    -0.12796368962740468, 0.7654140519453434, 0.9550736440472949, 0.26664293235993725, -0.6669380616522619,
    -0.9873392775238264, -0.39998531498835127, 0.5551133015206257, 0.9998433086476912, 0.5253219888177297,
    -0.4321779448847783, -0.9923354691509287, -0.6401443394691997, 0.3005925437436371, 0.9649660284921133,
    0.7421541968137826, -0.16299078079570548, -0.9182827862121189, -0.8293098328631502, 0.022126756261955736,
    0.853220107722584, 0.8998668269691938, 0.11918013544881928, -0.7710802229758452, -0.9524129804151563,
    -0.25810163593826746, 0.6735071623235862, 0.9858965815825497
]

"Computes values of amu where nang > 64"
@memoize function compute_amu(nang::Int64)::Vector{Float64}
    cos.((HALF_PI / (nang - 1)) .* (63:nang-2))
end

"Returns amu values generated as efficiently as possible for a given nang value"
function get_amu(nang::Int64)::Vector{Float64}
    if nang <= 64
        return @views PRECOMPUTED_AMU[1:nang]
    else
        return @views vcat(PRECOMPUTED_AMU, compute_amu(nang))
    end
end

"Helper function for calculating om"
@memoize function _om_calc(wavel::Float64)::Float64
    PI_5_996E3 / wavel
end

"Helper function for calculating om0"
@memoize function _om0r_calc(om0::Float64, fv::Float64, radcor::Float64)::Float64
    om0 + (fv / (radcor * 1e7))
end

"Helper function for calculating om^2"
@memoize function _om_sq_calc(om::Float64)::Float64
    om^2
end

"Helper function for calculating om0^2"
@memoize function _om0_sq_calc(om0::Float64)::Float64
    om0^2
end

"Helper function for calculating omp^2"
@memoize function _omp_sq_calc(omp::Float64)::Float64
    omp^2
end

"Helper function for calculating om0r^2"
@memoize function _om0r_sq_calc(om0r::Float64)::Float64
    om0r^2
end

" function that performs manipulations to om necessary for efficient mfp"
function _om_manipulations(wavel::Float64, fv::Float64, radcor::Float64,
    omp::Float64, om0::Float64)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64}

    om = _om_calc(wavel)
    om0r = _om0r_calc(om0, fv, radcor)
    om_sq = _om_sq_calc(om)
    om0_sq = _om0_sq_calc(om0)
    omp_sq = _omp_sq_calc(omp)
    om0r_sq = _om0r_sq_calc(om0r)
    om_sq_plus_om0_sq = om_sq + om0_sq
    return (om, om0r, om_sq, om0_sq, omp_sq, om0r_sq, om_sq_plus_om0_sq)
end

"""
Corrects refractive index and attenuation coefficient values to account for
the mean free-path effect on the free conductance electrons in the particles
(translated from Haiss et als FORTRAN implementation)
"""
function mfp(fv::Float64, wavel::Float64, radcor::Float64, omp::Float64,
    om0::Float64, rn::Float64, rk::Float64)::Tuple{Float64,Float64,Float64}

    om, om0r, om_sq, _, omp_sq, om0r_sq, om_sq_plus_om0_sq =
        _om_manipulations(wavel, fv, radcor, omp, om0)

    b1 = rn^2 - rk^2 - (1 - (omp_sq / om_sq_plus_om0_sq))
    b2 = 2 * rn * rk - (omp_sq * om0 / (om * om_sq_plus_om0_sq))
    a1r = 1 - (omp_sq / (om_sq + om0r_sq))
    a2r = omp_sq * om0r / (om * (om_sq + om0r_sq))

    r_const = sqrt((a1r / 2 + b1 / 2)^2 + (a2r / 2 + b2 / 2)^2)
    rnr = sqrt((a1r + b1) / 2 + r_const)
    rkr = sqrt(-(a1r + b1) / 2 + r_const)
    return (rnr, rkr, om0r)
end

"""
Solves for Qext accounting for the effect of particle coatings on extinction
and scattering coefficients of particles (translated from Bohren/Huffman's
FORTRAN77 implementation for coated particles)
"""
function bhcoat(x::Float64, y::Float64, rfrel1::ComplexF64, rfrel2::ComplexF64
)::Tuple{Float64,Float64,Float64}
    x1::ComplexF64 = rfrel1 * x
    x2::ComplexF64 = rfrel2 * x
    y2::ComplexF64 = rfrel2 * y
    nstop::UInt64 = UInt64(round(y + 4 * cbrt(y) + 1))
    refrel::ComplexF64 = rfrel2 / rfrel1

    d0x1::ComplexF64 = cot(x1)
    d0x2::ComplexF64 = cot(x2)
    d0y2::ComplexF64 = cot(y2)
    psi0y::Float64, psi1y::Float64 = cos(y), sin(y)
    chi0y::Float64, chi1y::Float64 = -sin(y), cos(y)
    # xi0y::ComplexF64 = ComplexF64(psi0y, -chi0y)
    xi1y::ComplexF64 = ComplexF64(psi1y, -chi1y)
    chi0y2::ComplexF64, chi1y2::ComplexF64 = -sin(y2), cos(y2)
    chi0x2::ComplexF64, chi1x2::ComplexF64 = -sin(x2), cos(x2)

    qsca::Float64, qext::Float64, xback::ComplexF64 = 0.0, 0.0, ComplexF64(0.0, 0.0)
    iflag::Bool = false

    @inbounds for n in (1:nstop)::UnitRange{Int64}
        two_n_minus_1::Float64 = 2 * n - 1
        psiy::Float64 = two_n_minus_1 * psi1y / y - psi0y
        chiy::Float64 = two_n_minus_1 * chi1y / y - chi0y
        xiy::ComplexF64 = ComplexF64(psiy, -chiy)
        d1y2::ComplexF64 = 1.0 / (n / y2 - d0y2) - n / y2

        if iflag == false
            d1x1 = x1 / (n - (d0x1 * x1)) - n / x1
            d1x2 = x1 / (n - (d0x2 * x2)) - n / x2
            chix2 = two_n_minus_1 * chi1x2 / x2 - chi0x2
            chiy2 = two_n_minus_1 * chi1y2 / y2 - chi0y2
            chipx2::ComplexF64 = chi1x2 - n * chix2 / x2
            chipy2 = chi1y2 - n * chiy2 / y2
            rack_denominator = chix2 * d1x2 - chipx2
            ancap::ComplexF64 = ((refrel * d1x1 - d1x2) / (refrel * d1x1 * chix2 - chipx2)) / rack_denominator
            brack = ancap * (chiy2 * d1y2 - chipy2)
            bncap::ComplexF64 = ((refrel * d1x2 - d1x1) / (refrel * chipx2 - d1x1 * chix2)) / rack_denominator
            crack = bncap * (chiy2 * d1y2 - chipy2)
            amess1::ComplexF64 = brack * chipy2
            amess2::ComplexF64 = brack * chiy2
            amess3::ComplexF64 = crack * chipy2
            amess4::ComplexF64 = crack * chiy2
            if !((abs(amess1) > BHCOAT_DEL * abs(d1y2))
                 || (abs(amess2) > BHCOAT_DEL)
                 || (abs(amess3) > BHCOAT_DEL * abs(d1y2))
                 || (abs(amess4) > BHCOAT_DEL))
                brack, crack = ComplexF64(0.0, 0.0), ComplexF64(0.0, 0.0)
                iflag = true
            end
        end
        dnbar::ComplexF64 = (d1y2 - brack * chipy2) / (1.0 - brack * chiy2)
        gnbar::ComplexF64 = (d1y2 - crack * chipy2) / (1.0 - crack * chiy2)
        n_over_y::Float64 = n / y
        xiy_minus_xi1y::ComplexF64 = xiy - xi1y
        psiy_minus_psi1y::Float64 = psiy - psi1y
        an_mult::ComplexF64 = dnbar / rfrel2 + n_over_y
        an::ComplexF64 = (an_mult * psiy_minus_psi1y) / (an_mult * xiy_minus_xi1y)
        bn_mult::ComplexF64 = rfrel2 * gnbar + n_over_y
        bn::ComplexF64 = (bn_mult * psiy_minus_psi1y) / (bn_mult * xiy_minus_xi1y)
        two_n_plus_1::Float64 = 2 * n + 1
        qsca += two_n_plus_1 * (abs(an)^2 + abs(bn)^2)
        xback_sign::Float64 = if iseven(n)
            1.0
        else
            -1.0
        end
        xback += xback_sign * two_n_plus_1 * (an - bn)
        qext += two_n_plus_1 * (real(an) + real(bn))
        psi0y, psi1y = psi1y, psiy
        chi0y, chi1y = chi1y, chiy
        xi1y = psi1y - chi1y * im
        chi0x2, chi1x2 = chi1x2, chix2
        chi0y2, chi1y2 = chi1y2, chiy2
        d0x1, d0x2, d0y2 = d1x1, d1x2, d1y2
    end
    y_sq = y^2
    two_over_y_sq = 2.0 / y_sq
    qsca = two_over_y_sq * qsca
    qext = two_over_y_sq * qext
    qback = real((1 / y_sq) * (xback * conj(xback)))
    return (qext, qsca, qback)
end


"""
Finds scattering parameters for uncoated particles. Heavily modified, but originally
based on a combination of the original FORTRAN77 BHMIE algorithm and a Python
implementation attributed to Herbert Kaiser hosted on ScatterLib (http://scatterlib.wikidot.com/mie)
"""
function bhmie(x::Float64, refrel::ComplexF64, nang::Int64
)::Tuple{Float64,Float64,Float64,Array{ComplexF64,1},Array{ComplexF64,1}}
    y::ComplexF64 = x * refrel
    nstop::Int64 = round(x + 4.0 * cbrt(x) + 2.0)
    nn::Int64 = round(max(nstop, abs(y)) + 14)
    amu::Vector = get_amu(nang)

    y_inv = inv(y)
    d = Vector{ComplexF64}(undef, nn)
    d[nn] = ComplexF64(0.0, 0.0)
    for n in nn:-1:2
        n_over_y = y_inv * n
        d[n-1] = n_over_y - inv(d[n] + n_over_y)
    end

    psi0 = cos(x)
    psi1 = sin(x)
    qsca = 0.0
    # xi0 = ComplexF64(psi0, psi1)
    xi1 = ComplexF64(psi1, -psi0)
    chi1 = psi0
    chi0 = -psi1
    s1_1 = zeros(ComplexF64, nang)
    s1_2 = zeros(ComplexF64, nang)
    s2_1 = zeros(ComplexF64, nang)
    s2_2 = zeros(ComplexF64, nang)
    pi0 = zeros(Float64, nang)
    pi1 = ones(Float64, nang)
    tau = Vector{Float64}(undef, nang)

    @inbounds for n in 1:nstop
        nm1x2 = Float64(2 * n - 1)
        psi = nm1x2 * psi1 / x - psi0
        chi = nm1x2 * chi1 / x - chi0

        xi = psi - (chi * im)
        dn = d[n]
        an_mult::ComplexF64 = dn / refrel + n / x
        an::ComplexF64 = (an_mult * psi - psi1) / (an_mult * xi - xi1)
        bn_mult::ComplexF64 = refrel * dn + n / x
        bn::ComplexF64 = (bn_mult * psi - psi1) / (bn_mult * xi - xi1)

        pi_ = pi1
        np1 = n + 1
        _2np1 = 2 * n + 1

        qsca += _2np1 * ((abs(an)^2) + (abs(bn)^2))

        fn = _2np1 / (n * (np1))
        pi1 = (_2np1 * amu .* pi_ - (n + 1) .* pi0) / n
        
        tau = n .* amu .* pi_ .- np1 .* pi0
        anpi = an * pi_
        antau = an .* tau
        bnpi = bn * pi_
        bntau = bn .* tau
        s1_1 .+= fn .* (anpi .+ bntau)
        s2_1 .+= fn .* (bnpi .+ antau)
        
        #* Optimisation: branchless if statement toggling between .+= and .-= for s1_2 and s2_2
        isodd_factor = Float64(isodd(n)) - Float64(!isodd(n))
        s1_2 .+= isodd_factor .* fn .* (anpi .- bntau)
        s2_2 .+= isodd_factor .* fn .* (bnpi .- antau)

        psi0, chi0 = psi1, chi1
        psi1, chi1 = psi, chi
        xi1 = ComplexF64(psi1, -chi1)
        pi0 = pi_
    end

    @inbounds s1::Vector{ComplexF64} = @views vcat(s1_1, reverse(s1_2))
    @inbounds s2::Vector{ComplexF64} = @views vcat(s2_1, reverse(s2_2))

    x_sq::Float64 = x^2
    four_over_x_sq::Float64 = 4.0 / x_sq
    qsca *= four_over_x_sq / 2
    qext = four_over_x_sq * real(s1[1])
    qback = four_over_x_sq * (abs(s1[2*nang-1])^2)

    return (qext, qsca, qback, s1, s2)
end

"Helper function to calculate the apparent cross section"
@memoize function _calc_apparent_cross_section(radcore::Float64, refmed::Float64)::Float64
    TWO_PI * radcore * refmed
end

# The following is a messy setup to allow for caching of the interpolation objects
# The reason this has been implemented this way is because Interpolations.jl has
# very confusing behavior when it comes to typing and caching. This is a workaround.
struct InterpolationPair
    rn::Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Vector{Float64},Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}},Tuple{Vector{Float64}}},Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}},Interpolations.Throw{Nothing}}
    rk::Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Vector{Float64},Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}},Tuple{Vector{Float64}}},Interpolations.Gridded{Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}},Interpolations.Throw{Nothing}}
end

const _cache = Dict{UInt64,InterpolationPair}()

"Helper function to get the interpolation objects for rn and rk"
function _get_rnrk_interp_objects(refcore::Array{Float64,2})::InterpolationPair
    refcore_hash = hash(refcore)
    if !haskey(_cache, refcore_hash)
        _cache[refcore_hash] = InterpolationPair(LinearInterpolation(refcore[:, 1], refcore[:, 2]), LinearInterpolation(refcore[:, 1], refcore[:, 3]))
    end
    return _cache[refcore_hash]
end

"Helper function to calculate the delta wavelength for a given wavelength range and number of values"
@memoize function _calc_delta_wl(wavel1::Float64, wavel2::Float64, numval::UInt32)::Float64
    (wavel2 - wavel1) / (numval - 1)
end

"Helper function to calculate the size parameter"
@memoize function _calc_size_parameter(apparent_cross_section::Float64, wavel::Float64)::Float64
    apparent_cross_section / wavel
end

"Wrapper for partial dispatch of the bhmie function"
struct _bhmie
    nang::UInt32
end

function (p::_bhmie)(x::Float64, refrel::ComplexF64)::Float64
    first(bhmie(x, refrel, p.nang))::Float64
end

"Wrapper for partial dispatch of the mfp function"
struct _mfp
    fv::Float64
    radcore::Float64
    omp::Float64
    om0::Float64
end

function (p::_mfp)(wavel::Float64, rn::Float64, rk::Float64)::ComplexF64
    refre1, refim1, _ = mfp(p.fv, wavel, p.radcore, p.omp, p.om0, rn, rk)
    return ComplexF64(refre1, refim1)::ComplexF64
end

"""
Calculates the expected extinction coefficient spectrum for particles of a given
size for uncoated particles (modified/translated from Haiss et als FORTRAN implementation)
"""
function qbare(wavel1::Float64, wavel2::Float64, numval::UInt32, scangles::UInt32,
    refmed::Float64, radcore::Float64, omp::Float64, om0::Float64,
    fv::Float64, refcore::Array{Float64,2})::Matrix{Float64}

    # calculate the new wavelengths
    qarray::Matrix{Float64} = let wavelengths::StepRangeLen{Float64} = wavel1:_calc_delta_wl(wavel1, wavel2, numval):wavel2
        hcat(wavelengths, Vector{Float64}(undef, numval))
    end
    interp_pair::InterpolationPair = _get_rnrk_interp_objects(refcore)
    let interp_ref_rn = interp_pair.rn, interp_ref_rk = interp_pair.rk, map_fn = _bhmie(scangles), mfp_fn = _mfp(fv, radcore, omp, om0)
        @inbounds @threads for i in 0x00000001:numval
            wl = qarray[i, 1]
            qarray[i, 2] = map_fn(
                _calc_apparent_cross_section(radcore, refmed) / wl,
                mfp_fn(wl, interp_ref_rn(wl), interp_ref_rk(wl)) / refmed
            )
        end
    end

    return qarray
end

"""
Calculates the expected extinction coefficient spectrum for particles of a given
size for coated particles (translated from Haiss et als FORTRAN implementation)
"""
function qcoat(wavel1::Float64, wavel2::Float64, numval::Int64, refmed::Float64,
    radcor::Float64, refre1::Float64, refim1::Float64,
    radcot::Float64, refre2::Float64, refim2::Float64,
    omp::Float64, om0::Float64, fv::Float64, refcore::Array{Float64,2})
    # TODO: Add static output!
    # refcore is an array containing the spectrum of refractive index vs particle size
    tau::Float64 = (1 / om0) * 1e-14
    fmpinf::Float64 = fv * tau
    delta::Float64 = (wavel2 - wavel1) / (numval - 1)

    wavrnrk::Vector{Tuple{Float64,Float64,Float64}} = Vector{Tuple{Float64,Float64,Float64}}()
    k::UInt32 = 1
    @simd for i in 1:numval
        wavel::Float64 = wavel1 + (i - 1) * delta
        # find values for refractive index of particles using linear interpolation:
        if (wavel != refcore[1, 1])
            while (refcore[k, 1] < wavel)
                k += 1
            end
            push!(wavrnrk, (wavel,
                (refcore[k-1, 2] + (wavel - refcore[k-1, 1]) / (refcore[k, 1] - refcore[k-1, 1]) * (refcore[k, 2] - refcore[k-1, 2])),
                (refcore[k-1, 3] + (wavel - refcore[k-1, 1]) / (refcore[k, 1] - refcore[k-1, 1]) * (refcore[k, 3] - refcore[k-1, 3])))
            )
        end
    end
    qarray::Array{Float64} = Array{Float64}(0, 4)
    @inbounds @simd for vals in wavrnrk
        wavel::Float64, rn::Float64, rk::Float64 = vals
        # adjust for "mean free path effect" on free electrons in metal core
        refre1, refim1, om0r = mfp(fv, wavel, radcor, omp, om0, rn, rk)
        # calculate relative size parameters x and y for core and coating respectively
        two_pi_times_refractive_ratio::Float64 = TWO_PI * (refmed / refre1)
        x = radcor * two_pi_times_refractive_ratio
        y = radcot * two_pi_times_refractive_ratio
        # calculate ComplexF64 refractive indices:
        rfrel1 = refre1 + (refim1 * im) / refmed
        rfrel2 = refre2 + (refim2 * im) / refmed
        # appends data to qarray as a row in the format [wavel qext qsca qback]
        # bhcoat outputs are adjusted to be relative to core instead of coating
        qarray = vcat(qarray, [wavel (bhcoat(x, y, rfrel1, rfrel2) .* ((radcot / radcor)^2))...])
    end
    return qarray
end

end
