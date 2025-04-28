module ATRSim
using Unzip
using Integrals

export ATRConfig, Thinfilm, Bulk, E₀², t², ξF, dp, ξ, χ, deff, prefactor

const RealsC = Union{Real, Vector{<:Union{Missing, Real}}}
const RealsA = Union{Real, AbstractVector{<:Union{Missing, Real}}}

"""
`Thinfilm` and `Bulk` are structures that store information about the ATR configuration. Their type information is used for method dispatch. When a vector is passed for n₁ or n₃, vectors such as ν and n₂ passed to its methods must be of the same length.
"""
abstract type ATRConfig end

## Thin film
struct Thinfilm <:ATRConfig
    N::Integer
    aIRE::Float64
    θ::Float64
    n₁::RealsC
    n₃::RealsC
end

## Semi-infinite bulk
struct Bulk <:ATRConfig
    N::Integer
    aIRE::Float64
    θ::Float64
    n₁::RealsC
end

ATRConfig(N::Integer, aIRE::Float64, θ::Float64, n₁::RealsA, n₃::RealsA) = 
    Thinfilm(N, aIRE, θ, n₁, n₃)

ATRConfig(N::Integer, aIRE::Float64, θ::Float64, n₁::RealsA) = 
    Bulk(N, aIRE, θ, n₁)

#=
n₁ can have missing elements (if generated from SellmeierEqn)
=#

## -----------------------------------------------------------------------------

unpack(x) = isa(x, Tuple) ? x : unzip(x)

## -----------------------------------------------------------------------------

"""
    dp(λ::RealsA, n₂::RealsA, config::ATRConfig)
    dp(θ::Real, λ::Real, n₁::Real, n₂::Real)

Penetration depth as a function of wavelength.
"""
function dp end

function dp(λ::RealsA, n₂::RealsA, config::ATRConfig)
    (; θ, n₁) = config
    dp.(θ, λ, n₁, n₂)
end

function dp(θ::Real, λ::Real, n₁::Real, n₂::Real)
    λ₁ = λ / n₁
    n₂₁ = n₂ / n₁
    sin(θ) < n₂₁ ? NaN : λ₁ / (2π * sqrt(sin(θ)^2 - n₂₁^2))
end

"""
    E₀²(n₂::RealsA, config::Thinfilm)
    E₀²(n₂::RealsA, config::Bulk)
    E₀²(θ::Real, n₁::Real, n₂::Real, n₃::Real)
    E₀²(θ::Real, n₁::Real, n₂::Real)

Normalized electric field intensity at the IRE-sample interface (Fresnel transmission coefficients). Notation of Harrick and Mirabella are used (as adopted by Arangio et al.).
"""
function E₀² end

function E₀²(n₂::RealsA, config::Thinfilm)
    (; θ, n₁, n₃) = config
    (Eₛ, Eₚ) = unpack(E₀².(θ, n₁, n₂, n₃))
    Eₛ² = @. abs2(Eₛ)
    Eₚ² = @. abs2(Eₚ)
    @. (Eₛ² + Eₚ²) / 2
end

function E₀²(n₂::RealsA, config::Bulk)
    (; θ, n₁) = config
    (Exo, Eyo, Ezo) = unpack(E₀².(θ, n₁, n₂))
    Eₛ² = @. abs2(Eyo)
    Eₚ² = @. abs2(Exo) + abs2(Ezo)
    @. (Eₛ² + Eₚ²) / 2
end

function E₀²(θ::Real, n₁::Real, n₂::Real, n₃::Real)
    n₂₁ = n₂ / n₁
    n₃₁ = n₃ / n₁
    n₃₂ = n₃ / n₂
    #=
    (s) perpendicular and (p) parallel polarization
    =#
    Eₛ = 2 * cos(θ) / sqrt(1 - n₃₁^2)
    Eₚ = Eₛ * sqrt((1 + n₃₂^4) * sin(θ)^2 - n₃₁^2) / sqrt((1 + n₃₁^2) * sin(θ)^2 - n₃₁^2)
    (Eₛ, Eₚ)
end

function E₀²(θ::Real, n₁::Real, n₂::Real)
    n₂₁ = n₂ / n₁
    #=
    Making n₂₁ complex yields unreasonable values.
    =#
    if sin(θ) < n₂₁ return ntuple(_ -> NaN, 3) end
    #=
    (s) perpendicular and (p) parallel polarization
    zo: parallel to plane of incidence and perpendicular to plane of the surface
    xo: parallel to plane of incidence and parallel to plane of the surface
    =#
    Eyo = 2 * cos(θ) / sqrt(1 - n₂₁^2)
    Exo = Eyo * sqrt(sin(θ)^2 - n₂₁^2) / sqrt((1 + n₂₁^2) * sin(θ)^2 - n₂₁^2)
    Ezo = Eyo * sin(θ) / sqrt((1 + n₂₁^2) * sin(θ)^2 - n₂₁^2)
    (Exo, Eyo, Ezo)
end

## -----------------------------------------------------------------------------

"""
    function t²(n₂::RealsA, config::Thinfilm)
    function t²(n₂::RealsA, config::Bulk)
    function t²(θᵢ::Real, nᵢ::Real, nₜ::Real)

Alternative calculation for the normalized electric field intensity (Fresnel transmission coefficients), derived with reference to original Fresnel equations. From Beattie et al., doi:10.1016/S0924-2031(00)00084-9, 2000.

### Parameters
n₂: refractive index of sample (also refered to as the rarer medium or transmitting medium)
"""
function t² end

function t²(n₂::RealsA, config::Thinfilm)
    (; θ, n₁, n₃) = config
    (tₚx, tₛy, tₚz) = unpack(t².(θ, n₁, n₃))
    @. tₚz = tₚz * (n₃ / n₂)^2
    @. (abs2(tₛy) + (abs2(tₚx) + abs2(tₚz))) / 2
end

function t²(n₂::RealsA, config::Bulk)
    (; θ, n₁) = config
    (tₚx, tₛy, tₚz) = unpack(t².(θ, n₁, n₂))
    @. (abs2(tₛy) + (abs2(tₚx) + abs2(tₚz))) / 2
end

function t²(θᵢ::Real, nᵢ::Real, nₜ::Real)
    θₜ = asin(complex(nᵢ / nₜ * sin(θᵢ)))
    tₛy = 2 * nᵢ * cos(θᵢ) / (nᵢ * cos(θᵢ) + nₜ * cos(θₜ))    
    tₚ = nᵢ * cos(θᵢ) / (nₜ * cos(θᵢ) + nᵢ * cos(θₜ))
    tₚx = tₚ * cos(θₜ)
    tₚz = tₚ * sin(θₜ)
    (tₚx, tₛy, tₚz)
end

## -----------------------------------------------------------------------------

"""
    deff(ν::RealsA, n₂::RealsA, config::Bulk)

Effective thickness for thick sample.
"""
function deff(ν::RealsA, n₂::RealsA, config::Bulk)
    (; θ, n₁) = config
    λ = 1e4 ./ ν
    n₂₁ = n₂ ./ n₁
    (n₂₁ .* E₀²(n₂, config) ./ cos(θ)) .* dp(λ, n₂, config) ./ 2
end

"""
    ξ(n₂::RealsA, config::Thinfilm)

Effective thickness divided by actual sample thickness (ξ = deff / d) for thin film.
"""
function ξ(n₂::RealsA, config::Thinfilm)
    #=
    Based on Harrick1979 and Mirabella1985
    To modify n₂₁, see
    doi:10.1039/c3sm52817k and
    doi:10.1016/S0924-2031(00)00084-9
    =#
    (; θ, n₁) = config
    n₂₁ = n₂ ./ n₁
    n₂₁ .* E₀²(n₂, config) ./ cos(θ)
end

"""
    χ(ν::RealsA, n₂::RealsA, R::Real, config::Thinfilm)

Particle modification factor. Equation S10 is reprogrammed solely in terms of radius `R` instead of function of film thickness `d`, eliminating the dependence on packing density.
"""
function χ(ν::RealsA, n₂::RealsA, R::Real, config::Thinfilm)
    λ = 1e4 ./ ν    
    dₚ = dp(λ, n₂, config)
    z = similar(dₚ)
    domain = (0, 2 * R)
    for i in eachindex(dₚ)
        if ismissing(dₚ[i])
            z[i] = missing
            continue
        elseif isnan(dₚ[i])
            z[i] = NaN
            continue            
        end
        f(z, p) = (1 - (R - z)^2 / R^2) * exp(-2 * z / dₚ[i]) / (4 / 3 * R)
        sol = solve(IntegralProblem(f, domain), QuadGKJL())
        z[i] = sol.u
    end
    z
end

"""
    prefactor(ν::RealsA, n₂::RealsA, R::Real, config::Thinfilm)

Prefactor is `z` in `z * α / ρ .* m`.
"""
function prefactor(ν::RealsA, n₂::RealsA, R::Real, config::Thinfilm)
    (; N, aIRE) = config
    N ./ aIRE .* ξ(n₂, config) .* χ(ν, n₂, R, config)
end

function prefactor(n₂::RealsA, config::Thinfilm)
    (; N, aIRE) = config
    N ./ aIRE .* ξ(n₂, config)
end

end # ATRSim
