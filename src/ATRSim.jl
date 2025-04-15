module ATRSim
using Integrals

export ATRConfig, Thinfilm, Bulk, E₀², dp, ξ, χ, deff, prefactor

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

"""
    dp(λ::RealsA, config::Thinfilm)
    dp(λ::RealsA, n₂::RealsA, config::Bulk)

Penetration depth as a function of wavelength.
"""
function dp end

function dp(λ::RealsA, config::Thinfilm)
    (; θ, n₁, n₃) = config
    λ₁ = λ ./ n₁
    n₃₁ = n₃ ./ n₁
    f(λ₁, n₃₁) = !ismissing(n₃₁) && sin(θ) > n₃₁ ?
        (λ₁ / (2π * sqrt(sin(θ)^2 - n₃₁^2))) :
        missing # n₁ can have missing elements (if generated from SellmeierEqn)
    f.(λ₁, n₃₁)
end

function dp(λ::RealsA, n₂::RealsA, config::Bulk)
    (; θ, n₁) = config
    λ₁ = λ ./ n₁
    n₂₁ = n₂ ./ n₁
    f(λ₁, n₂₁) = !ismissing(n₂₁) && sin(θ) > n₂₁ ?
        λ₁ / (2π * sqrt(sin(θ)^2 - n₂₁^2)) :
        missing # n₁ can have missing elements (if generated from SellmeierEqn)
    f.(λ₁, n₂₁)
end

"""
    E₀²(n₂::RealsA, config::Thinfilm)
    E₀²(n₂::RealsA, config::Bulk)

Normalized electric field intensity at the IRE-sample interface (Fresnel transmission coefficients). Notation of Harrick and Mirabella are used (as adopted by Arangio et al.).
"""
function E₀² end

function E₀²(n₂::RealsA, config::Thinfilm)
    (; θ, n₁, n₃) = config
    n₂₁ = n₂ ./ n₁
    n₃₁ = n₂ ./ n₁
    n₃₂ = n₃ ./ n₂
    function f(n₃₁, n₃₂)
        if !ismissing(n₃₁) && sin(θ) > n₃₁
            #=
            (s) perpendicular and (p) parallel polarization
            =#
            Eₛ = 2 * cos(θ) / sqrt(1 - n₃₁^2)
            Eₚ = Eₛ * sqrt((1 + n₃₂^4) * sin(θ)^2 - n₃₁^2) /
                sqrt((1 + n₃₁^2) * sin(θ)^2 - n₃₁^2)
            E₀² = (Eₛ^2 + Eₚ^2) / 2
        else
            missing
        end
    end
    f.(n₃₁, n₃₂)
end

function E₀²(n₂::RealsA, config::Bulk)
    (; θ, n₁) = config
    n₂₁ = n₂ ./ n₁
    function f(n₂₁)
        if !ismissing(n₂₁) && sin(θ) > n₂₁
            #=
            (s) perpendicular and (p) parallel polarization
            zo: parallel to plane of incidence and perpendicular to plane of the surface
            xo: parallel to plane of incidence and parallel to plane of the surface
            =#
            Eyo=  2 * cos(θ) / sqrt(1 - n₂₁^2)
            Exo = Eyo * sqrt(sin(θ)^2 - n₂₁^2) / sqrt((1 + n₂₁^2) * sin(θ)^2 - n₂₁^2)
            Ezo = Eyo * sin(θ) / sqrt((1 + n₂₁^2) * sin(θ)^2 - n₂₁^2)
            Eₛ² = Eyo^2
            Eₚ² = Exo^2 + Ezo^2
            E₀² = (Eₛ² + Eₚ²) / 2
        else
            missing
        end
    end
    f.(n₂₁)
end

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
    χ(ν::RealsA, R::Real, config::Thinfilm)

Particle modification factor. Equation S10 is reprogrammed solely in terms of radius `R` instead of function of film thickness `d`, eliminating the dependence on packing density.
"""
function χ(ν::RealsA, R::Real, config::Thinfilm)
    λ = 1e4 ./ ν    
    dₚ = dp(λ, config)
    z = similar(dₚ)
    for i in eachindex(dₚ)
        if ismissing(dₚ[i])
            z[i] = missing
        else
            f(z, p) = (1 - (R - z)^2 / R^2) * exp(-2 * z / dₚ[i]) /
                (4 / 3 * R)
            sol = solve(IntegralProblem(f, (0, 2 * R)), QuadGKJL())
            z[i] = sol.u
        end
    end
    z
end

"""
    prefactor(ν::RealsA, n₂::RealsA, R::Real, config::Thinfilm)

Prefactor to α / ρ .* m. 
"""
function prefactor(ν::RealsA, n₂::RealsA, R::Real, config::Thinfilm)
    (; N, aIRE) = config
    N ./ aIRE .* ξ(n₂, config) .* χ(ν, R, config)
end

function prefactor(n₂::RealsA, config::Thinfilm)
    (; N, aIRE) = config
    N ./ aIRE .* ξ(n₂, config)
end

end # ATRSim
