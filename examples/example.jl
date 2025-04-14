using DataFrames
using DataFramesMeta
using CSV
using Plots
using Printf
using SellmeierEqn
using ATRSim

pathname = joinpath(@__DIR__, "examples", "Myers2020_ammsulf.csv")
(ν, n, k) = let ammsulf = @subset CSV.read(pathname, DataFrame; comment="#") @. (650 <= :wavenum <= 4000)
    @with ammsulf (:wavenum, :n, :k)
end

α = @. 1 / log(10) * 4 * π * ν * k * 1e2 # m^-1

λ = 1e4 ./ ν # μm
n₁ = refidx(λ, :ZnSe)
N = 10
aIRE = 8.0 * 1e-4 # m^2
θ = 45 / 180 * π

configfilm = ATRConfig(N, aIRE, θ, n₁, 1.0)
configfixed = ATRConfig(N, aIRE, θ, refidx(1e4 /2000, :ZnSe), 1.0)
configbulk = ATRConfig(N, aIRE, θ, n₁)

#| output: true
println([typeof(configfilm), typeof(configfixed), typeof(configbulk)])

d = 0.05 # μm
R = 0.08 # μm
n₂fixed = 1.5

#| output: true
plot(xflip = true, legend = :topleft, xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance")
plot!(ν, α .* d .* ξ(n, configfilm) .* χ(ν, R, configfilm) .* 1e-6, label = "varying n₂")
plot!(ν, α .* d .* ξ(n₂fixed, configfixed) .* χ(ν, R, configfixed) .* 1e-6, label = "fixed n₂")

#| output: true
plot(xflip = true, legend = :topleft, xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance")
plot!(ν, α .* d .* ξ(n₂fixed, configfixed) .* 1e-6, label = "thin film")
plot!(ν, α .* d .* ξ(n₂fixed, configfixed) .* χ(ν, R, configfixed) .* 1e-6, label = "particle")
plot!(ν, 0.1 .* α .* deff(ν, n₂fixed, configbulk) .* 1e-6, label = "0.1 × bulk")

ηₕ = π * √3 / 6
ρ = 1.77 * 1e3  # kg / m^3
m = 1:150 # μg
d = m .* 1e-9 ./ aIRE ./ ρ .* 1e6 # μm
R = d ./ (4 / 3 * ηₕ) # μm

#| output: true
plot(xlabel = "mass (μg)", ylabel = "sample height (μm)")
plot!(m, d, label = "film")
plot!(m, 2 * R, label = "particle (ηₕ)")

Req(m) = m ./ aIRE ./ ρ .* 1e-3 ./ (4 * ηₕ / 3) # μm
m = [10., 20., 50.] # μg
A = hcat(map(m -> prefactor(ν, n₂fixed, Req(m), configfixed) .* α ./ ρ .* m .* 1e-9, m)...)

#| output: true
plot(ν, A,
     label = reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance",
     xflip = true,
     legend = :topleft, legendtitle = "mass deposited")

#| output: true
plot(ν, A,
     label = reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance",
     xflip = true,
     legend = :topleft, legendtitle = "mass deposited")
xlims!(800, 1400)

#| output: true
pl = plot(layout = (2, 1))
plot!(pl[1], xflip = true, ylim = [0, 4],
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "dₚ (μm)", legend = :topleft)
plot!(pl[1], ν, dp(λ, configfixed), 
     label = "thin film")
plot!(pl[1], ν, dp(λ, n₂fixed, ATRConfig(N, aIRE, θ, refidx(1e4 / 2000, :ZnSe))), 
     label = "bulk")
d = range(20, 100, step=20) .* 1e-3 # μm
R = d ./ (4 * ηₕ / 3) # μm
plot!(pl[2], ν, hcat(map(R -> χ(ν, R, configfixed), R)...),
     label = reshape(map(d -> @sprintf("%.0f nm", d * 1e3), d), (1, :)), 
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "χ",
     xflip = true, ylim = [0.5, 1.0],
     legend = :bottomright, legendtitle = "thickness")

#| output: true
m = [10, 20, 50, 100, 120]
plot(ν, hcat(map(m -> χ(ν, Req(m), configfixed), m)...),
     label = reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "χ",
     xflip = true, ylim = [0.5, 1.0],
     legend = :bottomright, legendtitle = "mass deposited")
