# ATRSim.jl


[![Build
Status](https://github.com/stakahama/ATRSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stakahama/ATRSim.jl/actions/workflows/CI.yml?query=branch%3Amain)

ATR calculations from Arangio et
al. doi:[10.1177/0003702818821330](https://doi.org/10.1177/0003702818821330),
2019.

Following code is found in `examples/example.jl`.

## Installation

``` julia
using Pkg
Pkg.add(url = "https://github.com/stakahama/SellmeierEqn.jl")
Pkg.add(url = "https://github.com/stakahama/ATRSim.jl")
```

## Usage

``` julia
using DataFrames
using DataFramesMeta
using CSV
using Plots
using Printf
using SellmeierEqn
using ATRSim
```

``` julia
pathname = joinpath(@__DIR__, "examples", "Myers2020_ammsulf.csv")
(ν, n, k) = let ammsulf = @subset CSV.read(pathname, DataFrame; comment="#") @. (650 <= :wavenum <= 4000)
    @with ammsulf (:wavenum, :n, :k)
end
```

``` julia
α = @. 1 / log(10) * 4 * π * ν * k * 1e2 # m^-1
```

``` julia
λ = 1e4 ./ ν # μm
n₁ = refidx(λ, :ZnSe)
N = 10
aIRE = 8.0 * 1e-4 # m^2
θ = 45 / 180 * π
```

``` julia
configfilm = ATRConfig(N, aIRE, θ, n₁, 1.0)
configfixed = ATRConfig(N, aIRE, θ, refidx(1e4 /2000, :ZnSe), 1.0)
configbulk = ATRConfig(N, aIRE, θ, n₁)
```

``` julia
d = 0.05 # μm
R = 0.08 # μm
n₂fixed = 1.5
```

Comparison of particle spectra with different assumptions for n₂
(refractive index of sample).

``` julia
plot(xflip = true, legend = :topleft, xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance")
plot!(ν, α .* d .* ξ(n, configfilm) .* χ(ν, R, configfilm), label = "varying n₂")
plot!(ν, α .* d .* ξ(n₂fixed, configfixed) .* χ(ν, R, configfixed), label = "fixed n₂")
```

![](README_files/figure-commonmark/cell-8-output-1.svg)

Comparison of spectra with different ATR models.

``` julia
plot(xflip = true, legend = :topleft, xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance")
plot!(ν, α .* d .* ξ(n₂fixed, configfixed), label = "thin film")
plot!(ν, α .* d .* ξ(n₂fixed, configfixed) .* χ(ν, R, configfixed), label = "particle")
plot!(ν, 0.1 .* α .* deff(ν, n₂fixed, configbulk), label = "0.1 × bulk")
```

![](README_files/figure-commonmark/cell-9-output-1.svg)

Relation among deposited mass, film thickness, and equivalent particle
radius (for hexagonal circle packing assumption).

``` julia
ηₕ = π * √3 / 6
ρ = 1.77 * 1e3  # kg / m^3
m = 1:150 # μg
d = m .* 1e-9 ./ aIRE ./ ρ .* 1e6 # μm
R = d ./ (4 / 3 * ηₕ) # μm
```

``` julia
plot(xlabel = "mass (μg)", ylabel = "sample height (μm)")
plot!(m, d, label = "film")
plot!(m, 2 * R, label = "particle (ηₕ)")
```

![](README_files/figure-commonmark/cell-11-output-1.svg)

Simulation of absorbance.

``` julia
Req(m) = m ./ aIRE ./ ρ .* 1e-3 ./ (4 * ηₕ / 3) # μm
m = [10., 20., 50.] # μg
A = hcat(map(m -> prefactor(ν, n₂fixed, Req(m), configfixed) .* α ./ ρ .* m .* 1e-9, m)...)
```

``` julia
plot(ν, A,
     label = reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance",
     xflip = true,
     legend = :topleft, legendtitle = "mass deposited")
```

![](README_files/figure-commonmark/cell-13-output-1.svg)

``` julia
plot(ν, A,
     label = reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "Absorbance",
     xflip = true,
     legend = :topleft, legendtitle = "mass deposited")
xlims!(800, 1400)
```

![](README_files/figure-commonmark/cell-14-output-1.svg)

The following figure is slightly different from Figure S2 (bottom panel)
as the penetration depth below is assumed to not be influenced by the
sample medium n₂ (bulk), but instead n₃ (thin flim).

``` julia
pl = plot(layout = (2, 1))
plot!(pl[1], xflip = true, ylim = [0, 4],
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "dₚ (μm)", legend = :topleft)
plot!(pl[1], ν, dp(λ, configfixed), 
     label = "thinfilm")
plot!(pl[1], ν, dp(λ, n₂fixed, ATRConfig(N, aIRE, θ, refidx(1e4 / 2000, :ZnSe))), 
     label = "bulk")
d = range(20, 100, step=20) .* 1e-3 # μm
R = d ./ (4 * ηₕ / 3) # μm
plot!(pl[2], ν, hcat(map(R -> χ(ν, R, configfixed), R)...),
     label = reshape(map(d -> @sprintf("%.0f nm", d * 1e3), d), (1, :)), 
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "χ",
     xflip = true, ylim = [0.5, 1.0],
     legend = :bottomright, legendtitle = "thickness")
```

![](README_files/figure-commonmark/cell-15-output-1.svg)

``` julia
m = [10, 20, 50, 100, 120]
plot(ν, hcat(map(m -> χ(ν, Req(m), configfixed), m)...),
     label = reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel = "Wavenumber (cm⁻¹)", ylabel = "χ",
     xflip = true, ylim = [0.5, 1.0],
     legend = :bottomright, legendtitle = "mass deposited")
```

![](README_files/figure-commonmark/cell-16-output-1.svg)
