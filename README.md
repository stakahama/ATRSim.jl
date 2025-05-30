# ATRSim.jl


[![Build
Status](https://github.com/stakahama/ATRSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stakahama/ATRSim.jl/actions/workflows/CI.yml?query=branch%3Amain)

Attenuated total reflection (ATR) simulator with focus on
(semi-)quantitative analysis of deposited thin films using unpolarized
radiation (from Arangio et
al. doi:[10.1177/0003702818821330](https://doi.org/10.1177/0003702818821330),
2019). ATR model equations originally taken from the following sources:

- Harrick, N. J., *Internal Reflection Spectroscopy*, Harrick Scientific
  Corp.: Ossining, NY, 1979.
- Mirabella Jr., F. M., “Internal Reflection Spectroscopy”, *Applied
  Spectroscopy Reviews*, Vol. 21, No. 1-2, p. 45-178, Taylor & Francis,
  doi:[10.1080/05704928508060428](https://doi.org/10.1080/05704928508060428)
  1985.
- Milosevic M., *Internal Reflection and ATR Spectroscopy*, John Wiley &
  Sons: Hoboken, NJ, 2012.

Upon deposition, analytes can form a liquid or amorphous solid film, or
crystalline particles. In the last case, the spectral response may be
modified by a factor (χ) due to increased sample height of particles
(assumed spherical) formed instead of a homogeneous film.

To generate the .jl script from this file, use the following commands.

``` bash
quarto convert README.qmd 
jupyter nbconvert --to script README.ipynb 
mv README.txt README.jl
```

or

``` bash
jupytext --to jl README.qmd
```

## Installation

``` julia
using Pkg
Pkg.add(url="https://github.com/stakahama/SellmeierEqn.jl")
Pkg.add(url="https://github.com/stakahama/ATRSim.jl")
```

Additionally, [Integrals.jl](https://github.com/SciML/Integrals.jl) is a
dependency for calculating the particle modification factor (χ).

Products exported by `ATRSim`.

- `ATRConfig` type
  - `Thinfilm` subtype
  - `Bulk` subtype
- normalized intensity
  - `E₀²`, `E0sq` (non-unicode alias): from Harrick 1979
  - `t²`, `tsq` (non-unicode alias): from Beattie et al. 2000
- `dp`: penetration depth
- `deff`: effective thickness for thick sample
  (`n₂₁ * E₀² / cos(θ) * dp / 2`)
- `ξ`, `xi` (non-unicode alias): `deff / d` for thin sample
  (`n₂₁ * E₀² / cos(θ)`)
- `χ`, `chi` (non-unicode alias): particle factor
- `prefactor`: `N / aIRE * ξ * χ`

Arguments to many of these functions take the form of vectors `ν`, `λ`,
and refractive indices. They are assumed to be aligned a priori.

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

Complex refractive index for ammonium sulfate by Myers et al. *Appl.
Spectrosc.*,
doi:[10.1177/0003702820928358](https://doi.org/10.1177/0003702820928358),
2020. (The data file can be downloaded directly from the
[article](https://journals.sagepub.com/doi/suppl/10.1177/0003702820928358/suppl_file/sj-pdf-1-asp-10.1177_0003702820928358.pdf).)
`n` below is the same as `n₂` (refractive index of sample medium) in the
notation of Harrick et al.

``` julia
pathname = joinpath(@__DIR__, "examples", "Myers2020_ammsulf.csv")
(ν, n, k) = let ammsulf = @subset CSV.read(pathname, DataFrame; comment="#") @. (650 <= :wavenum <= 4000)
    @with ammsulf (:wavenum, :n, :k)
end
n₂fixed = 1.5 # alternate n₂
```

Decadic linear absorption coefficient.

``` julia
α = @. 1 / log(10) * 4 * π * ν * k * 1e2 # m^-1
```

Define some basic parameters for the HATR (Horizontal Attenuated Total
Reflectance) accessory with ZnSe IRE (infrared element).

``` julia
λ = 1e4 ./ ν # μm
n₁ = refidx(λ, :ZnSe)
N = 10
aIRE = 8.0 * 1e-4 # m^2
θ = 45 / 180 * π
```

The optional last argument defines the refractive index of the third
medium above the sample (air) when the sample in contact with the IRE is
prepared as a thin film. The data type of variable returned determines
the behavior of some functions (e.g., penetration depth, electric field
intensity at the interface).

``` julia
nₐᵢᵣ = 1.0
n₁fixed = refidx(1e4 / 2000, :ZnSe)
configfilm = ATRConfig(N, aIRE, θ, n₁, nₐᵢᵣ)
configfixed = ATRConfig(N, aIRE, θ, n₁fixed, nₐᵢᵣ)
configbulk = ATRConfig(N, aIRE, θ, n₁)
```

``` julia
println([typeof(configfilm), typeof(configfixed), typeof(configbulk)])
```

    DataType[Thinfilm, Thinfilm, Bulk]

Assume a thickness and particle size for the following examples.

``` julia
d = 0.05 # μm
R = 0.08 # μm
```

Comparison of particle spectra with different assumptions for n₂
(refractive index of sample) - either varying with ammonium sulfate, or
fixed at 1.5.

``` julia
plot(xflip=true, legend=:topleft, xlabel="Wavenumber (cm⁻¹)", ylabel="Absorbance")
plot!(ν, α .* 1e-6 .* d .* ξ(n, configfilm) .* χ(ν, n, R, configfilm), label="varying n₂")
plot!(ν, α .* 1e-6 .* d .* ξ(n₂fixed, configfixed) .* χ(ν, n₂fixed, R, configfixed), label="fixed n₂")
plot!(size=(400, 250))
```

![](README_files/figure-commonmark/cell-9-output-1.svg)

Proceeding with the fixed-parameter case, compare spectra with different
ATR models.

``` julia
plot(xflip=true, legend=:topleft, xlabel="Wavenumber (cm⁻¹)", ylabel="Absorbance")
plot!(ν, α .* 1e-6 .* d .* ξ(n₂fixed, configfixed), label="thin film")
plot!(ν, α .* 1e-6 .* d .* ξ(n₂fixed, configfixed) .* χ(ν, n, R, configfixed), label="particle")
plot!(ν, 0.05 .* α .* 1e-6 .* deff(ν, n₂fixed, configbulk), label="0.05 × bulk")
plot!(size=(400, 250))
```

![](README_files/figure-commonmark/cell-10-output-1.svg)

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
plot(xlabel = "mass (μg)", ylabel="sample height (μm)")
plot!(m, d, label="film")
plot!(m, 2 * R, label="particle (at ηₕ)")
plot!(size=(400, 250))
```

![](README_files/figure-commonmark/cell-12-output-1.svg)

Simulation of absorbance.

``` julia
Req(m) = m ./ aIRE ./ ρ .* 1e-3 ./ (4 * ηₕ / 3) # μm
```

``` julia
m = [10., 20., 50.] # μg
A = stack(map(m -> prefactor(ν, n₂fixed, Req(m), configfixed) .* α ./ ρ .* m .* 1e-9, m))
```

``` julia
plot(ν, A,
     label=reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel="Wavenumber (cm⁻¹)", ylabel="Absorbance",
     xflip=true,
     legend=:topleft, legendtitle="mass deposited",
     legendtitlefontsize=8)
plot!(size=(400, 250))
```

![](README_files/figure-commonmark/cell-15-output-1.svg)

``` julia
plot(ν, A,
     label=reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     xlabel="Wavenumber (cm⁻¹)", ylabel="Absorbance",
     xflip=true,
     legend=:topleft, legendtitle="mass deposited",
     legendtitlefontsize=8)
xlims!(800, 1400)
plot!(size=(400, 250))
```

![](README_files/figure-commonmark/cell-16-output-1.svg)

The following figure reproduces the magnitude of particle modification
factor from Figure S2 (bottom panel) of Arangio et al.

``` julia
pl = plot(layout=(2, 1), xflip=true, xlim=[500, 4000], 
      xlabel="Wavenumber (cm⁻¹)")
plot!(pl[1], ν, dp(λ, n₂fixed, configfixed), label=false,
      ylim=[0, 4], ylabel="dₚ (μm)")
d = range(20, 100, step=20) .* 1e-3 # μm
R = d ./ (4 * ηₕ / 3) # μm
plot!(pl[2], ν, stack(map(R -> χ(ν, n₂fixed, R, configfixed), R)),
      label=reshape(map(d -> @sprintf("%.0f nm", d * 1e3), d), (1, :)), 
      line_z=d', color=cgrad(:blues, rev=true), colorbar=false,
      ylabel="χ", ylim=[0.7, 1.0],
      yticks=0.7:0.05:1.0,
      legend=:bottomright, legendtitle="thickness", legendtitlefontsize=8)
plot!(size=(400, 500))
```

![](README_files/figure-commonmark/cell-17-output-1.svg)

Replot as a function of deposited mass assuming hexagonal packing, which
provides a lower bound on the sample height above IRE in particle form.

``` julia
m = [10, 20, 50, 100, 120]
plot(ν, stack(map(m -> χ(ν, n₂fixed, Req(m), configfixed), m)),
     label=reshape(map(m -> @sprintf("%.0f μg", m), m), (1, :)),
     line_z=m', color=cgrad(:blues, rev=true), colorbar=false,
     xlabel="Wavenumber (cm⁻¹)", ylabel="χ",
     xflip=true, ylim=[0.7, 1.0],
     legend=:bottomright, legendtitle="mass deposited",
     legendtitlefontsize=8)
plot!(size=(400, 250))
```

![](README_files/figure-commonmark/cell-18-output-1.svg)
