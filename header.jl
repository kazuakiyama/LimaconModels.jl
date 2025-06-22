using Pkg
Pkg.activate(@__DIR__)

import ComradeBase
using DocStringExtensions
using VLBISkyModels
import VLBISkyModels: AbstractImageTemplate,
    AbstractAzimuthal, AbstractRadial,
    AzimuthalUniform, AzimuthalCosine,
    paramtype, _getxy,
    ringphase, radial_profile, azimuthal_profile,
    intensity_point

# Define the radial profile of the curve
abstract type AbstractPolarCurve end

"""
    $(TYPEDEF)

A standard Limacon curve defined by r(ϕ)=1+λcosϕ

## Fields
$(FIELDS)
"""
struct LimaconCurve{T} <: AbstractPolarCurve
    """ Limacon Shape Parameter """
    λ::T
end


@inline function radial_point(m::LimaconCurve, p)
    @unpack_params λ = m(p)
    ϕ = p.ϕ
    return 1 + λ * cos(ϕ)
end


"""
    $(TYPEDEF)

A flexible limacon template that forms a limacon by taking the product of
a radial and azimuthal brightness profile.

A list of radial profiles is given by `subtypes(AbstractRadial)`

A list of azimuthal profiles is given by `subtypes(AbstractAzimuthal)`

## Examples
```julia-repl
julia> rad = RadialGaussian(0.1)
julia> azi = AzimuthalUniform()
julia> limacon = modify(LimaconTemplate(rad, azi), Stretch(10.0), Shift(1.0, 2.0))
```

## Fields
$(FIELDS)
"""
struct PolarTemplate{R<:AbstractRadial,A<:AbstractAzimuthal,P<:AbstractPolarCurve} <: AbstractImageTemplate
    """
    Radial profile of the limacon
    """
    radial::R
    """
    Azimuthal profile of the limacon
    """
    azimuthal::A
    """
    Limacon Shape Parameter
    """
    polarcurve::P
end

@inline function intensity_point(d::PolarTemplate, p)
    (; X, Y) = p
    @unpack_params radial, azimuthal, polarcurve = d(p)
    r = hypot(X, Y)
    ϕ = ringphase(X, Y)
    # radial curve
    pϕ = merge((; ϕ), p)
    rϕ = radial_point(polarcurve, pϕ)
    # radial profile
    pr = merge((; r, rϕ), p)
    fr = radial_profile(radial, pr)
    # azimuthal profile
    fϕ = azimuthal_profile(azimuthal, pϕ)
    return fr * fϕ
end


"""
    PolarGaussian(σ)

Create a Gaussian profile with a standard deviation `σ`.

## Notes
This is usually couple with a azimuthal profile to create a general ring template

```julia-repl
julia> rad = PolarGaussian(0.1)
julia> azi = AzimuthalUniform()
julia> t = RingTemplate(rad, azi)
```

## Arguments
  - `σ`: The standard deviation for the Gaussian ring.
"""
struct PolarGaussian{T} <: AbstractRadial
    σ::T
end

@inline @fastmath function radial_profile(d::PolarGaussian, p)
    @unpack_params σ = d(p)
    r = p.r
    rϕ = p.rϕ
    return exp(-(r - rϕ)^2 / (2 * σ^2))
end

radialextent(d::PolarGaussian) = 2 + 3 * d.σ


"""
    PolarDblGaussian(σ)

Create an asymmetric Gaussian profile with two different standard deviations inside and outside of the peak.

## Notes
This is usually couple with a azimuthal profile to create a general ring template

```julia-repl
julia> rad = RadialGaussian(0.1)
julia> azi = AzimuthalUniform()
julia> t = RingTemplate(rad, azi)
```

## Arguments
  - `σinner`: The standard deviation for the Gaussian ring for r < P(ϕ)
  - `σouter`: The standard deviation for the Gaussian ring for r > P(ϕ)

"""
struct PolarDblGaussian{T} <: AbstractRadial
    σinner::T
    σouter::T
end

@inline @fastmath function radial_profile(d::PolarDblGaussian, p)
    @unpack_params σinner, σouter = d(p)
    r = p.r
    rϕ = p.rϕ
    if r < rϕ
        return exp(-(r - rϕ)^2 / (2 * σinner^2))
    else
        return exp(-(r - rϕ)^2 / (2 * σouter^2))
    end
end

radialextent(d::PolarDblGaussian) = 2 + 3 * d.σouter


"""
    $(SIGNATURES)

Implements the gaussian limacon template.

## Notes
This is a convienence constructor that uses [`PolarTemplate`](@ref) under
the hood. To create this function your self do

```julia
PolarTemplate(
    PolarGaussian(σ),
    AzimuthalUniform(),
    LimaconCurver(λ)
)
```

## Arguments
- `λ` : Limacon shape parameter. The shape is given by r(ϕ) = 1 + λ*cosϕ
- `σ` : standard deviation of the Gaussian Limacone
"""
@inline GaussianLimacon(λ, σ) = PolarTemplate(
    PolarGaussian(σ),
    AzimuthalUniform(),
    LimaconCurve(λ)
)


"""
    $(SIGNATURES)

Implements the gaussian limacon template.

## Notes
This is a convienence constructor that uses [`PolarTemplate`](@ref) under
the hood. To create this function your self do

```julia
modify(
    GaussianLimacon(λ2, σ / λ1),
    Stretch(λ1),
    Shift(x0, y0)
)
```

## Arguments
- `λ1, λ2` : Limacon shape parameters. The shape is given by r(ϕ) = λ1(1 + λ2*cosϕ)
- `σ` : standard deviation of the Gaussian Limacone
- `ϕ` : orientation of the limacon
- `x0`: location of the ring center horizontally
- `y0`: location of the ring center vertically
"""
@inline GaussianLimacon(λ1, λ2, σ, ϕ, x0, y0) = modify(
    GaussianLimacon(λ2, σ / λ1),
    Stretch(λ1),
    Rotate(ϕ),
    Shift(x0, y0)
)


"""
    $(SIGNATURES)

Implements the double gaussian limacon template.

## Notes
This is a convienence constructor that uses [`PolarTemplate`](@ref) under
the hood. To create this function your self do

```julia
PolarTemplate(
    PolarDblGaussian(σinner, σouter),
    AzimuthalUniform(),
    LimaconCurve(λ)
)
```

## Arguments
- `λ` : Limacon shape parameter. The shape is given by r(ϕ) = 1 + λcosϕ
- `σinner, σouter` : standard deviations of the Gaussian Limacone inner/outer the peak
"""
@inline DblGaussianLimacon(λ, σinner, σouter) = PolarTemplate(
    PolarDblGaussian(σinner, σouter),
    AzimuthalUniform(),
    LimaconCurve(λ)
)


"""
    $(SIGNATURES)

Implements the double gaussian limacon template.

## Notes
This is a convienence constructor that uses [`PolarTemplate`](@ref) under
the hood. To create this function your self do

```julia
PolarTemplate(
    PolarDblGaussian(λ, σinner, σouter),
    AzimuthalUniform(),
    LimaconCurver(λ)
)
```

## Arguments
- `λ1, λ2` : Limacon shape parameters. The shape is given by λ1(1 + λ2*cosϕ)
- `σinner, σouter` : standard deviations of the Gaussian Limacone inner/outer the peak
- `ϕ` : orientation of the limacon
- `x0`: location of the Limacon center horizontally
- `y0`: location of the Limacon center vertically
"""
@inline DblGaussianLimacon(λ1, λ2, σinner, σouter, ϕ, x0, y0) = modify(
    DblGaussianLimacon(λ2, σinner / λ1, σouter / λ1),
    Stretch(λ1),
    Rotate(ϕ),
    Shift(x0, y0)
)