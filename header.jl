using Pkg
Pkg.activate(@__DIR__)

using CairoMakie
import ComradeBase
using VLBISkyModels
import VLBISkyModels: paramtype, _getxy

# Abstract Type for Limacon Models
abstract type AbstractLimaconModel <: VLBISkyModels.AbstractModel end
ComradeBase.visanalytic(::Type{<:AbstractLimaconModel}) = VLBISkyModels.NotAnalytic()
ComradeBase.imanalytic(::Type{<:AbstractLimaconModel}) = VLBISkyModels.IsAnalytic()

# Delta limacon model
struct MLimacon{T,V<:Union{AbstractVector{T},NTuple}} <: AbstractLimaconModel
    """
    Real Fourier mode coefficients
    """
    α::V
    """
    Imaginary Fourier mode coefficients
    """
    β::V
    """
    Limacon Shape Parameter
    """
    λ::T
end


"""
    MLimacon(c::Union{NTuple{N, <:Complex}, λ, AbstractVector{<:Complex}})

Construct an MRing geometric model from a complex vector `c`
that correspond to the real and imaginary (or cos and sin) coefficients
of the Fourier expansion. The `N` in the type defines the order of
the Fourier expansion.
"""
function MLimacon(c::Union{AbstractVector{<:Complex},NTuple{N,<:Complex}}, λ::Real) where {N}
    α = real.(c)
    β = imag.(c)
    αT, βT, λT = promote(α, β, λ)
    return MLimacon(αT, βT, λT)
end

function MLimacon(a::Real, b::Real, λ::Real)
    aT, bT, λT = promote(a, b, λ)
    return MLimacon((aT,), (bT,), λT)
end

VLBISkyModels.radialextent(::MLimacon{T,V}) where {T,V} = convert(paramtype(T), 3)

@inline function VLBISkyModels.intensity_point(m::MLimacon{D,V}, p) where {D,V}
    x, y = _getxy(p)
    T = paramtype(D)
    r = hypot(x, y)
    θ = atan(-x, y)
    dr = T(0.02)
    @unpack_params α, β, λ = m(p)
    if (abs(r - 1 - λ * cos(θ)) < dr / 2)
        acc = one(T)
        for n in eachindex(α, β)
            s, c = sincos(n * θ)
            acc += 2 * (α[n] * c - β[n] * s)
        end
        return acc / (2 * T(π) * dr)
    else
        return zero(T)
    end
end

@inline VLBISkyModels.flux(::MLimacon{T,V}) where {T,V} = one(T)
