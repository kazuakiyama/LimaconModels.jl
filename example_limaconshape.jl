include(joinpath(@__DIR__, "header.jl"))

using CairoMakie

# limacon shape
λ = 0.5

# Define the Limacon curve
limacon = LimaconCurve(λ)

# Define ϕ and calculate the radial points
ϕ = range(0, 2π, length=1000)
r = map(ϕ -> radial_point(limacon, (ϕ=ϕ,)), ϕ)

# Create the plot with Cairo Makie
fig = Figure()
ax = Axis(fig[1, 1], title="Limacon Curve", xlabel="x", ylabel="y")
lines!(ax, r .* cos.(ϕ), r .* sin.(ϕ), color=:blue, linewidth=2)
ax.aspect = DataAspect()
ax.xlabel = "x"
ax.ylabel = "y"
ax.title = "Limacon Curve with λ = $λ"
save("example_limaconshape.png", fig)