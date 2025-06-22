include(joinpath(@__DIR__, "header.jl"))

using CairoMakie

#
# Limacon Parameters
#

# limacon shape
λ1 = μas2rad(20)
λ2 = 0.5
# brightness profile
α = [0.2]
β = [0.0]
# rotation
PA = π / 2
# blurring kernel FWHM
σinner = μas2rad(1) / √(8 * log(2))
σouter = μas2rad(1.5) / √(8 * log(2))
# image fov and number of pixels
fov = μas2rad(100)
nx = 1000

#
# Create a m-Limacon model
#

# create the image grid
grid = imagepixels(fov, fov, nx, nx)

# create a m-limacon model
#.  resize Limacon, rotate by -90 deg, and then convolve
mlmodel = DblGaussianLimacon(λ1, λ2, σinner, σouter, PA, 0, 0)

# map out the m-limacon model to the image grid
iml = intensitymap(mlmodel, grid)

# plot
ax = imageviz(iml)
save("example_ring.png", ax)