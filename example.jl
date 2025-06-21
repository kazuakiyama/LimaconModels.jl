include(joinpath(@__DIR__, "header.jl"))

using Plots

#
# Limacon Parameters
#

# limacon shape
λ1 = μas2rad(20)
λ = 0.4 # λ2/λ1
# brightness profile
α = [0.2]
β = [0.0]
# rotation
PA = -π / 2
# blurring kernel FWHM
fwhm = μas2rad(1) / √(8 * log(2))
# image fov and number of pixels
fov = μas2rad(100)
nx = 1000

#
# Create a m-Limacon model
#

# create the image grid
grid = imagepixels(fov, fov, nx, nx)

# create bluring kernel
g = stretched(Gaussian(), fwhm, fwhm)

# create a m-limacon model
#.  resize Limacon, rotate by -90 deg, and then convolve
mlmodel = convolved(rotated(stretched(MLimacon(α, β, λ), λ1, λ1), PA), g)

# map out the m-limacon model to the image grid
iml = intensitymap(mlmodel, grid)

# plot
ax = imageviz(iml)
save("example_ring.png", ax)