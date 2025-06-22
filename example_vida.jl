include(joinpath(@__DIR__, "header.jl"))

using CairoMakie
using OptimizationBBO
using VIDA
import VIDA: Constant

#
# Limacon Parameters
#
truth = (
    # Limacon Shape
    λ1=μas2rad(20),
    λ2=0.5,
    # rotation
    PA=π / 2,
    # blurring kernel FWHM
    σinner=μas2rad(1) / √(8 * log(2)),
    σouter=μas2rad(1.5) / √(8 * log(2)),
)

#Define our bounds
lower = (
    # Limacon Shape
    λ1=μas2rad(0.01),
    λ2=0,
    # rotation
    PA=-π,
    # blurring kernel FWHM
    σinner=μas2rad(0.01) / √(8 * log(2)),
    σouter=μas2rad(0.01) / √(8 * log(2)),
)

upper = (
    # Limacon Shape
    λ1=μas2rad(100),
    λ2=1,
    # rotation
    PA=1 * π,
    # blurring kernel FWHM
    σinner=μas2rad(20) / √(8 * log(2)),
    σouter=μas2rad(20) / √(8 * log(2)),
)

# Template Model
template(θ) = DblGaussianLimacon(θ.λ1, θ.λ2, θ.σinner, θ.σouter, θ.PA, 0, 0)

# image fov and number of pixels
fov = μas2rad(100)
nx = 128

# --------
# create a groundtruth image
grid = imagepixels(fov, fov, nx, nx)
# create a limacon model on the groundturth parameters
gtmodel = template(truth)
# map out the truth limacon model to the image grid
imgt = intensitymap(gtmodel, grid)
# plot
ax = imageviz(imgt)
save("example_vida_groundtruth.png", ax)

# -------
# Build the divergence we want to fit
bh = Bhattacharyya(imgt)

# Make an analytic image model

# Define the fitting problem
prob = VIDAProblem(
    bh,
    template,
    lower,
    upper)

using OptimizationBBO
xopt, opt_temp, divmin = vida(prob, BBO_adaptive_de_rand_1_bin(); maxiters=50_000)

@info "Fitting Results (left: groundtruth, right: fitted values)"
@info "λ1: $(rad2μas(truth.λ1)) μas, $(rad2μas(xopt.λ1)) μas"
@info "λ2: $(truth.λ2), $(xopt.λ2)"
@info "PA: $(rad2deg(truth.PA)), $(rad2deg(xopt.PA))"
@info "σinner: $(rad2μas(truth.σinner)) μas, $(rad2μas(xopt.σinner)) μas"
@info "σouter: $(rad2μas(truth.σouter)) μas, $(rad2μas(xopt.σouter)) μas"

#plot the results
fig = triptic(imgt, opt_temp)
save("example_vida_fitted.png", fig)