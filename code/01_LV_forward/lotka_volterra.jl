using DifferentialEquations
using Plots
using Statistics
using Random

ENV["GKSwstype"] = "100"
Plots.default(show=false)

# Lotka-Volterra equations (predator-prey model)
#
#   dx/dt = α x - β x y    (prey)
#   dy/dt = δ x y - γ y    (predator)
#
# Parameters:
#   α  : prey birth rate
#   β  : predation rate
#   δ  : predator reproduction rate per prey eaten
#   γ  : predator death rate

function lotka_volterra!(du, u, p, t)
    x, y = u          # x = prey, y = predator
    α, β, δ, γ = p
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
end

# Parameters
α = 1.0   # prey birth rate
β = 0.1   # predation rate
δ = 0.075 # predator reproduction rate
γ = 1.5   # predator death rate
p = [α, β, δ, γ]

# Initial conditions: x₀ = 10 prey, y₀ = 5 predators
u0 = [10.0, 5.0]

# Time span
tspan = (0.0, 30.0)

# Define and solve the ODE problem
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.1)

# Ground truth timeseries with additive independent Gaussian noise
# Set σ_noise = 0.0 for a clean ground truth
σ_noise = 5.0
rng = MersenneTwister(42)
gt = Array(sol) .+ σ_noise .* randn(rng, size(Array(sol)))

# Plot time series (clean + noisy observations)
p1 = plot(sol.t, sol[1, :], label="Presas (x)", xlabel="Tiempo", ylabel="Población",
          title="Lotka-Volterra: Series de tiempo", linewidth=2, color=:blue)
plot!(p1, sol.t, sol[2, :], label="Depredadores (y)", linewidth=2, color=:red)
scatter!(p1, sol.t, gt[1, :], label="Presas (obs.)", markersize=2, color=:blue, alpha=0.4)
scatter!(p1, sol.t, gt[2, :], label="Depredadores (obs.)", markersize=2, color=:red, alpha=0.4)

# Phase portrait
p2 = plot(sol[1, :], sol[2, :], label="Trayectoria", xlabel="Presas (x)",
          ylabel="Depredadores (y)", title="Retrato de fase", linewidth=2, color=:purple)
scatter!(p2, [u0[1]], [u0[2]], label="Condición inicial", markersize=6, color=:green)

# Combined figure
fig = plot(p1, p2, layout=(1, 2), size=(900, 400), bottom_margin=10Plots.mm)
savefig(fig, "lotka_volterra.png")
println("Figura guardada en lotka_volterra.png")

# ---------------------------------------------------------------------------
# Loss landscape: MSE vs ground truth over 2D grids of parameter pairs
#
# For each pair, the other two parameters are fixed at their true values.
# The MSE minimum should lie at the true parameter values.
# ---------------------------------------------------------------------------

# True parameter vector: [α, β, δ, γ]
p_true = [α, β, δ, γ]

function mse_loss_pair(i, j, vi, vj)
    p_test = copy(p_true)
    p_test[i] = vi
    p_test[j] = vj
    prob_test = ODEProblem(lotka_volterra!, u0, tspan, p_test)
    sol_test = solve(prob_test, Tsit5(), saveat=0.1)
    length(sol_test.t) != length(sol.t) && return Inf
    return mean((Array(sol_test) .- gt) .^ 2)
end

# Parameter metadata: (index, symbol, label, range)
param_info = [
    (1, "α", "α (nac. presas)",   range(0.4,  1.6, length=80)),
    (2, "β", "β (depredación)",   range(0.04, 0.16, length=80)),
    (3, "δ", "δ (repr. depr.)",   range(0.03, 0.12, length=80)),
    (4, "γ", "γ (muerte depr.)",  range(0.8,  2.2,  length=80)),
]

# All 6 pairs
pairs = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]

for (i, j) in pairs
    idx_i, sym_i, lbl_i, range_i = param_info[i]
    idx_j, sym_j, lbl_j, range_j = param_info[j]

    local loss_grid = [mse_loss_pair(idx_i, idx_j, vi, vj)
                       for vj in range_j, vi in range_i]

    local loss_clipped = min.(loss_grid, 200.0)

    fig_pair = heatmap(range_i, range_j, loss_clipped,
                       xlabel=lbl_i, ylabel=lbl_j,
                       title="MSE: $(sym_i) vs $(sym_j)",
                       color=cgrad([:white, :royalblue]), colorbar_title="MSE", clims=(0, 200))

    contour!(fig_pair, range_i, range_j, loss_clipped,
             levels=10, color=:black, linewidth=0.8, alpha=0.5, label="")

    scatter!(fig_pair, [p_true[idx_i]], [p_true[idx_j]],
             markersize=8, markershape=:star5, color=:cyan,
             label="Parámetros verdaderos")

    fname = "loss_$(sym_i)_$(sym_j).png"
    savefig(fig_pair, fname)
    println("Guardado: $fname")
end
