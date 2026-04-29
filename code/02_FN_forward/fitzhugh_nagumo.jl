using DifferentialEquations
using Plots
using Statistics
using Random

ENV["GKSwstype"] = "100"
Plots.default(show=false)

# FitzHugh-Nagumo model
#
#   dV/dt = c * (V - V³/3 - R + I_exp)
#   dR/dt = (V + a - b·R) / c
#
# Variables:
#   V      : membrane potential (activator)
#   R      : recovery variable (inhibitor)
#
# Parameters:
#   a      : offset in recovery nullcline
#   b      : sensitivity of recovery to potential
#   c      : time scale ratio (couples fast/slow dynamics)
#   I_exp  : external input current

function fitzhugh_nagumo!(du, u, p, t)
    V, R = u
    a, b, c, I_exp = p
    du[1] = c * (V - (V^3) / 3 - R + I_exp)
    du[2] = (V - a - b * R) / c
end

# Parameters
a     = 0.2
b     = 0.3
c     = 3.0
I_exp = 0.0
p     = [a, b, c, I_exp]

# Initial conditions
u0    = [-1.0, 1.0]
tspan = (0.0, 20.0)

# Solver: use Euler() with a coarse dt to introduce numerical error and wiggles
# in the loss landscape. Decrease dt for a more accurate (smoother) solution.
dt = 0.01

# Solve
prob = ODEProblem(fitzhugh_nagumo!, u0, tspan, p)
sol  = solve(prob, Euler(), dt=dt, saveat=0.1)

# Ground truth with Gaussian noise
σ_noise = 0.0
rng     = MersenneTwister(42)
gt      = Array(sol) .+ σ_noise .* randn(rng, size(Array(sol)))

# Time series
p1 = plot(sol.t, sol[1, :], label="V (potencial)", xlabel="Tiempo", ylabel="Amplitud",
          title="FitzHugh-Nagumo: Series de tiempo", linewidth=2, color=:blue)
plot!(p1, sol.t, sol[2, :], label="R (recuperación)", linewidth=2, color=:red)
scatter!(p1, sol.t, gt[1, :], label="V (obs.)", markersize=2, color=:blue, alpha=0.3)
scatter!(p1, sol.t, gt[2, :], label="R (obs.)", markersize=2, color=:red, alpha=0.3)

# Phase portrait
p2 = plot(sol[1, :], sol[2, :], label="Trayectoria", xlabel="V", ylabel="R",
          title="Retrato de fase", linewidth=2, color=:purple)
scatter!(p2, [u0[1]], [u0[2]], label="Condición inicial", markersize=6, color=:green)

fig = plot(p1, p2, layout=(1, 2), size=(900, 400), bottom_margin=10Plots.mm)
savefig(fig, "fitzhugh_nagumo.png")
println("Figura guardada en fitzhugh_nagumo.png")

# ---------------------------------------------------------------------------
# Loss landscape: MSE vs ground truth over 2D grids of parameter pairs
# ---------------------------------------------------------------------------

p_true = [a, b, c, I_exp]

function mse_loss_pair(i, j, vi, vj)
    p_test    = copy(p_true)
    p_test[i] = vi
    p_test[j] = vj
    prob_test = ODEProblem(fitzhugh_nagumo!, u0, tspan, p_test)
    sol_test  = solve(prob_test, Euler(), dt=dt, saveat=0.1)
    length(sol_test.t) != length(sol.t) && return Inf
    return mean((Array(sol_test) .- gt) .^ 2)
end

# Parameter metadata: (index, symbol, label, range)
param_info = [
    (1, "a",     "a (nullcline)",     range(-1.0, 1.0,  length=300)),
    (2, "b",     "b (sensibilidad)",  range(-1.0, 1.0,  length=300)),
    (3, "c",     "c (escala tiempo)", range(1.0,  6.0,  length=300)),
    (4, "I_exp", "I_exp (corriente)", range(-1.0, 1.0,  length=300)),
]

pairs = [(1,2)]

for (i, j) in pairs
    idx_i, sym_i, lbl_i, range_i = param_info[i]
    idx_j, sym_j, lbl_j, range_j = param_info[j]

    local loss_grid    = [mse_loss_pair(idx_i, idx_j, vi, vj)
                          for vj in range_j, vi in range_i]
    local loss_max    = quantile(filter(isfinite, vec(loss_grid)), 0.95)
    local loss_clipped = min.(loss_grid, loss_max)

    fig_loss = heatmap(range_i, range_j, loss_clipped,
                       xlabel=lbl_i, ylabel=lbl_j,
                       title="MSE: $(sym_i) vs $(sym_j)",
                       color=cgrad([:white, :royalblue]),
                       colorbar_title="MSE", clims=(0, loss_max))

    contour!(fig_loss, range_i, range_j, loss_clipped,
             levels=10, color=:black, linewidth=0.8, alpha=0.5, label="")

    scatter!(fig_loss, [p_true[idx_i]], [p_true[idx_j]],
             markersize=8, markershape=:star5, color=:cyan,
             label="Parámetros verdaderos")

    # For the a vs b pair, add 1D cross-sections through the true parameter values
    if (i, j) == (1, 2)
        # Slice along a, fixing b = b_true
        b_true_idx = argmin(abs.(collect(range_j) .- p_true[idx_j]))
        slice_i = loss_clipped[b_true_idx, :]
        p_slice_i = plot(range_i, slice_i, xlabel=lbl_i, ylabel="MSE",
                         title="Corte: MSE vs $(sym_i)  ($(sym_j) = $(p_true[idx_j]))",
                         linewidth=2, color=:royalblue, label="")
        vline!(p_slice_i, [p_true[idx_i]], color=:cyan, linestyle=:dash, label="verdadero")

        # Slice along b, fixing a = a_true
        a_true_idx = argmin(abs.(collect(range_i) .- p_true[idx_i]))
        slice_j = loss_clipped[:, a_true_idx]
        p_slice_j = plot(range_j, slice_j, xlabel=lbl_j, ylabel="MSE",
                         title="Corte: MSE vs $(sym_j)  ($(sym_i) = $(p_true[idx_i]))",
                         linewidth=2, color=:royalblue, label="")
        vline!(p_slice_j, [p_true[idx_j]], color=:cyan, linestyle=:dash, label="verdadero")

        fig_combined = plot(fig_loss, p_slice_i, p_slice_j,
                            layout=@layout([a{0.6h}; b c]),
                            size=(1200, 1200), dpi=150, bottom_margin=10Plots.mm)
        savefig(fig_combined, "loss_$(sym_i)_$(sym_j).png")
        println("Guardado: loss_$(sym_i)_$(sym_j).png")
    else
        fname = "loss_$(sym_i)_$(sym_j).png"
        savefig(fig_loss, fname)
        println("Guardado: $fname")
    end
end
