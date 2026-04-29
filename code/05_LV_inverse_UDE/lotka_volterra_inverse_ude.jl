using ComponentArrays
using DifferentialEquations
using LineSearches
using Lux
using Optimization, OptimizationOptimisers, OptimizationOptimJL
using Plots
using SciMLSensitivity
using Statistics
using Zygote
using Random

Random.seed!(42)

ENV["GKSwstype"] = "100"
Plots.default(show=false)

# ---------------------------------------------------------------------------
# Note: this example is based on the SciML tutorial "Automatically Discover
# Missing Physics by Embedding Machine Learning into Differential Equations":
# https://docs.sciml.ai/Overview/dev/showcase/missing_physics/
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# True Lotka-Volterra system (same parameters as 01_LV_forward)
#
#   dx/dt = α x - β x y    (prey)
#   dy/dt = δ x y - γ y    (predator)
#
# Parameters:
#   α = 1.0   : prey birth rate
#   β = 0.1   : predation rate
#   δ = 0.075 : predator reproduction rate per prey eaten
#   γ = 1.5   : predator death rate
# ---------------------------------------------------------------------------

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
end

p_true = [1.3, 0.9, 0.8, 1.8]   # [α, β, δ, γ]
u0     = [3.0, 1.0]
tspan  = (0.0, 5.0)
t_obs  = range(0.0, tspan[2], step=0.2)

prob_true = ODEProblem(lotka_volterra!, u0, tspan, p_true)
sol_true  = solve(prob_true, Tsit5(), saveat=t_obs)
target    = Array(sol_true)   # clean trajectories used as training target

# ---------------------------------------------------------------------------
# Universal Differential Equation (UDE)
#
# Same structure as 03_LV_forward_UDE: the interaction terms are replaced by
# a small neural network nn(u):
#
#   dx/dt = α x + nn(u)[1]
#   dy/dt = -γ y + nn(u)[2]
#
# In the true system:  nn(u) = [-β x y,  δ x y]
# ---------------------------------------------------------------------------

α = p_true[1]
γ = p_true[4]

# rbf(x) = exp.(-(x .^ 2))
# f_activation = rbf
f_activation = tanh

nn = Chain(
    Dense(2 => 5, f_activation),
    Dense(5 => 5, f_activation),
    Dense(5 => 5, f_activation),
    Dense(5 => 2),
    # WrappedFunction(x -> [-softplus(x[1]), softplus(x[2])])
)

rng_nn        = MersenneTwister(666)
nn_ps_, nn_st = Lux.setup(rng_nn, nn)
nn_ps         = ComponentArray(nn_ps_)

function lotka_volterra_ude!(du, u, p, t, p_true)
    x, y = u
    interaction, _ = nn([x, y], p, nn_st)
    du[1] =  p_true[1] * x + interaction[1]
    du[2] = -p_true[4] * y + interaction[2]
end

# Closure with the known parameters
nn_dynamics!(du, u, p, t) = lotka_volterra_ude!(du, u, p, t, p_true)

# ---------------------------------------------------------------------------
# predict / loss
# ---------------------------------------------------------------------------

prob_ude = ODEProblem(nn_dynamics!, u0, tspan, nn_ps)

function predict(ps)
    _prob = remake(prob_ude, u0=u0, tspan=tspan, p=ps)
    Array(solve(_prob, Vern7(), saveat=t_obs,
                abstol=1e-6, reltol=1e-6,
                sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))))
end

loss_history = Float64[]

function loss(ps)
    ŷ = predict(ps)
    size(ŷ, 2) != length(t_obs) && return Inf
    return mean((ŷ .- target) .^ 2)
end

# ---------------------------------------------------------------------------
# Training: Phase 1 — Adam, Phase 2 — BFGS
# ---------------------------------------------------------------------------

# Tunable training hyperparameters
adam_iters = 1500
bfgs_iters = 1000

# Optimization function compatible with Optimization.jl
optf = OptimizationFunction((ps, _) -> loss(ps), AutoZygote())

# Phase 1: Adam
println("Fase 1: Adam ($adam_iters iters)")
optprob1 = OptimizationProblem(optf, nn_ps)
res1 = Optimization.solve(optprob1, OptimizationOptimisers.Adam(),
                           callback = (ps, l) -> begin
                               push!(loss_history, l)
                               length(loss_history) % 50 == 0 &&
                                   println("  [Adam] iter $(length(loss_history))  —  Loss: $(round(l, digits=6))")
                               false
                           end,
                           maxiters = adam_iters)
println("  Loss final Adam: $(round(res1.objective, digits=6))\n")

# Phase 2: BFGS
println("Fase 2: BFGS ($bfgs_iters iters)\n")
optprob2 = OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(optprob2,
                           OptimizationOptimJL.BFGS(linesearch=BackTracking()),
                           callback = (ps, l) -> begin
                               push!(loss_history, l)
                               n = length(loss_history) - adam_iters
                               n % 10 == 0 &&
                                   println("  [BFGS] iter $n  —  Loss: $(round(l, digits=6))")
                               false
                           end,
                           maxiters = bfgs_iters)
println("  Loss final BFGS: $(round(res2.objective, digits=6))\n")

nn_ps = res2.u
println("Entrenamiento finalizado. Loss final: $(round(loss_history[end], digits=6))")

# ---------------------------------------------------------------------------
# Evaluate the trained UDE
# ---------------------------------------------------------------------------

prob_trained = ODEProblem(nn_dynamics!, u0, tspan, nn_ps)
sol_trained  = solve(prob_trained, Vern7(), saveat=0.1, abstol=1e-6, reltol=1e-6)
sol_gt_fine  = solve(prob_true, Tsit5(), saveat=0.1)

# ---------------------------------------------------------------------------
# Plot 1: trajectories — true vs UDE
# ---------------------------------------------------------------------------

p1 = plot(sol_gt_fine.t, sol_gt_fine[1, :], label="Presas (verdad)", linewidth=2,
          color=:blue, xlabel="Tiempo", ylabel="Población",
          title="UDE inverso: ajuste de trayectorias", legend=false)
plot!(p1, sol_gt_fine.t, sol_gt_fine[2, :], label="Depredadores (verdad)",
      linewidth=2, color=:red)
scatter!(p1, collect(t_obs), target[1, :], label="Presas (obs.)", markersize=3,
         color=:blue, alpha=0.5)
scatter!(p1, collect(t_obs), target[2, :], label="Depredadores (obs.)", markersize=3,
         color=:red, alpha=0.5)
plot!(p1, sol_trained.t, sol_trained[1, :], label="Presas (UDE entrenada)",
      linewidth=2, color=:dodgerblue, linestyle=:dash)
plot!(p1, sol_trained.t, sol_trained[2, :], label="Depredadores (UDE entrenada)",
      linewidth=2, color=:orangered, linestyle=:dash)

# Plot 2: training loss curve
p2 = plot(1:length(loss_history), loss_history, label="Loss (MSE)",
          xlabel="Época", ylabel="MSE", title="Curva de entrenamiento",
          linewidth=2, color=:purple, yscale=:log10, legend=false)

# Side legend panel (invisible axes, labels only)
p_legend = plot([], [], label="Presas (verdad)", color=:blue, linewidth=2,
                framestyle=:none, title="")
plot!(p_legend, [], [], label="Depredadores (verdad)", color=:red, linewidth=2)
scatter!(p_legend, [], [], label="Presas (obs.)", color=:blue, alpha=0.5, markersize=3)
scatter!(p_legend, [], [], label="Depredadores (obs.)", color=:red, alpha=0.5, markersize=3)
plot!(p_legend, [], [], label="Presas (UDE)", color=:dodgerblue, linewidth=2, linestyle=:dash)
plot!(p_legend, [], [], label="Depredadores (UDE)", color=:orangered, linewidth=2, linestyle=:dash)
plot!(p_legend, legend=:left, legendfontsize=8)

fig = plot(p1, p2, p_legend, layout=@layout([a b c{0.2w}]),
           size=(1200, 420), bottom_margin=10Plots.mm)
savefig(fig, "lotka_volterra_inverse_ude.png")
println("Figura guardada en lotka_volterra_inverse_ude.png")

# ---------------------------------------------------------------------------
# Plot 3: recovered interaction terms vs true interaction terms
#
# True:  f₁(x,y) = -β x y  and  f₂(x,y) = δ x y
# These are the terms the neural network should have learned.
# ---------------------------------------------------------------------------

# Sample (x, y) pairs along the true trajectory
X_traj = sol_gt_fine[1, :]
Y_traj = sol_gt_fine[2, :]

nn_out = reduce(hcat,
    [nn([x, y], nn_ps, nn_st)[1] for (x, y) in zip(X_traj, Y_traj)])

true_f1 = @. -p_true[2] * X_traj * Y_traj   # -β x y
true_f2 = @.  p_true[3] * X_traj * Y_traj   #  δ x y

p3 = plot(sol_gt_fine.t, true_f1, label="f₁ verdadera (−βxy)", linewidth=2, color=:blue,
          xlabel="Tiempo", ylabel="Término de interacción",
          title="Términos de interacción aprendidos")
plot!(p3, sol_gt_fine.t, true_f2, label="f₂ verdadera (δxy)", linewidth=2, color=:red)
plot!(p3, sol_gt_fine.t, nn_out[1, :], label="f₁ (red neuronal)",
      linewidth=2, color=:dodgerblue, linestyle=:dash)
plot!(p3, sol_gt_fine.t, nn_out[2, :], label="f₂ (red neuronal)",
      linewidth=2, color=:orangered, linestyle=:dash)

savefig(p3, "lotka_volterra_inverse_ude_interactions.png")
println("Figura de interacciones guardada en lotka_volterra_inverse_ude_interactions.png")

# ---------------------------------------------------------------------------
# Plot 4: NN outputs as 2D heatmaps over the (x, y) state space
#
# f₁(x,y) = nn(x,y)[1]  — should approximate -β x y
# f₂(x,y) = nn(x,y)[2]  — should approximate  δ x y
#
# Seismic colormap: blue = negative, white = zero, red = positive (symmetric).
# ---------------------------------------------------------------------------

# Grid tightly around the trajectory with a small margin
margin = 0.2
x_range = range(minimum(sol_gt_fine[1, :]) - margin, maximum(sol_gt_fine[1, :]) + margin, length=80)
y_range = range(minimum(sol_gt_fine[2, :]) - margin, maximum(sol_gt_fine[2, :]) + margin, length=80)

# Mask: NaN for grid points whose nearest trajectory point is beyond threshold
traj_x = sol_gt_fine[1, :]
traj_y = sol_gt_fine[2, :]
mask_threshold = 0.3   # distance in state-space units

function masked_nn_output(x_range, y_range, output_idx)
    [begin
        d = minimum((x - tx)^2 + (y - ty)^2 for (tx, ty) in zip(traj_x, traj_y))
        sqrt(d) < mask_threshold ? nn([x, y], nn_ps, nn_st)[1][output_idx] : NaN
    end
    for y in y_range, x in x_range]
end

nn_f1 = masked_nn_output(x_range, y_range, 1)
nn_f2 = masked_nn_output(x_range, y_range, 2)

cmap1 = cgrad([:blue, :white])   # f₁: negative → zero
cmap2 = cgrad([:white, :red])    # f₂: zero → positive

clim1 = (-6.0, 0.0)
clim2 = (0.0, 6.0)

p4 = heatmap(x_range, y_range, nn_f1,
             xlabel="x (presas)", ylabel="y (depredadores)",
             title="NN output 1  [≈ −βxy]",
             color=cmap1, clims=clim1)
plot!(p4, sol_gt_fine[1, :], sol_gt_fine[2, :],
      color=:black, linewidth=1.5, label="", alpha=0.7)
scatter!(p4, [u0[1]], [u0[2]], color=:black, markersize=5, label="")

p5 = heatmap(x_range, y_range, nn_f2,
             xlabel="x (presas)", ylabel="y (depredadores)",
             title="NN output 2  [≈ δxy]",
             color=cmap2, clims=clim2)
plot!(p5, sol_gt_fine[1, :], sol_gt_fine[2, :],
      color=:black, linewidth=1.5, label="", alpha=0.7)
scatter!(p5, [u0[1]], [u0[2]], color=:black, markersize=5, label="")

fig_nn = plot(p4, p5, layout=(1, 2), size=(1000, 420), bottom_margin=8Plots.mm)
savefig(fig_nn, "lotka_volterra_inverse_ude_nn_heatmap.png")
println("Figura de heatmaps guardada en lotka_volterra_inverse_ude_nn_heatmap.png")
