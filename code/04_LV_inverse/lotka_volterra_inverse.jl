using DifferentialEquations
using Optim
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
    x, y = u
    α, β, δ, γ = p
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
end

# -----------------------------------------------------------------------
# Ground truth: solve with known parameters and add Gaussian noise
# -----------------------------------------------------------------------

p_true = [1.0, 0.1, 0.075, 1.5]   # [α, β, δ, γ]
u0     = [10.0, 5.0]               # initial conditions
tspan  = (0.0, 30.0)
t_obs  = range(0.0, 30.0, step=0.5) # observation times

prob_true = ODEProblem(lotka_volterra!, u0, tspan, p_true)
sol_true  = solve(prob_true, Tsit5(), saveat=t_obs)

σ_noise = 1.0
rng     = MersenneTwister(42)
observations = Array(sol_true) .+ σ_noise .* randn(rng, size(Array(sol_true)))

# -----------------------------------------------------------------------
# Loss function: MSE between model prediction and noisy observations
# -----------------------------------------------------------------------

# History recorded on every loss evaluation
param_history = Vector{Vector{Float64}}()
loss_history  = Vector{Float64}()

function loss(p)
    any(p .<= 0) && return Inf

    prob = ODEProblem(lotka_volterra!, u0, tspan, p)
    sol  = solve(prob, Tsit5(), saveat=t_obs)

    length(sol.t) != length(t_obs) && return Inf

    val = mean((Array(sol) .- observations) .^ 2)
    push!(param_history, copy(p))
    push!(loss_history, val)
    return val
end

# -----------------------------------------------------------------------
# Optimization with Nelder-Mead
# -----------------------------------------------------------------------

# Initial guess: perturbed from the true values
p0 = p_true .* [1.3, 0.7, 1.5, 0.8]
println("Parámetros verdaderos:  α=$(p_true[1])  β=$(p_true[2])  δ=$(p_true[3])  γ=$(p_true[4])")
println("Estimación inicial:     α=$(p0[1])  β=$(p0[2])  δ=$(p0[3])  γ=$(p0[4])")
println("Loss inicial: $(loss(p0))\n")

result = optimize(loss, p0, NelderMead(),
                  Optim.Options(iterations=200, show_trace=false))

p_opt = Optim.minimizer(result)
println("Optimización finalizada: $(Optim.converged(result) ? "convergió" : "no convergió")")
println("Loss final: $(Optim.minimum(result))\n")
println("Parámetros recuperados:")
println("  α: $(p_opt[1])  (verdadero: $(p_true[1]))")
println("  β: $(p_opt[2])  (verdadero: $(p_true[2]))")
println("  δ: $(p_opt[3])  (verdadero: $(p_true[3]))")
println("  γ: $(p_opt[4])  (verdadero: $(p_true[4]))")

# -----------------------------------------------------------------------
# Plot: ground truth, noisy observations, and recovered solution
# -----------------------------------------------------------------------

sol_opt = solve(ODEProblem(lotka_volterra!, u0, tspan, p_opt), Tsit5(), saveat=0.1)
sol_gt  = solve(prob_true, Tsit5(), saveat=0.1)

p1 = plot(sol_gt.t, sol_gt[1, :], label="Presas (verdad)", linewidth=2, color=:blue,
          xlabel="Tiempo", ylabel="Población", title="Ajuste de trayectorias: Lotka-Volterra")
plot!(p1, sol_gt.t, sol_gt[2, :], label="Depredadores (verdad)", linewidth=2, color=:red)
scatter!(p1, collect(t_obs), observations[1, :], label="Presas (obs.)", markersize=3,
         color=:blue, alpha=0.5)
scatter!(p1, collect(t_obs), observations[2, :], label="Depredadores (obs.)", markersize=3,
         color=:red, alpha=0.5)
plot!(p1, sol_opt.t, sol_opt[1, :], label="Presas (ajustado)", linewidth=2,
      color=:dodgerblue, linestyle=:dash)
plot!(p1, sol_opt.t, sol_opt[2, :], label="Depredadores (ajustado)", linewidth=2,
      color=:orangered, linestyle=:dash)

p2 = bar(["α", "β", "δ", "γ"],
         [p_opt[i] / p_true[i] for i in 1:4],
         xlabel="Parámetro", ylabel="Estimado / Verdadero",
         title="Recuperación de parámetros", legend=false,
         color=[:blue, :red, :green, :orange], ylims=(0, 2))
hline!(p2, [1.0], color=:black, linewidth=2, linestyle=:dash, label="Valor exacto")

fig = plot(p1, p2, layout=(1, 2), size=(1000, 420), bottom_margin=10Plots.mm)
savefig(fig, "lotka_volterra_inverse.png")
println("\nFigura guardada en lotka_volterra_inverse.png")

# -----------------------------------------------------------------------
# Animation: trajectory approximation across iterations
# -----------------------------------------------------------------------

# Build best-so-far history: at each evaluation take the params with lowest loss seen
best_loss_so_far, best_params = let
    best_loss = Inf
    best_p    = param_history[1]
    bls = Vector{Float64}(undef, length(loss_history))
    bps = Vector{Vector{Float64}}(undef, length(param_history))
    for i in eachindex(param_history)
        if loss_history[i] < best_loss
            best_loss = loss_history[i]
            best_p    = param_history[i]
        end
        bls[i] = best_loss
        bps[i] = best_p
    end
    bls, bps
end

# Subsample to keep the animation manageable
n_frames  = min(80, length(best_params))
frame_idx = round.(Int, range(1, length(best_params), length=n_frames))

sol_gt_fine = solve(prob_true, Tsit5(), saveat=0.1)
t_fine      = sol_gt_fine.t

anim = @animate for i in frame_idx
    p_i    = best_params[i]
    loss_i = best_loss_so_far[i]

    sol_i = solve(ODEProblem(lotka_volterra!, u0, tspan, p_i), Tsit5(), saveat=0.1)

    frame = plot(t_fine, sol_gt_fine[1, :], label="Presas (verdad)", linewidth=2,
                 color=:blue, xlabel="Tiempo", ylabel="Población",
                 title="Iter $(frame_idx[findfirst(==(i), frame_idx)])  |  Loss = $(round(loss_i, digits=3))",
                 ylims=(0, 40), legend=:topright)
    plot!(frame, t_fine, sol_gt_fine[2, :], label="Depredadores (verdad)",
          linewidth=2, color=:red)
    scatter!(frame, collect(t_obs), observations[1, :], label="", markersize=2,
             color=:blue, alpha=0.4)
    scatter!(frame, collect(t_obs), observations[2, :], label="", markersize=2,
             color=:red, alpha=0.4)

    if length(sol_i.t) == length(t_fine)
        plot!(frame, sol_i.t, sol_i[1, :], label="Presas (modelo)", linewidth=2,
              color=:dodgerblue, linestyle=:dash)
        plot!(frame, sol_i.t, sol_i[2, :], label="Depredadores (modelo)", linewidth=2,
              color=:orangered, linestyle=:dash)
    end

    frame
end

gif(anim, "lotka_volterra_training.gif", fps=5)
println("Animación guardada en lotka_volterra_training.gif")
