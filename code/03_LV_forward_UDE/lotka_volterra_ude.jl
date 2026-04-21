using DifferentialEquations
using Lux
using Optimisers
using Plots
using Statistics
using Zygote
using Random

ENV["GKSwstype"] = "100"
Plots.default(show=false)

# ---------------------------------------------------------------------------
# Flag: set to true to pretrain the neural network before the UDE forward pass
# ---------------------------------------------------------------------------
PRETRAIN = false

# ---------------------------------------------------------------------------
# True Lotka-Volterra system (used to generate the reference solution)
# ---------------------------------------------------------------------------

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
end

α = 1.0;   β = 0.1;   δ = 0.075;   γ = 1.5
p_true = [α, β, δ, γ]
u0     = [10.0, 5.0]
tspan  = (0.0, 30.0)

prob_true = ODEProblem(lotka_volterra!, u0, tspan, p_true)
sol       = solve(prob_true, Tsit5(), saveat=0.1)

# ---------------------------------------------------------------------------
# Universal Differential Equation (UDE)
#
# We replace the interaction terms with a small neural network nn(u):
#
#   dx/dt = α x + nn(u)[1]
#   dy/dt = -γ y + nn(u)[2]
#
# In the true system, nn(u) = [-β x y,  δ x y], so the network should
# learn the nonlinear predator-prey coupling from data.
# ---------------------------------------------------------------------------

# Small neural network: R² → R²
# Input: [x, y]  |  Output: interaction terms for each equation
nn = Chain(
    Dense(2 => 5, tanh),
    Dense(5 => 5, tanh),
    Dense(5 => 2)
)

rng_nn       = MersenneTwister(0)
nn_ps, nn_st = Lux.setup(rng_nn, nn)

# ---------------------------------------------------------------------------
# Pretraining: supervise the network to approximate f(x,y) = (-xy, xy)
#
# We sample (x, y) pairs from the range seen in the true solution and
# train the network to match the target interaction terms before embedding
# it in the ODE.
# ---------------------------------------------------------------------------

if PRETRAIN
    nn_ps = let nn_ps = nn_ps
        X_mat = reduce(hcat, [[x, y] for x in range(0.0, 60.0, length=30)
                                      for y in range(0.0, 60.0, length=30)])
        Y_mat = reduce(hcat, [[-0.1 * x * y, 0.075 * x * y]
                               for x in range(0.0, 60.0, length=30)
                               for y in range(0.0, 60.0, length=30)])

        pretrain_loss(ps) = mean((nn(X_mat, ps, nn_st)[1] .- Y_mat) .^ 2)

        opt_state = Optimisers.setup(Adam(1e-3), nn_ps)
        println("Pre-entrenando la red neuronal...")
        for epoch in 1:10000
            grads = Zygote.gradient(pretrain_loss, nn_ps)[1]
            opt_state, nn_ps = Optimisers.update(opt_state, nn_ps, grads)
            if epoch % 100 == 0
                println("  Época $epoch — Loss: $(round(pretrain_loss(nn_ps), digits=6))")
            end
        end
        println("Preentrenamiento finalizado. Loss final: $(round(pretrain_loss(nn_ps), digits=6))\n")
        nn_ps
    end
else
    println("Preentrenamiento omitido — usando pesos aleatorios.\n")
end

# ---------------------------------------------------------------------------
# UDE forward pass
# ---------------------------------------------------------------------------

function lotka_volterra_ude!(du, u, p, _)
    α, γ, nn_ps, nn_st = p
    x, y = u
    interaction, _ = nn([x, y], nn_ps, nn_st)
    du[1] = α * x + interaction[1]
    du[2] = -γ * y + interaction[2]
end

p_ude    = (α, γ, nn_ps, nn_st)
prob_ude = ODEProblem(lotka_volterra_ude!, u0, tspan, p_ude)
sol_ude  = solve(prob_ude, Tsit5(), saveat=0.1)

fig_ude = plot(sol.t, sol[1, :], label="Presas (verdad)", linewidth=2, color=:blue,
               xlabel="Tiempo", ylabel="Población",
               title="UDE — " * (PRETRAIN ? "pesos preentrenados" : "pesos aleatorios"))
plot!(fig_ude, sol.t, sol[2, :], label="Depredadores (verdad)", linewidth=2, color=:red)
plot!(fig_ude, sol_ude.t, sol_ude[1, :], label="Presas (UDE)", linewidth=2,
      color=:dodgerblue, linestyle=:dash)
plot!(fig_ude, sol_ude.t, sol_ude[2, :], label="Depredadores (UDE)", linewidth=2,
      color=:orangered, linestyle=:dash)

savefig(fig_ude, "lotka_volterra_ude.png")
println("Figura UDE guardada en lotka_volterra_ude.png")
