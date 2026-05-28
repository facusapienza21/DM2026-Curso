---
title: No9 - PINNs Pt2
---

# Physics-Informed Neural Networks (PINNs)

**Fecha:** 11/05/2026

:::{iframe} https://www.youtube.com/embed/fr6iQbl_XfA
:width: 100%
:::

:::{seealso}
Esta clase continúa directamente desde la [Clase 8](clase8.md), donde se introdujeron las PINNs.
:::

## Tres tipos de algoritmos para incorporar ecuaciones diferenciales

Cuando queremos estimar o aproximar una ecuación diferencial con modelos de aprendizaje automático, hay tres enfoques principales:

1. **Emuladores:** se estima la ecuación diferencial con otro modelo de ML.
Por ejemplo, cuando hay ecuaciones costosas (Navier-Stokes, etc.)

2. **Restricciones suaves** (vía Lagrangiano): el modelo se entrena minimizando una función de costo empírica $\mathcal{L}_{\text{emp}}$ más un término de penalización $D[x(\theta)]$ que _incentiva_ satisfacer la ecuación diferencial.
Por ejemplo, las {term}`PINN`s:

$$
\min_{\theta,\, x} \mathcal{L}_{\text{emp}}(y, x, \theta) + \lambda \| D[x(\theta)] \|.
$$

3. **Restricciones fuertes** (se satisfacen _estrictamente_): el modelo debe satisfacer la ecuación diferencial $D[x(\theta)] = 0$ en todo momento.
Son las **NODE** y **UDE**:

$$
\min_{\theta,\, x} \mathcal{L}_{\text{emp}}(y, x(\theta), \theta) \quad \text{s.a. } D[x(\theta)] = 0.
$$

## Dualidad Lagrangiana

Un resultado fundamental de la optimización con restricciones establece que, si relajamos la restricción fuerte $D[x(\theta)] = 0$ por una restricción aproximada $\| D[x(\theta)] \| \leq \varepsilon$, existe un $\lambda(\varepsilon)$ tal que la solución del problema con restricción suave (penalizado con $\lambda(\varepsilon)$) coincide con la solución del problema con restricción fuerte relajada {cite}`boyd2004convex`:

$$
\lambda(\varepsilon) = \lambda \quad \Longleftrightarrow \quad \text{restricción suave} = \text{restricción fuerte relajada}
$$

:::{note} Hard PINN
Elegir $\lambda \to \infty$ en la formulación suave es equivalente a hacer $\varepsilon \to 0$, es decir, a forzar la restricción estrictamente.
Las Hard PINNs imponen la ecuación diferencial de forma exacta reparametrizando la solución directamente, en lugar de penalizarla {cite}`Lu_Pestourie_Yao_Wang_Verdugo_Johnson_2021`.
:::

## Dificultades al entrenar PINNs

En la práctica, entrenar PINNs tiene problemas importantes de optimización.
La función de pérdida puede tener múltiples términos: $\mathcal{L}_{\text{emp}}$ (ajuste a los datos) y uno o mas términos de penalización $\lambda_i \|D_i[x(\theta)]\|$, uno por cada restricción (condición inicial, condición de borde, residuo de la ecuación diferencial, etc.).

$$
\mathcal{L}_{\text{PINN}} = \mathcal{L}_{\text{emp}} + \sum_{i=1}^{n} \lambda_i \|D_i[x(\theta)]\|
$$

Cuando $\lambda$ es muy grande, el problema se vuelve **mal condicionado**: los gradientes del término de penalización dominan, entonces el optimizador tiene problemas para elegir la dirección óptima.

```{figure} ./figures/no9_pinns_mal_condicionadas.png
:width: 75%
:align: center

Curvas de nivel de $\mathcal{L}_{\text{PINN}}$ (rosa) alrededor de la variedad donde $D[x(\cdot,\theta)]=0$ (negro). En un entorno de $D[x(\cdot,\theta)]=0$ es probable que la PINN esté mal condicionada. 
```

:::{note} Condición de descenso de gradiente
Una forma de cuantificar cuán mal condicionado está el problema de optimización es analizar el {term}`número de condición <Número de condición>` de la matriz Hessiana $H = \nabla^2 \mathcal{L}$:

$$
\kappa(H) = \frac{\lambda_{\max}(H)}{\lambda_{\min}(H)}
$$

Un $\kappa(H)$ grande indica que la función de pérdida tiene direcciones con curvatura muy alta y otras con curvatura muy baja, lo que se traduce en curvas de nivel elongadas (celeste en el diagrama).
:::

### Selección de $\lambda$

En general, el vector de hiperparámetros $\lambda = (\lambda_1, \ldots, \lambda_n)$ controla cuánto pesa cada restricción. Hay tres estrategias para elegirlos:

* **Adimensionalización:** normalizar todos los términos para que operen en la misma escala:

$$
\mathcal{L}_{\text{emp}}^k \approx \lambda_1^k \mathcal{L}_{CI} \approx \lambda_2^k \mathcal{L}_{\partial} \approx \cdots
$$

* **Elección democrática — valor absoluto:** elegir los $\lambda^k$ tal que las magnitudes de todos los términos sean iguales en cada iteración:
$$
\left| \mathcal{L}_{\text{emp}}^k \right| \approx \left| \lambda_1^k \mathcal{L}_{CI} \right| \approx \left| \lambda_2^k \mathcal{L}_{\partial} \right| \approx \cdots
$$

* **Elección democrática — norma de gradientes:** igualar las normas de los gradientes de cada término, de forma que ninguno domine la actualización:

$$
\| \nabla_\theta \mathcal{L}_{\text{emp}}^k \| \approx \| \lambda_1^k \nabla_\theta \mathcal{L}_{CI} \| \approx \| \lambda_2^k \nabla_\theta \mathcal{L}_{\partial} \| \approx \cdots
$$


## Implementación

### Problema directo

Vemos una implementación de PINN directo.

```julia
t0, T = tspan
scale_layer_inp = WrappedFunction(t -> (t .- t0) ./ (T - t0))
scale_layer_out = WrappedFunction(x -> 10.0 .* x)
```

:::{note}
Se aplica escalado a la red porque las redes con bias tienden a ajustar mejor funciones de alta frecuencia que de baja frecuencia.
El escalado busca corregir este {term}`sesgo espectral <Sesgo espectral>`.
:::

Definimos la red neuronal:
```julia
nn = Chain(
    scale_layer_inp,
    Dense(1 => 5, tanh),
    Dense(5 => 10, tanh),
    Dense(10 => 10, tanh),
    Dense(10 => 2),
    scale_layer_out
)
```

**Puntos de colocación:**
En una PINN hay que elegir puntos en el dominio temporal donde evaluar el residuo de la ecuación diferencial $D[x(\theta)]$.
```julia
N_col = 100
# t_col = collect(range(tspan[1], tspan[2], length=N_col))
# t_col = 10.0.^collect(range(log10(tspan[1] + 1e-10), log10(tspan[2]), length = N_col))
t_col = [
    10.0.^collect(range(log10(tspan[1] + 1e-10), 
	log10(tspan[2]), length = div(N_col, 2)));
    collect(range(tspan[1], tspan[2], length=div(N_col, 2)))
]
```
```{figure} ../code/06_LV_forward_PINN/lv_pinn_solution.png
:width: 75%
:align: center
:figclass: text-center

Las cruces son los puntos de colocación: donde la ecuación diferencial se está evaluando.
```
Este es un ejemplo que claramente no convergió bien, probablemente el resultado mejoraría ejecutando el código por mas tiempo, pero.. ¿Qué mas podemos hacer?


#### Técnicas para mejorar la convergencia

**Puntos de colocación:**
La ubicación de los puntos de colocación no es arbitraria. En este ejemplo, exploramos solo dos opciones:

- **Uniforme:** puntos distribuidos uniformemente en el dominio.
- **Escala logarítmica desde $t_0$:** mayor densidad de puntos cerca de la condición inicial. Una solución espuria que la red puede encontrar es arrancar en $u_0$ e irse inmediatamente a la trayectoria nula, lo que satisface trivialmente las ecuaciones de Lotka-Volterra. Concentrar puntos al principio fuerza al modelo a seguir la trayectoria correcta desde $t_0$.

**Forzar la condición inicial exactamente:**
Vemos cómo implementamos la función de costo:
```julia
function pinn_loss_components(ps)
    u_ic    = nn([tspan[1]], ps, st)[1]
    ic_loss = sum((u_ic .- u0) .^ 2)

    phys_residuals = map(t_col) do t_i
        u    = nn([t_i], ps, st)[1]
        dudt = nn_dudt(t_i, ps, st)
        rhs  = lv_rhs(u)
        sum((dudt .- rhs) .^ 2)
    end
    phys_loss_weighted = λ_physics * mean(phys_residuals)

    return ic_loss, phys_loss_weighted
end

pinn_loss(ps) = sum(pinn_loss_components(ps))
```

Estamos definiendo a la condición inicial como una restricción suave. No es ideal. ¿Cómo podemos imponerla como restricción fuerte? Una reparametrización directa es:

$$
u(t) = \text{NN}_\theta(t)(t - t_0) + u_0
$$

En $t = t_0$ se tiene $u(t_0) = \text{NN}_\theta(t_0) \cdot 0 + u_0 = u_0$ independientemente de los parámetros. Sin embargo, el factor $(t - t_0)$ crece linealmente, obligando a la red a desaprender un trend lineal.

Una mejor alternativa es usar una función $\phi$ acotada:

$$
u_\theta(t) = u_0 + \phi(t)\, \text{NN}_\theta(t)
$$

donde $\phi(t_0) = 0$. Por ejemplo, $\phi(t) = \frac{t - t_0}{t - T_0}$ con $T_0$ cualquier número finito satisface $\phi(t_0) = 0$ y $\phi(t) \to 1$ cuando $t \to \infty$.


Podemos aplicar esto al ejemplo:
```julia
function u_ansatz(t_val, ps, st)
    raw = nn([t_val], ps, st)[1]   # NN(t) ∈ ℝ²
    return u0 .+ t_val .* raw      # u₀ + t·NN(t)
end

# Derivada temporal del ansatz por diferencias finitas centrales
function u_ansatz_dudt(t_val, ps, st; h=1e-4)
    u_fwd = u_ansatz(t_val + h, ps, st)
    u_bwd = u_ansatz(t_val - h, ps, st)
    return (u_fwd .- u_bwd) ./ (2h)
end
function pinn_loss(ps)
    phys_residuals = map(t_col) do t_i
        u    = u_ansatz(t_i, ps, st)
        dudt = u_ansatz_dudt(t_i, ps, st)
        rhs  = lv_rhs(u)
        sum((dudt .- rhs) .^ 2)
    end
    return mean(phys_residuals)
end
```

```{figure} ../code/07_LV_forward_PINN_hard_IC/lv_pinn_hardIC_solution.png
:width: 75%
:align: center
:figclass: text-center

Mejora... pero no es ideal.
```

Aún podría ser mejorable este resultado. Pero no es lo que nos interesa, no buscamos que una PINN sea excelente resolviendo una ecuación diferencial. Donde verdaderamente brillan las PINN es en el problema inverso.


### Problema inverso

En el problema inverso no solo se conoce la condición inicial sino también observaciones de la trayectoria. El objetivo deja de ser únicamente ajustar $u(t_0) = u_0$ y pasa a ser ajustar la trayectoria completa. En este caso se aplica a Lotka-Volterra sin conocer los parámetros del sistema.


#### Generación de datos observados

Se resuelve la ODE con los parámetros verdaderos (desconocidos para la PINN) y se submuestran $M = 40$ puntos uniformes:

```julia
const α_true = 1.0;   const β_true = 0.1
const δ_true = 0.075; const γ_true = 1.5
const u0    = [10.0, 5.0]
const tspan = (0.0, 15.0)

prob_true = ODEProblem(lotka_volterra!, u0, tspan, [α_true, β_true, δ_true, γ_true])
sol_true  = solve(prob_true, Tsit5(), saveat=0.01)

M_obs = 40
t_obs = collect(range(tspan[1], tspan[2], length=M_obs))
u_obs = hcat([sol_true(t) for t in t_obs]...)   # (2, M_obs)
```

#### Parámetros entrenables

Los parámetros de la ODE se inicializan lejos de los valores verdaderos y se agrupan junto con los pesos de la red en un único `NamedTuple` que `Optimisers.jl` puede diferenciar recursivamente:

```julia
p_ode_init = [0.5, 0.05, 0.04, 1.0]   # [α̂, β̂, δ̂, γ̂]

theta = (nn = nn_ps, ode = copy(p_ode_init))
```

#### Función de pérdida

```julia
const λ_phys = 1.0

function inverse_loss_components(theta)
    nn_ps = theta.nn
    α̂, β̂, δ̂, γ̂ = theta.ode

    data_residuals = map(1:M_obs) do j
        u_pred = nn([t_obs[j]], nn_ps, st)[1]
        sum((u_pred .- u_obs[:, j]) .^ 2)
    end
    data_loss = mean(data_residuals)

    phys_residuals = map(t_col) do t_i
        u    = nn([t_i], nn_ps, st)[1]
        dudt = nn_dudt(t_i, nn_ps, st)
        x, y = u[1], u[2]
        rhs  = [α̂*x - β̂*x*y,
                δ̂*x*y - γ̂*y]
        sum((dudt .- rhs) .^ 2)
    end
    phys_loss = λ_phys * mean(phys_residuals)

    return data_loss, phys_loss
end

inverse_loss(theta) = sum(inverse_loss_components(theta))
```


#### Resultados

```{figure} ../code/08_LV_inverse_PINN/lv_inverse_pinn_params.png
:width: 75%
:align: center

Convergencia de los parámetros estimados $\hat{p}$ hacia los valores verdaderos a lo largo del entrenamiento.
```

El error puntual muestra que la PINN ajusta exactamente en los puntos de colocación y en los puntos de observación, pero el error crece en las regiones intermedias.

```{figure} ../code/08_LV_inverse_PINN/lv_inverse_pinn_error.png
:width: 75%
:align: center
```

:::{important}
En este ejemplo se usó una grilla fija de puntos de colocación. En la práctica, conviene resamplear la grilla en cada iteración. Hay dos estrategias:

- **Resampleo uniforme:** nuevos puntos aleatorios en cada paso.
- **Important sampling:** concentrar puntos donde el residuo de la ecuación diferencial es mayor, es decir, donde la PINN está errando más.
:::
