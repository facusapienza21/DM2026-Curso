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

:::{note}
Se aplica escalado a la red porque las redes con bias tienden a ajustar mejor funciones de alta frecuencia que de baja frecuencia.
El escalado busca corregir este {term}`sesgo espectral <Sesgo espectral>`.
:::

#### Técnicas para mejorar la convergencia

**Puntos de colocación:**
En una PINN hay que elegir puntos en el dominio temporal donde evaluar el residuo de la ecuación diferencial $D[x(\theta)]$. En este ejemplo, exploramos solo dos opciones de hacer esto:

- **Uniforme:** puntos distribuidos uniformemente en el dominio.
- **Escala logarítmica desde $t_0$:** mayor densidad de puntos cerca de la condición inicial. Una solución espuria que la red puede encontrar es arrancar en $u_0$ e irse inmediatamente a la trayectoria nula, lo que satisface trivialmente las ecuaciones de Lotka-Volterra. Concentrar puntos al principio fuerza al modelo a seguir la trayectoria correcta desde $t_0$.

**Forzar la condición inicial exactamente:**
La condición inicial puede imponerse como restricción fuerte en lugar de suave. Una reparametrización directa es:

$$
u(t) = \text{NN}_\theta(t)(t - t_0) + u_0
$$

En $t = t_0$ se tiene $u(t_0) = \text{NN}_\theta(t_0) \cdot 0 + u_0 = u_0$ independientemente de los parámetros. Sin embargo, el factor $(t - t_0)$ crece linealmente, obligando a la red a desaprender un trend lineal.

Una mejor alternativa es usar una función $\phi$ acotada:

$$
u_\theta(t) = u_0 + \phi(t)\, \text{NN}_\theta(t)
$$

donde $\phi(t_0) = 0$. Por ejemplo, $\phi(t) = \frac{t - t_0}{t - T_0}$ con $T_0$ cualquier número finito satisface $\phi(t_0) = 0$ y $\phi(t) \to 1$ cuando $t \to \infty$.

### Problema inverso

En el problema inverso no solo se conoce la condición inicial sino también observaciones de la trayectoria. El objetivo deja de ser únicamente ajustar $u(t_0) = u_0$ y pasa a ser ajustar la trayectoria completa. En este caso se aplica a Lotka-Volterra sin conocer los parámetros del sistema.

```{figure} ../code/08_LV_inverse_PINN/lv_inverse_pinn_params.png
:width: 75%
:align: center
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
