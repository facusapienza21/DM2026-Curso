---
title: Clase No4 - NODEs y Solvers Numericos para ODEs
---

# NODEs y Solvers Numericos para ODEs

**Fecha:** 20/04/2026

---

## Video de la clase

:::{iframe} https://www.youtube.com/embed/VIDEO_ID
:width: 100%
:::

---

## Resolucion numerica de ODEs

Queremos resolver numericamente, para un $\theta$ fijo:

$$
\dfrac{du}{dt} = f(u,t)
$$

Hay dos metodos generales que engloban al resto:

- [Multi-step](#metodos-multi-step)
- [Runge Kutta](#metodos-runge-kutta)

:::{note}
Se debe discretizar el eje temporal, $\Delta t$ se puede fijar o determinar dinamicamente por un algoritmo, pero ninguno de los metodos siguientes definen formas de discretizacion temporal.
:::

:::{note}
Se nota $u^m = u(t^m)$
:::

### Metodos Multi-step

$$
\sum_{i=1}^{d_1} \alpha_{i} u^{m+1} = \sum_{j=1}^{d_2} \beta_{i} f(u^{m+j}, t^{m+j})
$$

con $d_1, d_2 \in \mathbb{N}$ hiperparametros del solver.

#### Ejemplo: Adams-Bashforth

$$
u^{m+1} = u^m + \dfrac{\Delta t}{2} \left(3 f(u^m, t^m) - f(u^{m-1}, t^{m-1}) \right)
$$

### Metodos Runge-Kutta

$$
u^{m+1} &= u^m + \sum_{i=1}^{s} b_i k_i 
$$
$$
k_i &= f\left(u^m + \sum_{j=1}^{s} a_{ij} k_j, t^i + c_i \Delta t \right)
$$

Con $s \in \mathbb{N}$ hiperparametro del solver.

#### Ejemplo: Metodo del punto medio

Para $s=2$:

$$
k_1 &= f(u^m, t^m) \\
k_2 &= f\left(u^m + \dfrac{\Delta t}{2} k_1, t^m + \dfrac{\Delta t}{2} \right) \\
u^{m+1} &= u^m + \Delta t k_2
$$

### Euler explicito

Euler explicito es un metodo de [Runge-Kutta](#metodos-runge-kutta) con $s=1$, o un metodo [multi-step](#metodos-multi-step) con $d_1=1, d_2=0$:
$$
u^{m+1} = u^m + \Delta t f(u^m, t^m)
$$