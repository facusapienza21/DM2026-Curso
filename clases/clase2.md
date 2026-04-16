---
title: Clase No2 - ODEs y NODEs
---

Algunos comentarios...

# ODEs y NODEs

**Fecha:** 13/04/2026

:::{iframe} https://www.youtube.com/embed/kgRSMKC8Rrg
:width: 100%
:::

% Incluir espacio

## Ecuaciones diferenciales ordinarias (ODEs)

Un ejemplo de como escribir una {term}`ODE`:

$$
\frac{du}{dt} = f(u, t, \theta), \quad u(t_0) = u_0
$$

donde $u(t) \in \mathbb{R}^n$ es el estado del sistema, $\theta$ son los {term}`parámetros <Parámetro>` del modelo, y $u_0$ es la {term}`condición inicial <Condición inicial>`. La función $f$ es el {term}`campo vectorial <Campo vectorial>`.

:::{note} Ecuaciones diferenciales con derivadas de mayor orden
...
:::

:::{note} Ecuaciones diferenciales en derivadas parciales
...
:::

### Ejemplos

% Lotka-Volterra

## Inferencia estadística

% Modelo para el estado, modelo para el ruido
% Errores gaussianos
% Ajuste de trayectorias (trajectory matching)

### Ejemplo (continuado)

% Como se ve esto en el caso del sistema Lotka-Volterra?


## Ecuaciones diferenciales ordinarias neuronales (NODEs)

Las {term}`NODE`s fueron introducidas por {cite}`chen2018neural` ...

### Ejemplo (implementación computacional)

En Myst, podemos incluir codigo!

```julia
using DifferentialEquations
using Plots
```
