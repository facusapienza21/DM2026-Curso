---
title: Clase No4 - Metodos Numericos
---

# Metodos Numericos

**Fecha:** 20/04/2026

---

## Video de la clase

:::{iframe} https://www.youtube.com/watch?v=Rdvz8KRA1JQ 
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
Se debe discretizar el eje temporal, y asi obtener $t^0, t^1, \ldots, t^M$, donde $\Delta t$ se puede fijar o determinar dinamicamente por un algoritmo, pero ninguno de los metodos siguientes definen formas de discretizacion temporal.
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

## Motivacion de las NODEs

En {cite}`chen2018neural`, se presenta la idea de usar una red neuronal para parametrizar la funcion $f$ de una EDO dada su similitud con el metodo de [Euler explicito](#euler-explicito).:

Si $g(h_i) = \sigma(W_i h_i + b)$ es una red neuronal con parametros $\theta = (W, b)$, entonces

$$
h_{i+1} = h_i + g(h_i) 	\rightsquigarrow  h_{i+1} + \Delta t \tilde{g}(h_i)
$$

Y aqui se observa la similitud con el metodo de [Euler explicito](#euler-explicito), lo que motiva la idea de "pasar al continuo" y definir una NODE como la solucion de la siguiente EDO:
$$
\dfrac{dh}{dt} = \tilde{g}(h(t)) \\
\dfrac{dh}{dt} = NN_\theta(h(t))
$$

## Pre-Entrenamiento

A veces resulta necesario pre-entrenar la red neuronal que define la funcion $f$ de la EDO, para "acercar" la solucion a un espacio donde la solucion tiene sentido (fisico).

Si tenemos una idea de la fisica subyacente dada por $\tilde{F}$, entonces podemos usar esta informacion para guiar el entrenamiento de la red neuronal optimizando primero la siguiente funcion de perdida:

$$
\mathcal{L}_{pre}(\theta) = \sum_{i=1}^{N} \|\tilde{F}(X) - NN_\theta(X)\|^2
$$

Esto es mucho menos costoso que optimizar la funcion de perdida que se obtiene al comparar la solucion de la EDO con los datos.
Observemos que no hay datos en esta etapa, solo la informacion fisica dada por $\tilde{F}$.

:::{bibliography}
:::