---
title: No15 - Generalizaciones
---

# Generalizaciones

**Fecha:** 10/06/2026

:::{iframe} https://www.youtube.com/embed/joZAfB1-1Ns
:width: 100%
:::

## SINDy: Identificación sparse de dinámicas no lineales

La idea de **SINDy** (*Sparse Identification of Nonlinear Dynamics*) es descubrir las
ecuaciones que gobiernan un sistema dinámico a partir de datos. Suponemos que el estado
$u(t) \in \mathbb{R}^n$ evoluciona según

```{math}
\frac{du}{dt} = f(u, \theta),
```

y postulamos que cada componente de $f$ se puede escribir como una **combinación lineal de
funciones candidatas** tomadas de un diccionario $\{1,\, u_1,\, u_2,\, u_3,\, u_1 u_2,\, u_1 u_3,\, \ldots\}$:

```{math}
f(u, \theta) =
\begin{bmatrix}
\theta_{11} + \theta_{12}\, u_1 + \theta_{13}\, u_2 + \theta_{14}\, u_3 + \theta_{15}\, u_1 u_2 + \theta_{16}\, u_1 u_3 + \cdots \\[4pt]
\theta_{21} + \theta_{22}\, u_1 + \theta_{23}\, u_2 + \cdots \\[4pt]
\vdots
\end{bmatrix}.
```

Los coeficientes forman una matriz

```{math}
\theta \in \mathbb{R}^{n \times M},
```

donde $M$ es el **tamaño del diccionario de funciones** (la cantidad de términos candidatos).

Notar que $f$ es en general **no lineal en $u$** (el diccionario contiene productos como
$u_1 u_2$, y podría contener $\sin(u)$, $u^3$, etc.), pero es **lineal en los coeficientes
$\theta$**. Esta estructura se puede escribir en forma compacta: si llamamos

```{math}
\Theta(u) = \begin{bmatrix} 1 & u_1 & u_2 & u_3 & u_1 u_2 & u_1 u_3 & \cdots \end{bmatrix}
\in \mathbb{R}^{1 \times M}
```

al vector fila con el diccionario evaluado en $u$, y $\Xi \in \mathbb{R}^{M \times n}$ a la
matriz de coeficientes (cada columna de $\Xi$ es una fila de la matriz $\theta$ de arriba),
entonces la ecuación de la foto es exactamente

```{math}
\frac{du}{dt} \approx \Theta(u)\, \Xi .
```

Es decir: la fila $k$ de $f$ en la forma expandida,
$\theta_{k1} + \theta_{k2} u_1 + \theta_{k3} u_2 + \cdots$, es el producto del vector de
funciones $\Theta(u)$ con la columna $k$ de $\Xi$. No estamos asumiendo nada nuevo al escribir
$\Theta(u)\Xi$: es la misma combinación lineal, reorganizada como producto matricial para que
el problema de encontrar los coeficientes sea una **regresión lineal**.

### Regresión sparse

Dado un muestreo de la trayectoria en tiempos $t_i$, evaluamos el diccionario sobre los
estados estimados $\hat{u}(t_i)$ para construir la **biblioteca** $\Theta(\hat{u}(t_i))$, y
buscamos los coeficientes $\Xi$ que reconstruyen la derivada observada. Para favorecer una
descripción **sparse** (que solo unos pocos términos del diccionario sobrevivan) agregamos una
penalización $\ell_1$:

```{math}
\min_{\Xi} \; \sum_i \left\| \frac{d\hat{u}}{dt}(t_i) - \Theta(\hat{u}(t_i))\, \Xi \right\|_2 + \lambda \, \| \Xi \|_1 .
```

:::{note} Sobre la notación
En la forma expandida de arriba, la matriz de coeficientes a aprender se anota $\theta$. En el
problema de regresión, $\Theta(\hat{u})$ denota la **biblioteca** de funciones evaluada en los
datos y $\Xi$ es la matriz de coeficientes sparse (el mismo rol que $\theta$, traspuesta). El
término $\lambda \| \Xi \|_1$ controla cuántos términos del diccionario quedan activos.
:::
