---
title: No8 - PINNs
---

# Physics-Informend Neural Networks (PINNs)

**Fecha:** 06/05/2026

:::{iframe} https://www.youtube.com/embed/hK3IgpZeTIc
:width: 100%
:::

# Optimización con restricciones y dualidad lagrangiana

Muchos de los problemas que hasta ahora trabajamos parecieran tratar de minimizar la función de costo **sin restricciones**; por ejemplo, cuadrados mínimos:

$$
\min_{\theta} \mathcal{L} (\theta,y)=\min_{\theta}\sum_{i=1}^{N} \left\|\left\| y_i - x(t_i,\theta) \right\| \right\|_2^2 \qquad (1)
$$

donde $\mathcal{L}$ es la función de costo.

Podemos reescribir esto como un problema **con restricciones**, dejando a $x$ libre pero exigiendo que satisfaga una ecuación diferencial:

$$\min_{\theta,x}
\sum_{i=1}^{N}
\left\|\left\| y_i - x(t_i) \right\| \right\|_2^2
\quad
\text{sujeto a}
\quad
\begin{cases}
\dfrac{dx}{dt} = f(x,t,\theta) \\
x(t_0)=x_0
\end{cases}$$

Es decir, estamos convirtiendo

$$
\min_{\theta} f(x(\theta))
$$

en

$$
\min_{\theta,x} f(x,\theta)
\quad
\text{sujeto a}
\quad
G(x,\theta)=0
$$

Si uno puede invertir $G(x,\theta)$ para obtener $x=x(\theta)$ (ya sea analítica o numéricamente), volvemos al problema sin restricciones.

En este caso,

$$
G(x,\theta)=
\begin{bmatrix}
\dfrac{du}{dt} - f(u,t,\theta) \\
u(t_0)-u_0
\end{bmatrix}
=0
$$

> Esto se hace con el solver numérico.

---

# Forma general del problema de optimización

La forma más general del problema de optimización puede escribirse como:

$$
\min_{\theta} f(\theta)
$$

sujeto a

$$
\begin{cases}
g(\theta)=0 \\
h(\theta)\leq 0
\end{cases}
$$

La diferencia entre resolver estos problemas con y sin restricciones es debido a la dualidad lagrangiana.

---

# Dualidad lagrangiana

La dualidad lagrangiana toma un problema con restricciones y lo transforma en uno sin restricciones mediante mediante el metodo de los multiplicadores de Lagrange.

Se define el lagrangiano:

$$ \mathcal{L}(\theta,\lambda,\nu)
= f(\theta)
+
\lambda g(\theta)
+
\nu h(\theta)$$

donde:

- $\lambda$: multiplicadores de Lagrange asociados a restricciones de igualdad.
- $\nu$: multiplicadores asociados a restricciones de desigualdad.

---

Imaginemos un problema sin restricciones $h$. La dualidad lagrangiana nos dice, bajo ciertas hipótesis (dualidad fuerte) ({cite}`boyd2004convex`):

$$\max_{\lambda}\min_{\theta}\mathcal{L}(\theta,\lambda)=\min_{\theta}\max_{\lambda}\mathcal{L}(\theta,\lambda)$$

# Ejemplo: Lasso

Consideremos el problema de optimización:

$$\min_{\theta}
\left\| \left\| y - x \theta \right\| \right\|_2^2
+
\lambda \left\| \left\|  \theta \right\| \right\|_1$$

donde:

- el primer término corresponde al error de ajuste,
- el segundo término penaliza la complejidad del modelo (en este caso, prefiere soluciones esparsas).

Recordemos que

$$\left\| \left\|  \theta \right\| \right\|_1=\sum_{i=1}^{p} \left\| \left\| \theta_i \right\| \right\|$$

Entonces puede mostrarse que el $\theta^\ast$ que minimiza:

$$\theta^\ast=
\arg\min_{\theta}
\left\| \left\| y- x\theta \right\| \right\|_2^2
+
\lambda \left\| \left\| \theta  \right\| \right\|_1$$

puede reinterpretarse como

$$\theta^\ast = \arg\min_{\theta} \left\| \left\| y-x\theta \right\| \right\|_2^2 \quad \text{sujeto a} \quad \left\| \left\| \theta \right\| \right\|_1 \leq C(\lambda)$$

Por lo cual, la dualidad lagrangiana nos asegura que un problema de minimización de un lagrangiano puede reescribirse como un problema de optimización con restricciones.

Notar que si nuestro $\lambda \rightarrow 0$ entonces la constante $C \rightarrow \infty$.
Esto puede interpretarse como que si no restringimos demasiado nuestras soluciones a que cumplan nuestros vínculos, entonces habrá más libertad para los valores $\theta$ que se puedan obtener, logrando que $C$ sea muy grande; análogo ocurre para el caso extremo opuesto.
Es decir, notamos que la dependencia entre $\lambda$ y $C$ es inversa.

---

# Problema relajado

Hay posibilidad de "relajar" la función de costo con un factor $\varepsilon$. Por ejemplo, tomemos el problema (1) y escribámoslo como:

$$\min_{\theta,x} \mathcal{L}(\theta,x)= \qquad (3)$$

sujeto a

$$\left\| \left\| \frac{dx}{dt} - f(x,t,\theta) \right\| \right\|_2 \leq \varepsilon \qquad \forall t$$

y

$$\left\| \left\| x(t_0)-u_0 \right\| \right\|_2 \leq \varepsilon$$

donde recuperamos el problema (1) si $\varepsilon=0$ (es decir, si la ecuación diferencial se cumple exactamente). Una pregunta natural sería ¿de qué lagrangiano proviene esta restricción? Por dualidad lagrangiana, puede mostrarse que:

$$\min_{\theta,x} \mathcal{L}(\theta,x)+
\lambda_{\varepsilon}
\int_{t_0}^{t_1}
\left\| \left\|
\frac{dx}{dt}
-f(x,t,\theta) \right\| \right\|_2^2 dt \qquad (4)$$

El segundo término actúa como un término de regularización con derivadas (en el caso clásico de estadística esto corresponde a un “profiling”). Asimismo, vemos que si $\varepsilon \rightarrow 0$ entonces $\lambda_\varepsilon \rightarrow \infty$ y viceversa. Esta forma tiene la misma (o similar) que resolver mediante ["splines" o "smooth splines"](./clase6.md), solo que allí se utilizan las derivadas segundas.

La ecuación (4) constituye el punto de partida para una PINN (*Physics-Informed Neural Network*).


# Physics-Informed Neural Networks (PINNs)

Las PINNs fueron introducidas recientemente por Raissi, Perdikaris y Karniadakis en el siguiente trabajo: {cite}`Raissi_Perdikaris_Karniadakis_2019`.

## Caso ODE

La idea es agarrar el "profiling" y nuestras incógnitas a optimizar ($x(t)$) usando una red neuronal. Es decir, cambiamos nuestra $x(t)$ por:

$$
x_\beta(t)
$$

donde $\beta$ representa el conjunto de parámetros de la red neuronal:

$$
\beta = [W_1,\dots,W_n,b_1,\dots,b_n]
$$

por lo cual, con esto, lo ponemos en la ecuación diferencial y optimizamos sobre los parámetros de la red neuronal. Por ejemplo, una red neuronal puede escribirse como:

$$
x(t) = \sigma \left( W_3 \sigma \left( W_2 \sigma(W_1 t)+b_2\right) +b_3 \right) $$

---

Para construir la función de costo necesitamos:

1. Poder evaluar $x_\beta(t)$ para todo $t$.
2. Poder evaluar

$$
\frac{dx_\beta}{dt}\Big|_{t=s}
$$

para cualquier $s$.

---

# Modo 1: PINN forward/directo

En este caso, $\theta$ permanece fijo y no se optimiza.

Una vez obtenida $x(t)$, la introducimos en el segundo término de la ecuación (4).

Tomamos puntos de prueba:

$$
t_0 < z_1 < z_2 < \dots < z_k  < \dots < t_1
$$

y definimos la función de costo:

$$
\min_{\beta} \left\| \left\| x(t_0)-x_0\right\| \right\|_2^2 + \tilde{\lambda} \sum_{k=1}^{K} \left\| \left\| \left( \frac{dx_\beta}{dt} - f(x,t,\theta) \right)_{t=x_k} \right\| \right\|^2
$$

La suma actúa como una aproximación de la integral.

---

Esto es equivalente a utilizar una red neuronal como *solver* numérico.

---

# Ejemplo: ecuación del calor

Consideremos:

$$
x \in \mathbb{R}^n,
\qquad
t \in \mathbb{R},
\qquad
(n=1,2,3)
$$

La ecuación del calor es:

$$ \frac{\partial u}{\partial t} - D \nabla^2 u = 0
$$

donde $D$ es la difusividad y $\nabla^2 = \frac{\partial}{\partial x} + \frac{\partial}{\partial y} + \frac{\partial}{\partial z}$ es el laplaciano.

La ecuación debe satisfacerse para:

$$
\forall x\in\Omega,
\qquad
\forall t\in[t_0,t_1]
$$

---

con la condición de borde:

$$
u(x,t)=u_B(x,t)
\qquad
\forall x\in\partial\Omega
$$

---

y su condición inicial:

$$
u(x,t_0)=u_0(x)
$$

---

Vamos a considerar una parametrización para el borde espacial, otra para la condición inicial y finalmente para los puntos del interior. Podemos visualizar esto en el siguiente diagrama 
:::{figure} ../images/clase8.jpeg
:width: 100%
:align: center

Las cruces simbolizan la parametrización del borde espacial, los círculos el borde temporal, y los triángulos los puntos del interior.
:::

## Función de costo total

La función de costo a minimizar sobre los parámetros $\beta$ es:

$$\min_{\beta} \; \lambda_1 \sum_{i=1}^{K_1} \left\| \left\|  u_\beta(t_0, x_i^I) - u_0(x_i^I) \right\| \right\|^2 + \lambda_2 \sum_{j=1}^{K_2} \left\| \left\| u_\beta(t_j^B, x_j^B) - u_B(t_j^B, x_j^B) \right\| \right\|^2 + \lambda_3 \sum_{m=1}^{K_3} \left\| \left\| \mathcal{G}[u_\beta] \Big|_{t_M, x_M} \right\| \right\|_2^2$$

$$\equiv \mathcal{L}_{\text{inicial}} + \mathcal{L}_{\text{borde}} + \mathcal{L}_{\text{físico}}$$

> - Los **círculos** representan $u_0(x_i^I)$, con $I$ = Inicial.  
> - Las **cruces** a $u_B(t_j^B, x_j^B)$, y $B$ = Borde.
> - Los **triángulos** representan a $\mathcal{G}[u_\beta] \Big|_{t_M, x_M}$, donde $M$ = puntos de colocación.

---

## Arquitectura de la red y funciones de costo

La red neuronal recibe como entradas $t, x_1, \ldots, x_n$ y produce una salida escalar $u_\beta (t,x)$. A partir de esa salida se computan tres cantidades:

| Operación | Resultado | Función de costo asociada |
|---|---|---|
| $\text{Id}$ | $u_\beta$ | $\mathcal{L}_{\text{inicial}}$ |
| $\frac{\partial}{\partial t}$ | $\dfrac{\partial u_\beta}{\partial t}$ | $\mathcal{L}_{\text{borde}}$ |
| $\nabla$ | $\nabla u_\beta$ | $\mathcal{L}_{\text{físico}}$ |

Las tres funciones de costo se combinan en:

$$\mathcal{L}_{\text{TOTAL}} = \mathcal{L}_{\text{inicial}} + \mathcal{L}_{\text{borde}} + \mathcal{L}_{\text{físico}}$$

### Objetivo del modo Forward

$$\implies \min_{\beta} \, \mathcal{L}_{\text{TOT}}(\beta) \simeq 0 \quad \text{y} \quad u_\beta \text{solución a la ecuación}$$

> **Hasta aquí todo es modo Forward.**

---

## PINNs Modo 2: Inverso

¡Este es el uso verdadero de una PINN! El modo "forward" es mas de juguete que otra cosa.

### Problema

Dados los datos observados $\left[u^{\text{obs}}(t_i, x_i)\right]_{i=1}^{N}$, se busca **recuperar** $D = D(x)$: la difusividad como función de $x$.

### Idea

Se desarrolla $D$ como una **Red Neuronal** con parámetros $\beta_2$, y se considera una función de costo empírica adicional:

$$\mathcal{L}_{\text{emp}} = \lambda \sum_{i=1}^{N} \left\| \left\| u_{\beta_1}(x_i, t_i) - u^{\text{obs}}(t_i, x_i) \right\| \right\|_2^2$$

### Ahora la arquitectura involucra dos redes

Se tienen **dos redes** con parámetros distintos:

- **Red $\beta_1$**: aproxima la solución $u_{\beta_1}(t, x)$, igual que en el modo Forward.
- **Red $\beta_2$**: aproxima la difusividad $D_{\beta_2}(x) \to \mathbb{R}$ (salida escalar, solo depende de $x$).

| Operación sobre $u_{\beta_1}$ | Resultado | Función de costo asociada |
|---|---|---|
| $\text{Id}$ | $u_{\beta_1}$ | $\mathcal{L}_{\text{inicial, borde, empírica}}$ |
| $\frac{\partial}{\partial t}$ | $\dfrac{\partial u_{\beta_1}}{\partial t}$ | $\mathcal{L}_{\text{física}}$ |
| $\nabla$ | $\nabla u_{\beta_1}$ | $\mathcal{L}_{\text{física}}$ |

| Operación sobre $D_{\beta_1}$ | Resultado | Función de costo asociada |
|---|---|---|
| $\text{Id}$ | $D_{\beta_1}$ | $\mathcal{L}_{\text{física}}$ |


La función de costo empírica también se incorpora a $\mathcal{L}_{\text{TOTAL}}$.

### Parámetros a optimizar

$$\beta = [\beta_1, \beta_2]$$

Se minimiza la función de costo total conjuntamente sobre $\beta_1$ y $\beta_2$, de modo que al final $u_{\beta_1}$ aproxima la solución **y** $D_{\beta_2}$ recupera la difusividad desconocida.
