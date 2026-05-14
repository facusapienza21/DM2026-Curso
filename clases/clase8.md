---
title: No8 - PINNs
---

# Physics-Informend Neural Networks (PINNs)

**Fecha:** 06/05/2026

:::{iframe} https://www.youtube.com/embed/hK3IgpZeTIc
:width: 100%
:::

# Optimización con restricciones y dualidad lagrangiana

Este es un problema de optimización **sin restricciones**:

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

Si uno puede invertir $G(x,\theta)$ para obtener $x=x(\theta)$, volvemos al problema sin restricciones.

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

---

# Dualidad lagrangiana

La dualidad lagrangiana toma un problema con restricciones y lo transforma en uno sin restricciones mediante multiplicadores de Lagrange.

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

Imaginemos un problema sin restricciones $h$. La dualidad lagrangiana nos dice:

$$\min_{\lambda}\max_{\theta}\mathcal{L}(\theta,\lambda)=\max_{\theta}\min_{\lambda}\mathcal{L}(\theta,\lambda)$$

# Ejemplo: LASSO

Consideremos el problema de optimización:

$$\min_{\theta}
\left\| \left\| y - x \theta \right\| \right\|_2^2
+
\lambda \left\| \left\|  \theta \right\| \right\|_1$$

donde:

- el primer término corresponde al error de ajuste,
- el segundo término penaliza la complejidad del modelo.

Recordemos que

$$\left\| \left\|  \theta \right\| \right\|_1=\sum_{i=1}^{p} |\theta_i|$$

Entonces,

$$\theta^\ast=
\arg\min_{\theta}
\left\| \left\| y- x\theta \right\| \right\|_2^2
+
\lambda \left\| \left\| \theta  \right\| \right\|_1$$

puede reinterpretarse como

$$\theta^\ast=
\arg\min_{\theta} \left\{
\begin{aligned}
& \left\| \left\| y-x\theta \right\| \right\|_2^2 \\
&\text{sujeto a } \left\| \left\| \theta  \right\| \right\|_1 \leq C(\lambda)
\end{aligned} \right.$$
---

# Problema relajado

Tomemos el problema (1) y escribámoslo como:

$$\min_{\theta,x} \mathcal{L}(\theta,x)=f(\theta,x)
\qquad (3)$$

sujeto a

$$\left\| \left\|
\frac{dx}{dt}
-
f(x,t,\theta)
\right\| \right\|_2
=
\varepsilon
\qquad \forall t$$

y

$$\left\| \left\| x(t_0)-u_0 \right\| \right\|_2 \leq \varepsilon$$

El problema (1) se recupera cuando $\varepsilon=0$.

---

# Dualidad lagrangiana

Aplicando dualidad lagrangiana:

$$\min_{\theta,x} \mathcal{L}(\theta,x)+
\lambda_{\varepsilon}
\int_{t_0}^{t_1}
\left\| \left\|
\frac{dx}{dt}
-f(x,t,\theta) \right\| \right\|_2^2 dt \qquad (4)$$

El segundo término actúa como un término de regularización con derivadas.

> En el caso clásico de estadística esto corresponde a un “profiling”.

---

Ahora no resolvemos explícitamente ninguna ecuación diferencial, sino que optimizamos directamente sobre $x$.

Este problema puede resolverse usando splines.

La ecuación (4) constituye el punto de partida para una PINN (*Physics-Informed Neural Network*).

# Physics-Informed Neural Networks (PINNs)

## Caso ODE

La idea es escribir $x(t)$ como una red neuronal:

$$
x_\beta(t)
$$

donde $\beta$ representa el conjunto de parámetros de la red neuronal:

$$
\beta = [W_1,\dots,W_n,b_1,\dots,b_n]
$$

Por ejemplo, una red neuronal puede escribirse como:

$$
x(t)
=
\sigma\!\left(
W_3\,
\sigma\!\left(
W_2\,
\sigma(W_1 t)+b_2
\right)
+b_3
\right)
$$

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
\min_{\beta}
\left\| \left\| x(t_0)-x_0\right\| \right\|_2^2
+
\tilde{\lambda}
\sum_{k=1}^{K}
\left\| \left\|
\left(
\frac{dx_\beta}{dt}
-
f(x,t,\theta)
\right)_{t=x_k}
\right\| \right\|^2
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

$$
\frac{\partial u}{\partial t}
-
D\nabla^2 u
=
0
$$

donde $D$ es la difusividad.

La ecuación debe satisfacerse para:

$$
\forall x\in\Omega,
\qquad
\forall t\in[t_0,t_1]
$$

---

## Condiciones de borde

$$
u(x,t)=u_B(x,t)
\qquad
\forall x\in\partial\Omega
$$

---

## Condición inicial

$$
u(x,t_0)=u_0(x)
$$

---

Vamos a considerar una parametrización para el borde espacial, otra para la condición inicial y finalmente para los puntos del interior. Podemos visualizar esto en el siguiente diagrama 
![Las cruces simbolizan la parametrización del borde espacial, los círculos el borde temporal, y los triangulos los puntos del interior.](images/clase8.jpeg)

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

### Problema

Dados los datos observados $\{u^{\text{obs}}(t_i, x_i)\}_{i=1}^{N}$, se busca **recuperar** $D = D(x)$: la difusividad como función de $x$.

### Idea

Se desarrolla $D$ como una **Red Neuronal** con parámetros $\beta_2$, y se considera una función de costo empírica adicional:

$$\mathcal{L}_{\text{emp}} = \lambda \sum_{i=1}^{N} \left\| \left\| u_{\beta_1}(x_i, t_i) - u^{\text{obs}}(t_i, x_i) \right\| \right\|_2^2$$

### Ahora la arquitectura involucra dos redes

Se tienen **dos redes** con parámetros distintos:

- **Red $\beta_1$**: aproxima la solución $u_{\beta_1}(t, x)$, igual que en el modo Forward.
- **Red $\beta_2$**: aproxima la difusividad $D_{\beta_2}(x) \to \mathbb{R}$ (salida escalar, solo depende de $x$).

| Operación sobre $u_{\beta_1}$ | Resultado | Función de costo asociada |
|---|---|---|
| $\text{Id}$ | $u_{\beta_1}$ | $\mathcal{L}_{\text{inicial}}$, $\mathcal{L}_{\text{borde}}$, $\mathcal{L}_{\text{empírica}}$ |
| $\frac{\partial}{\partial t}$ | $\dfrac{\partial u_{\beta_1}}{\partial t}$ | $\mathcal{L}_{\text{física}}$ |
| $\nabla$ | $\nabla u_{\beta_1}$ | $\mathcal{L}_{\text{físico}}$ (usa $D_{\beta_2}$) |

| Operación sobre $D_{\beta_1}$ | Resultado | Función de costo asociada |
|---|---|---|
| $\text{Id}$ | $D_{\beta_1}$ | $\mathcal{L}_{\text{física}}$ |


La función de costo empírica $\mathcal{L}_{\text{emp}}$ también se incorpora a $\mathcal{L}_{\text{TOTAL}}$.

### Parámetros a optimizar

$$\beta = [\beta_1, \beta_2]$$

Se minimiza $\mathcal{L}_{\text{TOTAL}}$ conjuntamente sobre $\beta_1$ y $\beta_2$, de modo que al final $u_{\beta_1}$ aproxima la solución **y** $D_{\beta_2}$ recupera la difusividad desconocida.
