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

## Ecuaciones diferenciales con derivadas de mayor orden

Muchos problemas físicos aparecen originalmente con derivadas de orden mayor que uno, pero en general pueden reescribirse como sistemas de primer orden.

Por ejemplo, si tenemos el oscilador armónico

$$
\frac{d^2x}{dt^2} + \omega^2 x = F(x,t),
$$

y definimos

$$
v = \frac{dx}{dt},
$$

podemos escribir

$$
\frac{d}{dt}
\begin{bmatrix}
x \\
v
\end{bmatrix}
$$

=

$$
\begin{bmatrix}
v \\
-\omega^2 x + F(x,t)
\end{bmatrix}
$$

=

$$
\begin{bmatrix}
0 & 1 \\
-\omega^2 & 0
\end{bmatrix}
\begin{bmatrix}
x \\
v
\end{bmatrix}
+
\begin{bmatrix}
0 \\
F(x,t)
\end{bmatrix}
$$

En la práctica, gran parte de la teoría y del cómputo numérico para ODEs se apoya en este tipo de formulaciones de primer orden.

## Ecuaciones diferenciales en derivadas parciales

Hay problemas en derivadas parciales que también pueden reducirse a sistemas de ODEs luego de discretizar el espacio.

Por ejemplo, si

$$
\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2},
\qquad
u = u(t,x),
$$

entonces una aproximación por diferencias finitas da

$$
\frac{\partial^2 u}{\partial x^2}(x)
\approx
\frac{u(x+h)-2u(x)+u(x-h)}{h^2}.
$$

De esta forma, para los valores nodales $u_0, u_1,\dots, u_n$, obtenemos un sistema de ODEs:

$$
\frac{d}{dt}
\begin{bmatrix}
u_0 \\
u_1 \\
\vdots \\
u_n
\end{bmatrix}
$$

=

$$
\frac{D}{h^2}
\begin{bmatrix}
-2 & 1 & 0 & \cdots & 0 \\
1 & -2 & 1 & \cdots & 0 \\
0 & 1 & -2 & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & 1 \\
0 & \cdots & 0 & 1 & -2
\end{bmatrix}
\begin{bmatrix}
u_0 \\
u_1 \\
\vdots \\
u_n
\end{bmatrix}
$$

Este procedimiento es una versión del método de las líneas: discretizamos las derivadas espaciales y dejamos el tiempo como variable continua. No es trivial fijar las condiciones de borde para este tipo de problemas, y es fundamental.

Nos va a interesar, a partir de datos, estimar el parámetro $\theta$ de $f(u,t,\theta) = \frac{du}{dt}$
Notación: usamos $u$ de "unknown", sino $x$ e $y$ para variables concretas.

### Ejemplos

Un ejemplo sencillo de dinámica de poblaciones consiste en modelar

-   $x(t)$: población de conejos;
-   $y(t)$: población de lobos.

Sin interacción, podríamos escribir

$$
\frac{dx}{dt} = \alpha x, \qquad \alpha > 0,
$$

y

$$
\frac{dy}{dt} = -\beta y, \qquad \beta > 0.
$$

La interpretación es que, sin depredadores, los conejos crecen; y sin
presas, los lobos decrecen.

Si ahora incorporamos interacción, aparece el sistema de
Lotka-Volterra, un "caballo de batalla" para modelar comportamientos entre especiaes en biología:

$$
\begin{cases}
\dfrac{dx}{dt} = \alpha x - \delta x y, \\
\dfrac{dy}{dt} = -\beta y + \gamma x y,
\end{cases}
\qquad
\alpha,\beta,\gamma,\delta > 0.
$$

El estado inicial se resume como

$$
\begin{bmatrix}
x \\
y
\end{bmatrix}(t_0)
$$

=

$$
\begin{bmatrix}
x_0 \\
y_0
\end{bmatrix},
$$

y el vector de parámetros puede escribirse como

$$
\theta =
\begin{bmatrix}
\alpha \\
\beta \\
\gamma \\
\delta
\end{bmatrix}
\in \mathbb{R}^4.
$$

Distintas elecciones de $\theta$ producen trayectorias oscilantes
cualitativamente distintas para $x(t;\theta)$ e $y(t;\theta)$. Estos parámetros involucran no-linealidades en el modelo.

## Inferencia estadística

Ahora, vamos a distinguir entre:

1.  el **modelo de estado**, que describe la dinámica del sistema;
2.  el **modelo observacional**, que conecta esa dinámica con los datos
    medidos.

### Modelo de estado

En el ejemplo anterior, el modelo de estado está dado por la ODE que
gobierna

$$
\begin{bmatrix}
x(t;\theta) \\
y(t;\theta)
\end{bmatrix}.
$$

Este modelo representa el mecanismo subyacente del sistema, aunque dicho
mecanismo no se observe de manera perfecta.

### Modelo observacional

Si observamos la variable $x$ en tiempos $t_1,\dots,t_N$, una forma
simple de escribir el modelo observacional es

$$
X_i^{\mathrm{obs}} = X(t_i;\theta) + \varepsilon_i.
$$

Si además se observara la población de lobos, tendría sentido considerar
también

$$
Y_i^{\mathrm{obs}} = Y(t_i;\theta) + \eta_i.
$$

Aquí $X_i^{\mathrm{obs}}$ y $Y_i^{\mathrm{obs}}$ son observaciones
discretas, mientras que $\varepsilon_i$ y $\eta_i$ representan ruido o
error observacional.

Una hipótesis frecuente es asumir

$$
\varepsilon_i \stackrel{\mathrm{iid}}{\sim} N(0,\sigma^2),
$$

aunque en la práctica esta es una simplificación: los errores pueden no
ser gaussianos, no ser independientes, o presentar correlación temporal.

La diferencia entre observación y predicción del modelo es el residuo.

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
