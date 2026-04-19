---
title: Clase No3 - UDEs
---

# Universal Differential Equations (UDEs)

**Fecha:** 22/04/2026

:::{iframe} https://www.youtube.com/embed/zkoWO0nnc3s
:width: 100%
:::


## Cosas que escribió Fran

Vimos el modelo de Lotka-Volterra, que tiene cuatro parámetros: $\alpha$, $\beta$, $\delta$ y $\gamma$. Según los valores de estos parámetros y dada una condición inicial, el sistema genera distintas trayectorias. A estas trayectorias se les puede agregar ruido gaussiano para representar los datos que observamos en la realidad.

Acá hacemos una distinción importante:

- El **estado del sistema** es la solución de las ecuaciones diferenciales. En la mayoría de los casos, este estado es desconocido.
- En la realidad observamos solo algunos puntos del sistema y no conocemos la dinámica subyacente. Además, las observaciones suelen tener ruido.

En muchos casos, el modelo observacional se construye a partir del modelo del estado. Por lo tanto, el modelo observacional dependerá tanto de la dinámica del sistema como del ruido de medición.

Hasta ahora, podemos resumir la situación de la siguiente manera:

$$
\theta \longrightarrow x(t;\theta) \longrightarrow y_i = x_{\text{obs}}(t_i;\theta) + \epsilon_i
$$

donde $\theta$ representa los parámetros del sistema, $x(t;\theta)$ el estado y $\epsilon_i$ el ruido observacional.

En este tipo de modelos, la evolución del estado no depende del ruido: el ruido aparece únicamente en las observaciones. En general, la intensidad del ruido no se conoce y puede terminar siendo otro parámetro a estimar.

Una forma simple de ajustar el sistema a los datos es mediante cuadrados mínimos:

$$
\min_{\theta} L(\theta) = \sum_i \left\| y_i - x(t_i;\theta) \right\|_2^2 .
$$

En el caso de Lotka-Volterra tenemos cuatro parámetros. Al graficar el valor de la función de pérdida en un mapa 2D, eligiendo alguna combinación de dos parámetros, ya se pueden observar problemas como la presencia de mínimos locales.

## NODEs

La idea de las redes neuronales para ecuaciones diferenciales ordinarias, o **Neural ODEs**, es reemplazar la función que describe la dinámica del sistema por una red neuronal:

$$
\frac{du}{dt} = f(u;\theta,t) \quad \Longrightarrow \quad \frac{du}{dt} = NN_{\theta}(u).
$$

¿Por qué hacer esto?

Porque las redes neuronales son aproximadores universales. Entonces, dados ciertos datos observacionales, la red puede aproximar el comportamiento del sistema incluso sin conocer explícitamente la ecuación que rige el fenómeno.

:::{note} Obs. 1: Sistemas autónomos
Un sistema autónomo es aquel en el que el tiempo no aparece explícitamente en las ecuaciones:

$$
\frac{du}{dt} = f(u;\theta).
$$

Si el sistema depende explícitamente del tiempo,

$$
\frac{du}{dt} = f(u,t;\theta),
$$

se puede transformar en un sistema autónomo agregando el tiempo como una variable más:

$$
\tilde{u} = (u,t) \in \mathbb{R}^{n+1},
$$

de modo que

$$
\frac{d\tilde{u}}{dt}
=
\begin{pmatrix}
f(u,t;\theta) \\
1
\end{pmatrix}.
$$
:::

Una red neuronal está compuesta, en general, por tres partes: una capa de entrada, una o más capas ocultas y una capa de salida. Cada capa tiene una cierta cantidad de neuronas, que depende del problema que se quiera resolver.

Las neuronas de distintas capas están conectadas mediante pesos y sesgos. Sin embargo, el poder de las redes neuronales aparece al introducir no linealidades mediante funciones de activación. Algunas de las más comunes son:

$$
\text{ReLU}(x) = \max(0,x),
$$

y

$$
\sigma(x) = \frac{1}{1+e^{-x}}.
$$

En nuestro caso, la red recibe como entrada el vector de estado $u$ y devuelve una aproximación de la función dinámica $f(u;\theta)$.

:::{note} Obs. 2: Construcción de NODEs

En las NODEs, los parámetros $\theta$ pasan a ser los pesos y sesgos de la red neuronal:

$$
\theta = \{W_i,b_i\}.
$$
:::


Volviendo al ejemplo de Lotka-Volterra, supongamos que conocemos una parte de la dinámica: sabemos que, en ausencia de interacción, la población de conejos crece y la de lobos decrece. Sin embargo, no sabemos exactamente cómo interactúan ambas poblaciones.

En ese caso, en lugar de reemplazar toda la dinámica por una red neuronal, podemos conservar la parte conocida del modelo y reemplazar solo la parte desconocida:

$$
\frac{dx}{dt} = \alpha x + NN_1(x,y),
$$

$$
\frac{dy}{dt} = -\beta y + NN_2(x,y).
$$

En este caso, los parámetros a ajustar son tanto los parámetros conocidos del modelo como los parámetros de las redes neuronales:

$$
\theta = \left[\alpha,\beta,W_1,\dots,W_n,b_1,\dots,b_n\right].
$$



% Pequeño comentario sobre NODEs

% Observación sobre sistemas autónomos

## Ejemplo: continuación

% Comentarios y continuación sobre el sistema de lotka-volterra. Minimos locales...


## Redes neuronales

## Calibración con respecto a la condición inicial



