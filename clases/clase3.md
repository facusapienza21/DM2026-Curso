---
title: Clase No3 - UDEs
---

# Universal Differential Equations (UDEs)

**Fecha:** 22/04/2026

:::{iframe} https://www.youtube.com/embed/zkoWO0nnc3s
:width: 100%
:::

Vimos el modelo de Lotka-Volterra, que tiene cuatro parámetros representados como $\theta$ = ($\alpha$, $\beta$, $\delta$ y $\gamma$).
Según los valores de estos parámetros y dada una condición inicial, el sistema genera distintas trayectorias.
A estas trayectorias se les puede agregar ruido gaussiano para representar los datos que observamos en la realidad.

A continuación tenemos dos implementaciones interactivas del modelo de Lotka-Volterra, una en Python y otra en Julia.
En ambos casos se pueden modificar los parámetros y explorar cómo cambian las trayectorias temporales y el retrato de fases.

- [Código interactivo en Python](./code/03_LV_forward_UDE/Lotka_volterra_interactivo.py)
- [Código interactivo en Julia](./code/03_LV_forward_UDE/lotka_volterra_ude.jl)

:::{figure} ./figures/no3_lv_interactivo_py.png
:width: 100%
:align: center

Interfaz interactiva en Python para explorar la dinámica del modelo de Lotka-Volterra variando sus parámetros.
:::

Acá hacemos una distinción importante:

- El **estado del sistema** es la solución de las ecuaciones diferenciales.
- En la mayoría de los casos, este estado es desconocido.
- En la realidad observamos solo algunos puntos del sistema y no conocemos la dinámica subyacente.
- Además, las observaciones suelen tener ruido.

En muchos casos, el modelo observacional se construye a partir del modelo del estado.
Por lo tanto, el modelo observacional dependerá tanto de la dinámica del sistema como del ruido de medición.

Hasta ahora, podemos resumir la situación de la siguiente manera:
$$
\theta \longrightarrow x(t;\theta) \longrightarrow y_i = x_{\text{obs}}(t_i;\theta) + \epsilon_i
$$
donde $\theta$ representa los parámetros del sistema, $x(t;\theta)$ el estado y $\epsilon_i$ el ruido observacional.

En este tipo de modelos, la evolución del estado no depende del ruido.
El ruido aparece únicamente en las observaciones.
En general, la intensidad del ruido no se conoce y puede terminar siendo otro parámetro a estimar.

Una estrategia para ajustar el sistema a los datos consiste en minimizar el cuadrado de los residuos entre las observaciones y la trayectoria predicha por el modelo.
Este procedimiento se conoce como **ajuste por trayectorias** o **cuadrados mínimos no lineales**.
Para una discusión más detallada de esta metodología, ver la {doc}`Clase N.º 2 <./clase2>`.
$$
\min_{\theta} L(\theta) = \sum_i \left\| y_i - x(t_i;\theta) \right\|_2^2 .
$$
En el caso de Lotka-Volterra tenemos cuatro parámetros.
Al graficar el valor de la función de pérdida en un mapa 2D, eligiendo alguna combinación de dos parámetros, ya se pueden observar problemas como la presencia de mínimos locales.

## NODEs

La idea de las redes neuronales para ecuaciones diferenciales ordinarias, o **Neural ODEs**, es reemplazar la función que describe la dinámica del sistema por una red neuronal:
$$
\frac{du}{dt} = f(u;\theta,t) \quad \Longrightarrow \quad \frac{du}{dt} = NN_{\theta}(u).
$$
¿Por qué hacer esto?

Porque las redes neuronales son aproximadores universales.
Entonces, dados ciertos datos observacionales, la red puede aproximar el comportamiento del sistema incluso sin conocer explícitamente la ecuación que rige el fenómeno.

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

Una red neuronal está compuesta, en general, por tres partes: una capa de entrada, una o más capas ocultas y una capa de salida.
Cada capa tiene una cierta cantidad de neuronas, que depende del problema que se quiera resolver.

Matemáticamente, una red neuronal feedforward puede escribirse como una composición de transformaciones afines y funciones de activación no lineales.
Si la entrada es $u \in \mathbb{R}^n$, una red con $L$ capas puede representarse como
$$
NN_\theta(u)=h_L \circ h_{L-1} \circ \cdots \circ h_1(u),
$$
donde cada capa tiene la forma
$$
h_k(z)=\sigma_k(W_k z+b_k).
$$
Aquí, $W_k$ y $b_k$ son los pesos y sesgos de la capa $k$, respectivamente.
Las funciones $\sigma_k$ introducen la no linealidad en cada capa.
El conjunto de parámetros de la red es entonces
$$
\theta=\{W_1,\dots,W_L,b_1,\dots,b_L\}.
$$
Algunas de las funciones de activación más comunes son
$$
\text{ReLU}(x)=\max(0,x),
$$
y
$$
\sigma(x)=\frac{1}{1+e^{-x}}.
$$
En nuestro caso, la red recibe como entrada el vector de estado $u$ y devuelve una aproximación de la función dinámica $f(u;\theta)$.

:::{note} Obs. 2: Construcción de NODEs
En las NODEs, los parámetros $\theta$ pasan a ser los pesos y sesgos de la red neuronal:
$$
\theta = \{W_i,b_i\}_{i=0}^K.
$$
:::

Volviendo al ejemplo de Lotka-Volterra, supongamos que conocemos una parte de la dinámica.
Sabemos que, en ausencia de interacción, la población de conejos crece y la de lobos decrece.
Sin embargo, no sabemos exactamente cómo interactúan ambas poblaciones.

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

## Calibración con respecto a la condición inicial

También vimos, al observar las curvas de MSE en distintas proyecciones del espacio de parámetros, la presencia de mínimos locales.
Esto nos indica que incluso en este problema canónico podemos caer en soluciones que no corresponden al mínimo global.

Si uno fija una condición inicial y aplica únicamente descenso por gradiente, es posible converger a uno de estos mínimos locales.
Para evitarlo, existen algoritmos de búsqueda global que permiten obtener una buena estimación inicial de los parámetros y así aumentar la probabilidad de alcanzar el mínimo absoluto.
