---
title: No6 - Función de costo
---

# Funciones de costo, modelo observacional, y más optimización

**Fecha:** 27/04/2026

:::{iframe} https://www.youtube.com/embed/insx5toVEus
:width: 100%
:::

En esta clase nos vamos a centrar en responder:

## Origen de la función de costo $\mathcal{L}$

La función de costo se deriva del modelo observacional que conecta la ecuación de estado de mi problema con los datos.

**Ejemplo 1: distribución Gaussiana**

En nuestro ejemplo característico de Lotka-Volterra (*Depredador-Presa* ver {doc}`Clase N.º 2 <./clase2>`) tenemos:
$$
\mathcal{L}(\theta,y) = \sum_{i=1}^{N} \left\| x(t_i; \theta) - y_i \right\|_{2}^{2} ,
$$
$$
y_i = x(t_i; \theta) + \varepsilon_i ,
$$
$$
\frac{dx}{dt} = f(x, t, \theta) ,
$$
donde $x(t_i; \theta)$ es la función de estado de mi sistema, que está descripta en este caso (*y en los que se va a enfocar este curso*) por una ecuación diferencial y un $\varepsilon_i$ que representa el ruido observacional.

Para Lotka-Volterra tenemos $\theta \in \mathbb{R}^{p}$ con $n = 2$ (dimensión del estado) y $p = 4$ (número de parámetros).

:::{note}
Cuando $p \gg 1$, este planteo nos lleva al régimen de las **NODE** (Neural ODEs), donde la dinámica está parametrizada por una red neuronal con muchos parámetros.
:::

Supongamos que estamos en una dimensión ($n=1$) y que el ruido $\varepsilon_i \sim \mathcal{N}(0,\sigma^{2})$ está caracterizado por una distribución gaussiana donde los $\varepsilon_i$ son independientes entre sí y de $x$ y están idénticamente distribuidos.

La probabilidad Gaussiana de observar $y_i$ dado $x_i$ y $\sigma$ de este modelo es:
$$
P(y_i|x_i,\sigma) = \frac{1}{\sqrt{2\pi}\sigma} \,e^{-\frac{(y_i-x_i)^{2}}{2\sigma^{2}}}
$$
y como las variables son independientes, la probabilidad de observar todos los puntos es:
$$
P(y_1,...,y_n)|x(t_i;\theta)) = \prod_{i=1}^{n}P(y_i|x_i,\sigma),
$$
donde $P(y_1,...,y_n)|x(t_i;\theta))$ es la **verosimilitud** y la llamaremos $L(\theta;y)$. 

Este va a ser nuestro modelo probabilístico para este ejemplo, que nos dice dadas nuestras distribuciones de probabilidad, cómo los datos desde $y_1$ hasta $y_n$ se generan aleatoriamente.

:::{tip}
Para facilitar las cuentas definimos:
$$
\ell(\theta;y) = \log L(\theta;y) = \sum_{i=1}^{N} \log P(y_i|x_i;\theta)
$$
Esto lo hacemos para sacarnos de encima los productos y tener todo descrito por sumas.
:::

Para poder hacer inferencia entre los parámetros del problema y nuestro modelo probabilístico vamos a utilizar el {term}`principio de máxima verosimilitud <Principio de Máxima Verosimilitud>`
$$
\ell(\theta;y_i) = -\sum_{i=1}^{N} \left( \frac{(y_i-x_i)^{2}}{2\sigma^{2}} + \log(\sqrt{2\pi}\,\sigma) \right),
$$
aplicamos el principio de máxima verosimilitud:
$$
\hat{\theta}_{MLE} = \arg\max_{\theta} \left[ -\sum_{i=1}^{N} \left( \frac{(y_i-x_i)^{2}}{2\sigma^{2}} + \log(\sqrt{2\pi}\,\sigma) \right) \right]
$$
Notemos que estamos maximizando sobre la variable $\theta$, entonces nos podemos "deshacer" del término $\log(\sqrt{2\pi}\,\sigma)$ ya que no depende de esta variable, y sacarlo no va a cambiar el resultado esperado, pero nos va a proporcionar una ecuación mucho más simple:
$$
\hat{\theta}_{MLE} = \arg\min_{\theta} \left[ \frac{1}{2\sigma^2} \sum_{i=1}^{N} (y_i-x_i)^{2} \right] = \arg\min_{\theta} \left[ \sum_{i=1}^{N} (y_i-x_i)^{2} \right]
$$
En conclusión, para este problema, el estimador de máxima verosimilitud es el que minimiza los residuos cuadráticos.

:::{important}
En muchos casos podemos mirar un problema de optimización como uno de máxima verosimilitud, mediante la siguiente identidad:
$$
\max_{\theta} L(\theta;y) = \min_{\theta} \left( -\log \mathcal{L}(\theta,y) \right)
$$
:::

**Ejemplo 2: distribución Laplaciana**

Asumimos $\varepsilon_i$ con una distribución de Laplace
$$
\varepsilon_i \sim \text{Lap}(0,b)
$$
Cuya probabilidad es:
$$
P(\varepsilon_i|b) = \frac{1}{2b}\, e^{-\frac{|\varepsilon_i|}{2b}} \quad \text{con } b>0
$$
Como en el ejemplo anterior, vamos a maximizar la variable $\theta$ para encontrar la función de costo.
$$
\hat{\theta}_{MLE} = \arg\max_{\theta} \left[ -\frac{1}{2b} \sum_{i=1}^{N} |y_i-x_i| \right] = \arg\min_{\theta} \left[ \sum_{i=1}^{N} |y_i-x_i| \right]
$$
Luego la función de costo es:
$$
\mathcal{L}(\theta,y) = \sum_{i=1}^{N} |y_i-x_i|
$$
:::{note}
Es una distribución que tiene colas más pesadas (los extremos decaen más lentamente) a diferencia de la Gaussiana, por eso se usa para hacer estadística más robusta.
:::

Hasta ahora vinimos haciendo máxima verosimilitud sólo sobre los parámetros $\theta$, en el siguiente ejemplo veremos qué pasa si $\sigma$ también es un parámetro.

**Ejemplo 3: distribución gaussiana con $\sigma_i \neq \text{cte}$**

Para una distribución gaussiana cuyo $\sigma_i$ ahora no es constante, tenemos:
$$
\mathcal{L}(\theta;y) = \sum_{i=1}^{N} \omega_i (y_i-x_i)^{2}
$$
Entonces nos queda una función de costo como una función de cuadrados mínimos pesados.

**Ejemplo 4: generalización de la distribución gaussiana**

Asumimos que los $\varepsilon_i$ están distribuidos de forma gaussiana, pero esta vez están correlacionados entre sí, esto quiere decir, que no son independientes entre sí.

Su matriz de covarianza $\Sigma_{ij}$ representa el valor medio $\mathbb{E} [\varepsilon_i,\varepsilon_j]$:
$$
\boldsymbol{\varepsilon} = \begin{bmatrix}
\varepsilon_1 \\
\vdots \\
\varepsilon_N
\end{bmatrix} \sim \mathcal{N}(\bar{0},\Sigma)
$$
Cuya probabilidad es:
$$
P(\boldsymbol{\varepsilon}|\Sigma) = \frac{1}{(2\pi)^{N/2}\,\left|\det(\Sigma)\right|^{1/2}}\, e^{-\frac{1}{2}\boldsymbol{\varepsilon}^{T}\Sigma^{-1} \boldsymbol{\varepsilon}}
$$
:::{note}
Esta matriz general $\Sigma_{ij}$ representa geométricamente una distribución gaussiana rotada, a esto se lo conoce como una distribución normal multivariada.
:::

Ahora la función de costo $\mathcal{L}(\theta)$ está dada de la siguiente forma:
$$
\mathcal{L}(\theta,y) = (x-y)^{T} \Sigma^{-1} (x-y) = \left\| y -x \right\|_{\Sigma}
$$
:::{note}
Esta norma está pesada en $\Sigma$ y se conoce como **norma de Mahalanobis**, notar que si sólo me queda la diagonal de esta matriz $\Sigma$ me devuelve la norma Euclídea.
:::

## Biyección entre $L(\theta,y)$ y $\mathcal{L}(\theta,y)$:

No siempre vamos a poder encontrar una biyección entre $L(\theta,y)$ y $\mathcal{L}(\theta,y)$, pero podemos hacer lo siguiente:

Dado $\mathcal{L}(\theta,y)$ queremos encontrar $L(\theta,y)$, para ello vamos a definirnos una función de probabilidad a la que llamaremos $L^{*}(\theta,y)$:
$$
L^{*}(\theta,y) =\frac{e^{-\mathcal{L}(\theta,y)}}{z(\theta)} 
$$
$$
z(\theta) = \int e^{-\mathcal{L}(\theta,y)}dy,
$$
donde $z(\theta)$ representa el factor de normalización, que puede depender del parámetro $\theta$, esto no nos pasaba en los ejemplos anteriores y como consecuencia, al recuperar la función de costo como veníamos haciendo, se nos va a agregar un término extra.
$$
\mathcal{L}^{*}(\theta,y) = -\log L^{*}(\theta,y) = \mathcal{L}(\theta,y) + \log z(\theta)
$$
definimos $R(\theta) = \log z(\theta)$ y lo llamaremos término de regularización.

Luego **el caso más general** de un problema de optimización va a tener esta forma:
$$
\mathcal{L}(\theta,y) =  \mathcal{L}_{\text{empírica}}(\theta,y) + R(\theta) ,
$$
donde $\mathcal{L}_{\text{empírica}}(\theta,y)$ depende tanto de los parámetros como de los datos y $R(\theta)$ sólo depende de los parámetros.

**Ejemplos:** vamos a ver distintas funciones de costo que suelen aparecer cotidianamente, con $y \in \mathbb{R}^{n}$, $x \in \mathbb{R}^{n \times p}$, $\theta \in \mathbb{R}^{p}$ y los $\lambda$ son hiperparámetros del problema.

:::{important}
Todos los ejemplos que vamos a mostrar tienen solución analítica exacta.
:::

**1) Regresión lineal Ridge** 
$$
\min_{\theta} \underbrace{\left\| y - x \theta \right\|^{2}_{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda \left\| \theta \right\|^{2}_{2}}_{\text{Término de Regularización}}
$$
El término de Regularización penaliza la norma dos del vector, y se lo llama **Ridge**, esto provoca que el parámetro $\theta$ no se mueva en demasía, es decir, que tienda a converger a cero, de manera tal que cuando ingresen nuevos datos en el programa la curva ajuste mejor.

**2) Regresión Lineal Lasso**
$$
\min_{\theta} \underbrace{\left\| y - x \theta \right\|^{2}_{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda \left\|\theta \right\|_{1}}_{\text{Término de Regularización}}
$$
El término de Regularización penaliza la norma uno del vector, y se lo llama **Lasso**, esto provoca esparcidad en las soluciones.

**3) Regresión lineal Elastic-Net**
$$
\min_{\theta} \underbrace{\left\| y - x \theta \right\|^{2}_{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda (\alpha \left\| \theta  \right\|_{1} + (1-\alpha) \left\| \theta \right\|^{2}_{2})}_{\text{Término de Regularización}} \ con \ \alpha \in [0,1]
$$
Combinación de los ejemplos 1 y 2.

**4) Smoothing Splines**

En el espacio que vamos a estar pensando es en el espacio funcional de dimensión infinita que contiene a todas las funciones que tienen hasta la segunda derivada continua.
$$
\min_{f} \underbrace{\sum_{i=1}^{N} (y_{i} - f(x_{i}))^{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda \int_{x_{0}}^{x_{1}}(f^{''}(x))^{2}dx}_{\text{Término de Regularización}}
$$
El término de Regularización penaliza la segunda derivada, lo que impone suavidad sobre las posibles soluciones.

Podemos observar qué pasa cuando variamos el $\lambda$:

- $\text{Si } \lambda \longrightarrow \infty \Rightarrow \text{regresión lineal}$

- $\text{Si } \lambda = 0 \Rightarrow \text{interpolación}$

## Introducción a la estadística Bayesiana:

Hasta ahora resolvimos nuestros problemas pensándolos con estadística frecuentista, y nos preguntamos cómo se relacionan estas ideas con la estadística Bayesiana (ver {doc}`Clase N.º 7 <./clase7>`).

En la estadística Bayesiana vamos a tener:

1) Verosimilitud, la probabilidad $\mathbb{P}(y|\theta)$ que nos dice cómo los datos están generados en función de los parámetros $\theta$.

2) $\mathbb{P}(\theta)$ que es una distribución de probabilidad sobre&nbsp;$\theta$, donde $\theta$ es **aleatorio**.

Utilizando la definición de la Probabilidad Condicional:
$$
\mathbb{P}_{post}(\theta,y) = \frac{\mathbb{P}(y,\theta) \mathbb{P}_{prior}(\theta)}{\mathbb{P}(y)}
$$
Si podemos calcular esta distribución, no solo vamos a obtener el $\theta$ que maximiza nuestro modelo, sino que también, nos va a dar una **noción de la incertidumbre** alrededor de ese $\theta$.

Buscando quien maximiza la distribución $\mathbb{P}_{post}(\theta,y)$, podemos recuperar la solución encontrada con la estadística frecuentista.
$$
\theta_{MAP} = \arg\max_{\theta}\mathbb{P}_{post}(\theta,y) = \arg\min_{\theta} \left[ -\underbrace{\log \mathbb{P}(y,\theta)}_{\ell(\theta,y)} - \underbrace{\log \mathbb{P}_{prior}(\theta)}_{R(\theta)} \right],
$$
donde en $\ell(\theta,y)$ está la Verosimilitud y $R(\theta)$ es el término de Regularización de la estadística frecuentista.

## Resumen

| Distribución del error | Función de costo | Método |
|---|---|---|
| Gaussiana i.i.d. | $\sum_i (y_i - \hat{y}_i)^2$ | Mínimos cuadrados (OLS) |
| Laplace i.i.d. | $\sum_i \lvert y_i - \hat{y}_i \rvert$ | Mínimas desviaciones absolutas |
| Gaussiana multivariada | $(y-x)^T \Sigma^{-1} (y-x)$ | Mínimos cuadrados generalizados (GLS) |
| Gaussiana&nbsp;+&nbsp;prior&nbsp;gaussiano&nbsp;sobre&nbsp;$\theta$ | OLS $+\ \lambda \lVert \theta \rVert_2^2$ | Ridge |
| Gaussiana&nbsp;+&nbsp;prior&nbsp;de&nbsp;Laplace&nbsp;sobre&nbsp;$\theta$ | OLS $+\ \lambda \lVert \theta \rVert_1$ | Lasso |


























