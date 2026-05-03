---
title: No6 - Función de costo
---

# Funciones de costo, modelos observacional, y más optimización

**Fecha:** 27/04/2026

:::{iframe} https://www.youtube.com/embed/insx5toVEus
:width: 100%
:::

En esta clase nos vamos a centrar en responder:

## ¿De donde proviene la función de costo $\mathcal{L}$?

La función de costo se deriva del modelo observacional que conecta la ecuación de estado de mi problema con los datos.

**Ejemplo 1: distribución Gaussiana**

En nuestro ejemplo característico de Lotka-Volterra (*Depredador-Presa* ver {doc}Clase N.º 2 <./clase2>) tenemos:

$$
\mathcal{L}(\theta,y) = \sum_{i=1}^{N} \left\| x(t_i; \theta) - y_i \right\|_{2}^{2}
$$
$$
y_i = x(t_i; \theta) + \varepsilon_i
$$
$$
\frac{dx}{dt} = f(x, t, \theta)
$$

Donde $x(t_i; \theta)$ es la función de estado de mi sistema, que esta descripta en este caso (*y en los que se va a enfocar este curso*) por una ecuación diferencial y un $\varepsilon_i$ que representa el ruido observacional.

Supongamos que estamos en una dimensión ($n=1$) y que el ruido $\varepsilon_i \sim N(0,\sigma^{2})$ esta caracterizado por una distribución gaussiana donde los $\varepsilon_i$ son independientes entre si y de $x$ y están idénticamente distribuidos.

La probabilidad Gaussiana de observar $y_i$ dado $x_i$ y $\sigma$ de este modelo es:

$$
P(y_i|x_i,\sigma) = \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{(y_i-x_i)^{2}}{2\sigma^{2}}}
$$

y como las variables son independientes, la probabilidad de observar todos los puntos es:

$$
P(y_1,...,y_n)|x(t_i;\theta)) = \prod_{i=1}^{n}P(y_i|x_i,\sigma)
$$

Donde $P(y_1,...,y_n)|x(t_i;\theta))$ es la **Verosimilitud** y la llamaremos $L(\theta;y)$, este va a ser nuestro modelo probabilístico para este ejemplo que nos dice dadas nuestras distribuciones de probabilidad, como los datos desde $y_1$ hasta $y_n$ se generan aleatoriamente.

:::{note}
Este modelo incorpora como estamos modelando el estado y las observaciones de nuestro sistema.
:::

:::{tip}
Para facilitar las cuentas definimos:

$$
l(\theta;y) = log (L(\theta;y)) = \sum_{i=1}^{N} log(P(y_i|x_i;\theta))
$$

Esto lo hacemos par sacarnos de encima los productos y tener todo descrito por sumas.
:::

Para poder hacer inferencia entre los parámetros del problema y nuestro modelo probabilístico vamos a utilizar:

**El Principio de Máxima Verosimilitud:**

Este principio busca estimar los parámetros del modelo que maximicen la Verosimilitud, esta ultima nos dice que tan probable es observar los parámetros $y_1$ hasta $y_N$ dada la trayectoria observada $x(t_i;\theta)$, si calculamos esta probabilidad para distintos $\theta$ la verosimilitud va a dar distintos valores, entonces lo que queremos estimar es cuales son los parámetros $\theta$ que la maximizan.

$$
\hat{\theta}_{MLE} = Arg\max_{\theta} (L(\theta;y)) = Arg\max_{\theta} (l(\theta;y))
$$

**Notar** que la ultima igualdad es valida porque el logaritmo es una función monótona creciente.

Volviendo a nuestro ejemplo:
 
$$
l(\theta;y_i) = -\sum_{i=1}^{N} (\frac{(y_i-x_i)^{2}}{2\sigma^{2}} + log(\sqrt{2\pi}\sigma))
$$

aplicando el principio de máxima Verosimilitud :

$$
\hat{\theta}_{MLE} = Arg\max_{\theta} \left[  -\sum_{i=1}^{N} (\frac{(y_i-x_i)^{2}}{2\sigma^{2}} + log(\sqrt{2\pi}\sigma))\right]
$$

Notemos que estamos maximizando sobre la variable $\theta$, entonces nos podemos "deshacer" del termino $log(\sqrt{2\pi}\sigma)$ ya que no depende de esta variable, y sacarlo no va a cambiar el resultado esperado, pero nos va proporcionar una ecuación mucho mas simple:

$$
\hat{\theta}_{MLE} = Arg\min_{\theta}\left[  \frac{1}{2\sigma^2} \sum_{i=1}^{N} (y_i-x_i)^{2}\right] = Arg\min_{\theta} \left[  \sum_{i=1}^{N} (y_i-x_i)^{2}\right]
$$

En conclusión, para este problema, el estimador de máxima verosimilitud es el que minimiza los residuos cuadráticos.

:::{important}
En general podemos mirar un problema de optimización como uno de máxima verosimilitud.

$$
Max (L(\theta;y)) = Min(- log (\mathcal{L}(\theta,y))
$$

:::

**Ejemplo 2: distribución Laplaciana**

Asumimos $\varepsilon_i$ con una distribución de Laplace

$$
\varepsilon_i = Lap(\theta,b)
$$

Cuya probabilidad es:

$$
P(\varepsilon_i|b) = \frac{1}{2b} e^{-\frac{|\varepsilon_i|}{2b}}  \ con \ b>0
$$

Como en el ejemplo anterior, vamos a maximizar la variable $\theta$ para encontrar la función de costo.

$$
\hat{\theta}_{MLE} = Arg\max_{\theta}\left[ - \frac{1}{2b} \sum_{i=1}^{N} |y_i-x_i|\right] = Arg\min_{\theta} \left[  \sum_{i=1}^{N} |y_i-x_i|\right] 
$$

Luego la función de costo es:

$$
\mathcal{L}(\theta,y) = \sum_{i=1}^{N} |y_i-x_i|
$$

:::{note}
Es una distribución que tiene colas mas pesadas (los extremos decaen mas lentamente) a diferencia de la Gaussiana, por eso se usa para hacer estadística mas robusta.
:::

Hasta ahora vinimos haciendo máxima verosimilitud solo sobre los parámetros $\theta$, en el siguiente ejemplo veremos que pasa si $\sigma$ también es un parámetro.

**Ejemplo 3: distribución gaussiana con $\sigma_i\neq cte$**

Para una distribución gaussiana cuyo $\sigma_i$ ahora no es constante, tenemos:

$$
\mathcal{L}(\theta;y) = \sum_{i=1}^{N} \omega_i (y_i-x_i)^{2}
$$

Entonces nos queda una función de costo como una función de cuadrados mínimos pesados.

**Ejemplo 4: generalización de la distribución gaussiana**

Asumimos que los $\varepsilon_i$ están distribuidos de forma gaussiana, pero esta vez estan correlacionados entre si, esto quiere decir, que no son independientes entre si.

Su matriz de covarianza $\Sigma_{ij}$ representa el valor medio $\mathbb{E} [\varepsilon_i,\varepsilon_j]$:

$$
\varepsilon_i = \begin{bmatrix}
\varepsilon_1 \\
 .\\
 .\\
 .\\
\varepsilon_N
\end{bmatrix} \sim  N(\bar{0},\Sigma_{ij})
$$

Cuya probabilidad es:

$$
P(\varepsilon_i|\Sigma) = \left(\frac{1}{(2\pi)^{\frac{N}{2}}\left|det(\Sigma)\right|^{\frac{1}{2}}}\right) e^{-\frac{1}{2}\varepsilon^{T}\Sigma^{-1} \varepsilon}
$$

:::{note}
Esta matriz general $\Sigma_{ij}$ representa geometricamente una distribución gaussiana rotada, a esto se lo conoce como una distribución normal multivariada.
:::

Ahora la función de costo $\mathcal{L(\theta)}$ esta dada de la siguiente forma:

$$
\mathcal{L}(\theta,y) = (x-y)^{T} \Sigma^{-1} (x-y) = \left\| y -x \right\|_{\Sigma}
$$

:::{note}
Esta norma esta pesada en $\Sigma$, notar que si solo me queda la diagonal de esta matriz $\Sigma$ me devuelve la norma euclidea.
:::

## ¿Podemos encontrar siempre una biyección entre $L(\theta,y)$ y $\mathcal{L}(\theta,y)$ como venimos haciendo?

La respuesta corta es NO, pero podemos hacer lo siguiente:

Dado $\mathcal{L}(\theta,y)$ queremos encontrar $L(\theta,y)$, para ello nos vamos a definirnos una función de probabilidad a la que llamaremos $L^{*}(\theta,y)$:

$$
L^{\*}(\theta,y) =\frac{e^{-\mathcal{L}(\theta,y)}}{z(\theta)} 
$$

$$
z(\theta) = \int e^{-\mathcal{L}(\theta,y)}dy
$$

Donde $z(\theta)$ representa el factor de normalización, que puede depender del parámetro $\theta$, esto no nos pasaba en los ejemplos anteriores y como consecuencia, al recuperar la función de costo como veníamos haciendo, se nos va a agregar un termino extra.

$$
\mathcal{L}^{\*}(\theta,y) = -log(L^{\*}(\theta,y)) = \mathcal{L}(\theta,y) + log(z(\theta))
$$

Definimos $R(\theta) = log(z(\theta))$ y lo llamaremos Termino de Regularización.

Luego **el caso mas general** de un problema de optimización va a tener esta forma:

$$
\mathcal{L}(\theta,y) =  \mathcal{L_{EMPIRICA}}(\theta,y) + R(\theta)
$$

donde $\mathcal{L_{EMPIRICA}}(\theta,y)$ depende tanto de los parámetros como de los datos y $R(\theta)$ solo depende de los parámetros.

**Ejemplos:** vamos a ver distintas funciones de costo que suelen aparecer cotidianamente, con $y\in \mathbb{R}^{n}$, $x\in \mathbb{R}^{n*p}$, $\theta\in \mathbb{R}^{p}$ y los $\lambda$ son hiperparametros del problema.

:::{important}
Todos los ejemplos que vamos a mostrar tienen solución analítica exacta.
:::

**1) Regresión lineal Ridge**   

$$
\min_{\theta} \underbrace{\left\| y - x \theta \right\|^{2}_{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda \left\| \theta \right\|^{2}_{2}}_{\text{Termino de Regularización}}
$$

El termino de Regularización penaliza la norma dos del vector, y se lo llama **Ridge**, esto provoca que el parámetro $\theta$ no se mueva en demasía, es decir, que tienda a converger a cero, de manera tal que cuando ingresen nuevos datos en el programa la curva ajuste mejor.

**2) Regresión Lineal Lasso**

$$
\min_{\theta} \underbrace{\left\| y - x \theta \right\|^{2}_{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda \left\|\theta \right\|_{1}}_{\text{Termino de Regularización}}
$$

El termino de Regularización penaliza la norma uno del vector, y se lo llama **Lasso**, esto provoca esparcidad en las soluciones.

**3) regresión lineal Flastic-Next**

$$
\min_{\theta} \underbrace{\left\| y - x \theta \right\|^{2}_{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda (\alpha \left\| \theta  \right\|_{1} + (1-\alpha) \left\| \theta \right\|^{2}_{2})}_{\text{Termino de Regularización}} \ con \ \alpha \in [0,1]
$$

Combinación de los ejemplos 1 y 2.

**4) Smoothin Splines**

En el espacio que vamos a estar pensando es en el espacio funcional de dimensión infinita que contiene a todas las funciones que tienen hasta la segunda derivada continua.

$$
\min_{f} \underbrace{\sum_{i=1}^{N} (y_{i} - f(x_{i}))^{2}}_{\text{Función de Costo Empírica}} + \underbrace{\lambda \int_{x_{0}}^{x_{1}}(f^{''}(x))^{2}dx}_{\text{Termino de Regularización}}
$$

El termino de Regularización penaliza la segunda derivada, lo que impone suavidad sobre las posibles soluciones.

Podemos observar que pasa cuando variamos el $\lambda$:

- $Si\ \lambda \longrightarrow \infty \ regresion \ Lineal $

- $Si\ \lambda = 0 \ Interpolacion$

## Hasta ahora resolvimos nuestros problemas pensándolos con estadística frecuentista, ¿Como la podemos conectar con la estadística Bayesiana?

En la estadística Bayesiana vamos a tener:

1) Verosimilitud, la probabilidad $\mathbb{P}(y|\theta)$ que nos dice como los datos están generados en función de los parámetros $\theta$.

2) $\mathbb{P}(\theta)$ que es una distribución de probabilidad sobre $\theta$, donde $\theta$ es **aleatorio**.

Utilizando la definición de la Probabilidad Condicional:

$$
\mathbb{P}_{post}(\theta,y) = \frac{\mathbb{P}(y,\theta) \mathbb{P}_{prior}(\theta)}{\mathbb{P}(y)}
$$

Si podemos calcular esta distribución, no solo vamos a obtener el $\theta$ que maximiza nuestro modelo, sino que también, nos va a dar una **noción de la incertidumbre** alrededor de ese $\theta$.

Buscando quien maximiza la distribución $\mathbb{P}_{post}(\theta,y)$, podemos recuperar la solución encontrada con la estadística frecuentista.

$$
\theta_{Mab} = \max_{\theta}\mathbb{P}_{post}(\theta,y) = \min_{\theta}-\underbrace{log(\mathbb{P}(y,\theta))}_{l(\theta,y)} - \underbrace{log(\mathbb{P}_{prior}(\theta))}_{R(\theta)}
$$

Donde en $l(\theta,y)$ esta la Verosimilidad y $R(\theta)$ es el termino de Regularización de la estadística frecuentista.


























