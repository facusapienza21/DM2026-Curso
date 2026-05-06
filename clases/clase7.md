---
title: No7 - Bayes
---

# Inferencia Bayesiana, Estimación del Ruido


**Fecha:** 04/05/2026

:::{iframe} https://www.youtube.com/embed/QDeZqV5sWTQ
:width: 100%
:::

# Clase anterior

Problemas de optimización sobre parámetros, función de costos y término de regularización. 

Caso general (frecuentista):

$$\min_{\theta} \mathcal{L}(\theta, y)+\mathcal{R}(\theta)$$

Aquí, $\mathcal{L}(\theta, y)$ es la función de costo, que tiene un término empírico y otro con parámetros, y $\mathcal{R}(\theta)$, el término de regularización. Regularizar es, en términos generales, agregar un bias inductivo, un sesgo para aportar información de manera intencional y condicionar el resultado a conocimientos previos. Por ejemplo, considerando los métodos vistos: obtener un vector esparso o de norma chica. El early stopping también es una forma de regularización. Todas ellas impiden el overfitting y permiten generalizar a datos nuevos.  

:::{figure} ./figures/clase_7_01.JPG 
:width: 100% 
:align: center 

Ubicación de $\theta_{MLE}$ respecto de su distribución 
:::

# Observación: cuantificación de incertidumbre

## Caso Bayesiano: es información que nos aporta el posterior

## Caso frecuentista: Bootstrap
**Bootstrap paramétrico:** 
Sabemos que $Y \curvearrowright \theta^*$ a través de la función likelihood. Luego, fijando la curva con $\theta^*$ y agregando ruido, obtenemos $\hat{Y}_1 \curvearrowright \theta^*_1, \hat{Y}_2 \curvearrowright \theta^*_2,...$ Estimadores de $Y$ con sus correspondientes $\theta^*$

:::{figure} ./figures/clase_7_02.JPG 
:width: 50% 
:align: center 

Elección de $\hat{\theta}$ entre los valores con ruido 
:::

**Bootstrap no paramétrico:**

$Y=\{Y_1, ..., Y_N\} \curvearrowright \theta^*$
Sampleando muestras con repetición, creamos "copias" de $Y$:

$\hat{Y}_1=\{Y^1_1, ..., Y^1_N\} \curvearrowright \theta^*_1$

.

.

.

$\hat{Y}_k=\{Y^k_1, ..., Y^k_N\} \curvearrowright \theta^*_k$

Si bootstrap se comporta bien, en ciertos casos particulares se obtiene el posterior

# Estimar el posterior

## Algortimos

**MCMC (Markov Chain Monte Carlo)**

Monte Carlo es un término poco preciso que se refiere a estimar o calcular algo por sampleo. Mientras tanto, Markov Chain representa cómo se da cada paso del algoritmo.

Luego, este algoritmo avanza distintos valores de $\theta$ explorando la densidad del mismo

$\theta_0 \curvearrowright \theta_1 \curvearrowright ... \curvearrowright \theta_k \curvearrowright \theta_{k+1}$

:::{figure} ./figures/clase_7_03.JPG 
:width: 50% 
:align: center 

Así exploran el espacio los algoritmos MCMC 
:::

Esta cadena se corta en un k determinado tal que el algoritmo está encaminado y se encuentra en el soporte de $\theta$

$\theta = \{\theta_k, \theta_{k+1},...,\theta_{k+m}\} \; k, m\in \mathbb{N}$ 

$\Theta \thicksim \mathbb{P}(\theta|Y)$

Con infinitos puntos, es posible describir la densidad.

# Aproximación de Laplace
Es aproximar el posterior con una distribución conocida. Hacen falta dos cosas:

* $\theta_k \curvearrowright$ proponer $\theta_{k+1}^*$

Por ejemplo, $\theta^*_{k+1} \thicksim N(\theta_k, \sigma_k^2)$, con ruido gaussiano que lo pereturba

* Aceptar o rechazar cada $\theta_k$
Para esto, definimos $$\alpha_k=\frac{\mathbb{P}(\theta^*_{k+1}|Y)}{\mathbb{P}(\theta^*_k|Y)} = \frac{\mathbb{P}(Y|\theta^*_{k+1}).\mathbb{P}(\theta^*_{k+1})}{\mathbb{P}(Y|\theta^*_{k}).\mathbb{P}(\theta^*_{k})}.\frac{\mathbb{P}(Y)}{\mathbb{P}(Y)}$$

Observemos que se cancela $\mathbb{P}(Y)$ que es lo difícil de calcular y quedan todas cosas conocidas. 

Si $\theta^*_{k+1} > \theta^*_{k}$, entonces $\alpha_k>1$

**Método de Metropolis-Hastings (M-H)**

Sampleamos $U_k \sim Unif([0,1])$ y aceptamos si $U_k \leq \alpha_k$

Se corre al infinito y sacamos el principio. Los puntos garantizan que la distribución marginal de los mismos se parece al posterior. 

Este método se parece mucho a la optimización por gradiente. Esto nos dice que la estadística Bayesiana y Frecuentista se parecen mucho. Sin embargo, la Bayesiana es más difícil de calcular porque el algoritmo de optimización busca explorar el espacio para ver la distribución de $\theta$.

**Observación** $\theta \in \mathbb{R}^p$, MCMC M-H funciona cuando $p\sim 1$. Si $p>>1$, se usa Hamiltonian-MCMC que usa el gradiente.

:::{figure} ./figures/clase_7_04.JPG 
:width: 50% 
:align: center

 Dispersión de los datos
:::

:::{figure} ./figures/clase_7_05.JPG 
:width: 50% 
:align: center 

Dispersión de los datos para el caso Lotka-Volterra 
:::