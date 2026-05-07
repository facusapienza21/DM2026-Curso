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

# Estadistica Bayesiana
Para el caso de Estadistica Bayesiana hay dos funciones importantes que debemos tener en cuenta:
* **Verosimilitud:** $\mathbb{P}(Y \mid \theta)$  (es lo que en Estadistica Frecuentista conocíamos como likelihood, es decir la función $\mathcal{L}(\theta \mid x)$ ) 
* **Prior:**  $\mathbb{P}(\theta)$, es la función del parámetro (es lo que en Estadistica Frecuentista conocíamos como el término de Regularización $R(\theta)$)

A continuación recordemos el **Teorema de Bayes**, el cual nos dice

$$P(\theta \mid Y) = \frac{P(Y \mid \theta)P(\theta)}{P(Y)}$$

Como la distribución de los datos no depende de $\theta$, luego podemos escribir a $\frac{1}{\mathbb{P}(Y)}$ como una constante $\alpha$ cuyo valor desconocemos. Por lo tanto, obtenemos

$$P(\theta \mid Y) = \alpha \cdot P(Y \mid \theta)P(\theta)$$

En este caso, queremos calcular $\mathbb{P}(\theta \mid Y)$, es decir la distribución del Posteriori. No estamos interesados en la estimación puntual de $\theta$ sino que nos interesa **Estimar una Distribución**

:::{figure} ./figures/clase_7_01.JPG 
:width: 100% 
:align: center 

Ubicación de $\theta_{MLE}$ respecto de su distribución 
:::

Queremos hallar el **Maximum a Posteriori**, es decir

$$\theta_{MAP} = \max_{\theta} P(\theta \mid Y) = \max_{\theta} P(Y \mid \theta) P(\theta)$$
Podemos tomar logaritmo pues es una función creciente y además, como sabemos que el problema de maximización puede reformularse como un problema de minimización considerando la función objetivo multiplicada por $- 1$, obtenemos finalmente
 $$\theta_{MAP} =  \min_{\theta} [ - \log(P(Y \mid \theta)) - \log(P(\theta)) ]$$

Observemos que podemos asociar:  $- \log(P(Y \mid \theta))$ con la **función de costo empírica** $\mathcal{L}_{EMP}(\theta \mid Y)$ y a su vez podemos asociar: $- \log(P(\theta))$ con el **término de Regularización** $R(\theta)$.


# Observación: Tanto en el caso Frecuentista como en el Bayesiano
Cuando el numero de observaciones $n \to \infty$ tenemos que 

$$\theta_{MAP} \to \theta_{0}$$ 
$$\theta_{MLE} \to \theta_{0}$$

Siendo $\theta_{0}$ el verdadero valor de $\theta$



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

# Estimación del Ruido
Tenemos el modelo 
$$Y_i = X(t_i; \theta) + \epsilon_i, \quad \epsilon_i \sim \mathcal{N}(0, \sigma^2)$$

Bajo este supuesto, la estimación por máxima verosimilitud (MLE) de los parámetros del modelo equivale a resolver:

$$\min_{\theta} \sum_{i=1}^{n} \left(Y_i - X(t_i;\theta)\right)^2$$

En general, $\sigma^2$ también es un parámetro desconocido del modelo y debe ser estimado. Su estimador de Máxima Verosimilitud es:

$$\hat{\sigma}^2_{MLE} = \frac{1}{n} \sum_{i=1}^{n} (Y_i - X(t_i;\theta))^2$$

:::{figure} ./figures/clase_7_04.JPG 
:width: 50% 
:align: center

 Dispersión de los datos
:::

## Observación: 
En **Estadística Bayesiana**, es necesario especificar una distribución a priori para todos los parámetros desconocidos. En este caso, debemos definir priors tanto para $\theta$ como para $\sigma^2$: $P(\theta)$ y $P(\sigma^2)$. En un escenario general, tenemos: $$Y \in \mathbb{R}^{n}, \quad Y_{ij} = X_{i}(t_{j}) + \epsilon_{ij}, \quad \epsilon_{ij} \sim \mathcal{N}(0, \sigma^2)$$ 

#### El caso de Lotka-Volterra 
:::{figure} ./figures/clase_7_05.JPG 
:width: 50% 
:align: center 

Dispersión de los datos para el caso Lotka-Volterra 
:::

En el modelo de **Lotka–Volterra**, los distintos componentes (población de presas vs. población de depredadores) pueden tener niveles de ruido diferentes. Esto significa que los $\sigma_i$ no son necesariamente iguales. Si los $\sigma_i$ son distintos, la log-verosimilitud nos lleva al siguiente problema de optimización para encontrar los parámetros:
 $$\min_{\theta, \sigma} \sum_{j} \sum_{i} \left[ \frac{1}{2\sigma_i^2} (Y_{ij} - X_i(t_j;\theta))^2 + \log(2 \pi \sigma_i^2) \right]$$


Donde $i$ es la componente y $j$ el paso temporal. Definimos los pesos: $$w_i = \frac{1}{2\sigma_i^2}$$ 

Esto significa que, cuanto menos ruidosa es la señal, más peso le damos en la función de costo. 

Sin embargo, surge un **problema**: los pesos $w_i$ dependen de los mismos parámetros $\sigma_i$ que queremos estimar. Por lo tanto, la estimación de los pesos y la de $\theta$ quedan **acopladas**.Esto nos induce a una **Estrategia de Optimización Alternada**, que se traduce en el siguiente algoritmo:


 **Algoritmo de Optimización Alternada** 
 1. **Paso $\theta$:** Dado $\hat{\sigma}_i$ actual, estimamos $\theta$ minimizando la suma pesada. 
 2. **Paso $\sigma$:** Dado $\hat{\theta}$ actual, estimamos los nuevos $\sigma_i$. 
 3. **Repetir** ambos pasos hasta alcanzar la convergencia. 




