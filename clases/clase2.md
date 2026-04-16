---
title: Clase No2 - ODEs y NODEs
---

# ODEs y NODEs

**Fecha:** 13/04/2026

:::{iframe} https://www.youtube.com/embed/kgRSMKC8Rrg
:width: 100%
:::

% Incluir espacio

Estas notas cubren la definición de Ecuaciones Diferenciales Ordinarias (ODEs), su aplicación en modelado dinámico, la conversión de otros tipos de ecuaciones a ODEs, y cómo se utilizan para la estimación de parámetros a partir de datos observacionales, culminando con la introducción a las Neural ODEs.



## Ecuaciones diferenciales ordinarias (ODEs)

Las {term}`ODE`s describen la evolución de un sistema en función de una única variable independiente, típicamente el tiempo ($t$).

### Definición y Parámetros

* **Variable de estado:** $x(t) \in \mathbb{R}^n$ representa el estado del sistema en el instante $t$.
* **Dinámica:** Está determinada por una función $f$, conocida como {term}`campo vectorial <Campo vectorial>`, y un conjunto de {term}`parámetros <Parámetro>` del modelo $\theta \in \mathbb{R}^p$.
* **Forma general:** $$\frac{dx}{dt} = f(x, t, \theta)$$
  donde $f: \mathbb{R}^n \times \mathbb{R} \times \mathbb{R}^p \rightarrow \mathbb{R}^n$.
* **Condición inicial:** $x(t_0) = x_0$. Define el estado de partida del sistema ({term}`condición inicial <Condición inicial>`) en el tiempo inicial $t_0$.

Estas soluciones no siempre son analíticas, de hecho, casi nunca son analíticas, pero no nos importa, en este curso vamos a intentar de encontrar herramientas para resolverlas numéricamente.

### Ejemplos de aplicación y transformación

#### Reducción de Orden: El Oscilador Armónico

Una ecuación diferencial de segundo orden, como la del oscilador armónico forzado:
$$\frac{d^2x}{dt^2} + \omega^2x = f(x,t)$$
no se presenta inicialmente en la forma estándar de una ODE de primer orden ($\dot{x} = f(x,t)$). Sin embargo, es posible transformar cualquier ecuación de orden superior en un sistema de {term}`ODE`s de primer orden mediante la definición de variables de estado adicionales.

**Conversión a sistema de primer orden**
Definimos la velocidad como una nueva variable de estado $v = \frac{dx}{dt}$. Esto nos permite descomponer la ecuación original en un sistema de dos ecuaciones de primer orden:
1. $\frac{dx}{dt} = v$
2. $\frac{dv}{dt} = -\omega^2x + f(x,t)$

**Representación matricial**
Para el vector de estado $\mathbf{u}(t) = \begin{pmatrix} x \\ v \end{pmatrix}$, el sistema se expresa de forma compacta:
$$\frac{d}{dt}\begin{pmatrix}
x \\ v
\end{pmatrix} = \begin{pmatrix}
v \\ -\omega^2x+f
\end{pmatrix}=\begin{pmatrix} 0 & 1 \\ -\omega^2 & 0\end{pmatrix}\begin{pmatrix}x \\ v \end{pmatrix}+\begin{pmatrix}0 \\ f \end{pmatrix}$$

La idea central es que siempre podemos llevar una ecuación de orden más alto a un sistema de ODEs de primer orden aumentando la dimensionalidad del vector de estado.

#### Ecuaciones Diferenciales Parciales (PDEs) y el Método de Líneas

El enfoque de reducir problemas a sistemas matriciales también puede aplicarse a Ecuaciones Diferenciales Parciales (PDEs) mediante la discretización espacial, una técnica conocida como el **Método de Líneas**.

**Ejemplo: Ecuación de Difusión**
Consideremos un campo $u(x,t)$ que evoluciona según la ecuación de difusión:
$$\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2}$$
Para resolver este sistema, también debemos especificar condiciones iniciales y de borde (por ejemplo, $u(0,t)=0$ y $u(1,t)=1$).

**Discretización espacial**
Para transformar esta PDE en un sistema de {term}`ODE`s, aproximamos la segunda derivada espacial utilizando diferencias finitas centradas con un paso espacial $h$:
$$\frac{\partial^2 u}{\partial x^2} \approx \frac{u(x+h) - 2u(x) + u(x-h)}{h^2}$$
Al aplicar esto, la derivada espacial desaparece y nos queda una ecuación que depende únicamente de una derivada temporal, convirtiéndose efectivamente en una ODE.

**Representación matricial**
Definiendo un vector de estado $\mathbf{u}$ con los valores discretizados, obtenemos el siguiente sistema matricial:

$$\frac{d}{dt} \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_n \end{pmatrix} = \frac{D}{h^2} \begin{pmatrix} -2 & 1 & 0 & \dots \\ 1 & -2 & 1 & \dots \\ 0 & 1 & -2 & \dots \\ \vdots & \vdots & \vdots & \ddots \end{pmatrix} \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_n \end{pmatrix}$$

Donde la matriz de coeficientes es tridiagonal, con $-2$ en la diagonal principal y $1$ en la sub y supra diagonal.

## 1. El Modelo Dinámico: Lotka-Volterra (Depredador-Presa)

Este modelo describe la dinámica temporal del estado de un sistema compuesto por dos poblaciones: conejos ($x$) y lobos ($y$).

**Intuición de la dinámica:**
* Si los conejos están solos, se reproducen exponencialmente: $\frac{dx}{dt}=\alpha x$.
* Si los lobos están solos, mueren exponencialmente: $\frac{dy}{dt}=-\beta y$.
* Para modelar la interacción (los lobos se comen a los conejos), se agrega un término no lineal $xy$.

**Ecuaciones del modelo:**
$$\frac{dx}{dt} = \alpha x - \gamma xy$$
$$\frac{dy}{dt} = -\beta y + \eta xy$$

**Componentes del sistema:**
* **Vector de estado:** $\begin{pmatrix} x \\ y \end{pmatrix}(t; \theta)$. Depende del tiempo $t$ y está parametrizado por $\theta$. Para cada conjunto de parámetros, las curvas temporales (trayectorias) lucirán diferentes. Suelen ser oscilatorias.
* **Condición inicial:** $x(t_0)=x_0$ e $y(t_0)=y_0$.
* **Vector de parámetros:** $\theta = \begin{pmatrix} \alpha \\ \gamma \\ \beta \\ \eta \end{pmatrix} \in \mathbb{R}^4$. Como hay 4 parámetros y términos cruzados, el sistema es altamente no lineal.
    * $\alpha > 0$: tasa de crecimiento de conejos.
    * $\gamma > 0$: tasa de depredación.
    * $\beta > 0$: tasa de mortalidad de lobos.
    * $\eta > 0$: tasa de crecimiento de lobos por interacción.

## 2. Modelo Observacional y Ruido

En la realidad, si uno tuviera conocimiento absoluto de la dinámica, conocería la trayectoria perfecta. Sin embargo, nunca se observan estas trayectorias puras. Uno observa datos (no continuos en el tiempo) que se asemejan a la trayectoria, pero contaminados con ruido aleatorio.

**Ecuación del modelo observacional:**
$$x_i^{\text{obs}} = x(t_i; \theta) + \varepsilon_i$$
* $x_i^{\text{obs}}$: Observación en el tiempo $t_i$.
* $x(t_i; \theta)$: Señal o realización perfecta del modelo dinámico.
* $\varepsilon_i$: Ruido estadístico observacional.

**Características del ruido estadístico ($\varepsilon_i$):**
* **Caso estándar:** Se asume que los ruidos son independientes e idénticamente distribuidos (i.i.d.) de forma Gaussiana: $\varepsilon_i \sim N(0, \sigma^2)$, con valor medio nulo y varianza constante para cada sitio muestreado. Generalmente se supone que este ruido no depende del valor de $x$ (aunque no siempre es cierto).
* **Ruido correlacionado:** Común en series de tiempo, donde la correlación entre dos errores es distinta de cero: $\mathbb{E}[\varepsilon_i, \varepsilon_j] \neq 0$ para $i \neq j$.
    * *Nota estadística:* Si dos distribuciones son Gaussianas y tienen correlación $0$, son independientes. Si no son Gaussianas, tener correlación $0$ no implica independencia.

## 3. Inferencia Estadística y Ajuste de Trayectorias

El problema central es: dadas las observaciones $x_i^{\text{obs}}$ e $y_i^{\text{obs}}$, ¿cómo estimamos los parámetros $\theta$ que mejor describen la dinámica subyacente?

**Trajectory Matching (Fiteo de Trayectorias):**
El objetivo es convertir esto en un problema de optimización, buscando minimizar una función de costo $\mathcal{L}(\theta; \text{DATOS})$ que compare las observaciones reales con las trayectorias generadas por el modelo $x(t; \theta)$.

**Mínimos Cuadrados No Lineales:**
A diferencia del fiteo lineal de una recta, acá se compara con una función que depende de $\theta$ de manera muy no trivial.
$$\min_{\theta} \sum_{i=1}^{N} (x_i^{\text{obs}} - x(t_i; \theta))^2$$
* Minimizar el **cuadrado de los errores** asume inherentemente un ruido de naturaleza Gaussiana.
* Minimizar la **norma** asume un ruido de naturaleza Laplaciana.

---

## 4. Generalización de Interacciones (Redes Neuronales)

¿Qué ocurre si la dinámica es más compleja y desconocemos la forma exacta de la interacción entre especies? Podemos parametrizar y generalizar las funciones de interacción en lugar de usar formas fijas:

$$\frac{dx}{dt} = \alpha x - f(x,y)$$
$$\frac{dy}{dt} = -\beta y + g(x,y)$$

**Enfoque por Diccionario de Funciones:**
Podemos parametrizar $f(x,y)$ y $g(x,y)$ como combinaciones lineales de una "base" o diccionario de funciones prescritas ($f_j, g_j$), buscando cubrir lo mejor posible el espacio de funciones:
$$f(x,y) = \sum_{j=1}^k c_j \cdot f_j(x,y)$$
El problema de estimación ahora se amplía para inferir los coeficientes $c_j$ junto a $\alpha$ y $\beta$.

**Enfoque de Redes Neuronales (Neural ODEs):**
Matemáticamente, agarramos nuestras observaciones en 2D y las mapeamos a un espacio "oculto" de mayor dimensión. 
Por ejemplo: 
$$h_1 = G_1\left(M_1 \begin{pmatrix} x \\ y \end{pmatrix}\right)$$
Donde $G_1$ es una función no lineal y $M_1$ es una matriz de pesos. A medida que sumamos capas en la red, mapeamos vectores iterativamente.

Si generalizamos por completo (sin suponer un término de crecimiento o muerte específico), obtenemos una **Neural ODE (NODE)**:
$$u \in \mathbb{R}^n$$
$$\frac{du}{dt} = \text{NN}(u; \theta)$$
*Advertencia:* Al llevar el modelo a este nivel de complejidad no lineal, un problema frecuente en la optimización es que es muy fácil caer en regiones del espacio de parámetros donde matemáticamente no existen soluciones al problema.

## Ecuaciones Diferenciales Ordinarias Neuronales (NODEs)

Las **Neural Ordinary Differential Equations (NODEs)**, introducidas formalmente por {cite}`chen2018neural`, representan un cambio de paradigma al fusionar el aprendizaje profundo con los sistemas dinámicos continuos. 

### 1. Generalización de Modelos Clásicos
Una forma intuitiva de entender las NODEs es viéndolas como una generalización de modelos dinámicos preexistentes. Por ejemplo, en el **modelo de Lotka-Volterra**, en lugar de imponer una forma matemática rígida para las interacciones biológicas, podemos reemplazar o aumentar esos términos utilizando redes neuronales:

$$\frac{dX}{dt} = \alpha X - \text{NN}_1(X, Y)$$
$$\frac{dY}{dt} = -\beta Y + \text{NN}_2(X, Y)$$

En este caso, $\text{NN}_1$ y $\text{NN}_2$ son redes neuronales. Sus parámetros internos (pesos y sesgos) pasan a formar parte del vector global de parámetros a estimar, $\theta$.

### 2. Definición Formal
En su forma más general y abstracta, una NODE parametriza la derivada del estado continuo de un sistema directamente a través de una red neuronal. Se define como:

$$\frac{du}{dt} = \text{NN}(u; v)$$

Donde:
* $u$ representa el estado del sistema en un tiempo dado.
* $v$ engloba los parámetros (pesos y sesgos) de la red neuronal.

### 3. Propiedades y Consideraciones Clave
* **Aproximación universal:** Al estar basadas en redes neuronales, las NODEs heredan la capacidad de ser aproximadores universales. Esto significa que tienen la flexibilidad necesaria para aprender y representar una gama casi ilimitada de dinámicas continuas a partir de datos empíricos.
* **Estabilidad y regularización:** Al entrenar una NODE, es fundamental incorporar técnicas de regularización. Esto no solo ayuda a prevenir el sobreajuste (*overfitting*), sino que es vital para asegurar la estabilidad de las soluciones y evitar que las trayectorias se vuelvan computacionalmente intratables para los *solvers* de las ecuaciones diferenciales.

<!-- ### Ejemplo (implementación computacional)

En Myst, podemos incluir codigo!

```julia
using DifferentialEquations
using Plots
``` -->
