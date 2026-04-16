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

## El oscilador armonico no es una ODE

La ecuacion del oscilador armonico forzado no es una ODE
$$\frac{d^2x}{dt^2}+\omega^2 x = f(x,t)$$
ya que tiene dos derivadas temporales. Pero podemos definir la velocidad $v=\frac{dx}{dt}$ y obtenemos
$$\frac{dv}{dt}+\omega^2x=f$$
entonces $$\frac{d}{dt}\begin{pmatrix}
x \\ v
\end{pmatrix} = \begin{pmatrix}
v \\ -\omega^2x+f
\end{pmatrix}=\begin{pmatrix} 0 & 1 \\ -\omega^2 & 0\end{pmatrix}\begin{pmatrix}x \\ v \end{pmatrix}+\begin{pmatrix}0 \\ f \end{pmatrix}$$.
La idea es que siempre podemos llevar una ecuacion de orden mas alto a un sistema de ODEs.

:::{note} Ecuaciones diferenciales con derivadas de mayor orden
...
:::

## Ecuaciones diferenciales con derivadas de mayor orden

Tambien podemos aplicar esto a ecs en diferencias parciales con el  metodo de lineas, tenemos que reemplazar las deriavdas espaciales con diferencias finitas
$$ \frac{\partial u}{\partial t}= D \frac{\partial^2u}{\partial x^2} \simeq D (\frac{u(x+h)-2u(x)+u(x+h)}{h^2})$$
y ahora esto solo tiene una derivada temporal, osea es una ODE aproximadamente. Tambien tenemos condiciones iniciales y de borde que habra que especificar.
Entonces nos queda un vector con $u$ discretizado, y tenemos una ecuacion matricial 
$$\frac{d}{dt}\begin{pmatrix}
u_{0}\\ u_{1} \\ \vdots \\ u_{n} 
\end{pmatrix} =  \frac{D}{h^2} \bar{\bar{M}} \begin{pmatrix}
u_{0}\\ u_{1} \\ \vdots \\ u_{n} 
\end{pmatrix}$$
donde $\bar{\bar{M}}$ es una matriz tridiagonal con 1s en la sub y supra diagonal, y -2 en la diagonal.

:::{note} Ecuaciones diferenciales en derivadas parciales
...
:::

### Ejemplos


$\frac{du}{dt}=f(u,t,\theta)$
tenemos $x$ la poblacion de conejos e $y$ la poblacion de lobos.
Si los conejos estan solos se reporucen exponencialmente $\frac{dx}{dt}= \alpha x$. Si los lobos se quedan solos, se mueren exponencialmente $\frac{dy}{dt}=-\beta y$
Si lo dejamos aca no hay interaccion, asi que agregamos un termino $xy$ para notar que los lobos se comen a los conejos, entonces tenemos
$$\frac{dx}{dt}=\alpha x - \gamma x y$$
$$\frac{dy}{dt}=-\beta x + \eta x y$$
y una condicion inicial $x(t_0)=x_0$ e $y(t_0)=y_0$.
Este es el caballito de batalla de los ecologistas, para probar cosas. El oscilador armonico...
Tenemos 4 parametros, lo cual hace que esto sea altamente no lineal. 
insert grafico de poblaciones
- Esto es un modelo para el estado de un sistema. Si tuviera conocimiento absoluto de la dinamica de algo, que es ese algo.
- El estado del sistema vamos a asumir que es el vector $(x,y)$ que depende del tiempo, pero tambien de los parametros.
- Vamos a usar notacion, el vector $\begin{pmatrix}x \\ y\end{pmatrix}(t ;\theta)$ depende de $t$, y tambien depende de los parametros $\theta$. Lo que viene despues del ; son parametros
- Las trayectorias dependen del tiempo, obvio. Pero tambien estas curvas van a estar parametrizadas por los valores de $\theta$. Por cada conjunto de parametros, las curvas temporales van a lucir diferentes. El estado del sistema es la parte temporal digamos
- La segunda parte de esto, es un modelo observacional que lo que te dice es que vos nunca observas estas trayectorias. Uno observa cosas parecidas, datos que se asemejan a la trayectoria, pero con una componente de ruido aleatorio. Ademas los datos no son continuos en el tiempo. 
- El modelo mas simple que uno puede pensar es que en el fondo, uno lo que tiene es un $x_{i}^{obs}=x(t_{i};\theta)+\epsilon_{i}$ donde el $\epsilon_{i}$ es un ruido estadistico observacional, y el $x(t_i;\theta)$ seria la realizacion perfecta, la senial perfecta.
El caso mas simple de un modelo observacional es, uno tiene una trayectoria, evalua en ciertos tiempos, que puede estar perfectamente sampleados o no, en este caso, uno simplemente observa la serie de tiempo en determinados puntos, y le agrega el ruido observacional.

$\epsilon_1, \cdots,\epsilon_N$ los ruidos observacionales. Normalmente asumimos que estos estan distribuidos de manera normal con un valor medio nulo, y con una varianza que puede o no ser constante. 

Para cada punto del sampleo tenemos el ruido estadistico $\epsilon_i$ que tiene una distribucion normal: $\epsilon_i \sim N(0,\sigma)$ independientemente para cada sitio.

Una suposicion es que los ruidos observacionales no dependen de $x$, que tampoco tiene que ser cierto necesariamente. 

Uno aca puede hacer lo que quiera dependiendo de lo que se quiera modelar.

Tambien puede pasar que la correlacion entre dos epsilons sea distinta de cero 
$\mathbb{E}[\epsilon_i,\epsilon_j]\neq0$ para $i\neq j$. Esto significa que estan correlacionados. 
- Si dos distribuciones GAUSSIANAS tiene correlacion 0, entonces son independientes.
- Si la distribucion no es gaussiana, entonces independencia no implica correlacion 0.
Puede ser que las distribuciones no sean gaussianas...
La pregunta es, cual es la dinamica que mejor describe estos casos.

## Inferencia estadística

Ahora si tenemos observaciones, y tenemos el estado del sistema. Lo mas simple es intentar de encontrar los parametros que mejor reproducen las observaciones.
Agarramos las observaciones, y compararlo con los $x(t)$ que observaria con el modelo dinamico
$$\min_{\theta}\sum_{i=1}^N(x_{i}^{obs}-x(t_{i};\theta))^2$$
y minimizo esta cantidad, me deberia dar la mejor combinacion de parametros $\theta$.
cuando minimizo el cuadrado asumo ruido gaussiano
cuando minimizo la norma asumo ruido poissoniano (H que sea poisson)
insert grafico 
Esto no es cuadrados minimos en el sentido de fitear una recta, estamos comparando con una funcion que depende de $\theta$ de una manera potencialmente muy no trivial. Esto le llamamos CUADRADOS MINIMOS NO-LINEAL. Vamos a requerir herramientas no lineales para resolver el problema.
La idea en el fondo es convertir este problema en un problema de optimizacion, donde tenemos una funcion de costo L que depende de los parametros del $\theta$ y los DATOS. AL final del dia queremos conseguir una funcion de costos 
$$\min_{\theta}\mathcal{L}(\theta;DATOS)$$
A esto le llamamos fiteo de trayectorias (o _trajectory matching_).
Uno tiene un set de datos, y nosotros queremos fitearle a esos puntos una treyectoria. Un poco en el nombre esta la intuicion fisica que atras tenemos un problema dinamico, entonces le decimos trayectoria, pero matematicamente se reduce a cuadrados minimos no lineales.

### Ejemplo (continuado)

Dadas las OBSERVACIONES $x_i^{obs}$ y $y_i^{obs}$ queremos encontrar los 4 parametros $\vec\theta= (\alpha,\beta,\gamma,\eta)$. Este es un sistema simple de solo 4 dimensiones, casi lo pordiamos hacer por fuerza bruta. 
Pero podriamos tenes un problea mas complicado que no sabemos resolver, por ejemplo si desconocemos la interaccion entre los lobos y los conejos

$$\frac{dx}{dt}=\alpha x-f(x,y)$$
$$\frac{dy}{dt}=-\beta y+g(x,y)$$
Aca el problema se hace mas complicado, y queremos optimizar las funciones. En el caso mas practico vamos a parametrizar las funciones de alguna manera. Formalmente pensamos que estamos minimizando sobre el espacio de funciones.
PARAMETRIZAMOS: $f(x,y)=\sum_j c_j f_j(x,y)$ donde las funciones $f_j$ estan prescriptas, son nuestra base de funciones. Base esta mal dicho porque no es formalmente una base, pero es nuestro DICCIONARIO de funciones. Queremos una familia de funciones que sean lo mas grandes posibles para poder cubrir lo mayor posible el espacio de funciones (ver E3 xd).

Aca entonces usamos nuestro caballito de batalla de redes neuronales.
En el fondo las NN son 
Matematicamente agarramos nuestras observaciones 2D, y la mapeamos a un espacio de mayor dimension
$h_1=G_{1}(M_{1}\begin{pmatrix}x \\ y\end{pmatrix})$
donde G es una fucion no lineal y tenemos una matriz $M_1$ que va a ser de $m_1 \times 2$ y nos mapea el vector x,y a un vector _hidden_ de mayor dimension. Mientras mas capas tiene la NN, tenemos mas funciones $G_j$ y matrices $M_j$ que paso a paso nos van mapeando vectores en vectores con diferentes pesos y qcyo.
Tambien podriamos olvidarnos de los conejos y los lobos, y tenemos una NEURAL ODE (NODE)

$$u\in \mathbb{R}^n$$ $$\frac{du}{dt}=NN(u;\theta)$$
Uno de los problemas es que es muy facil caer en regiones del espacio donde no existen soluciones al problema.


## Ecuaciones diferenciales ordinarias neuronales (NODEs)

Las {term}`NODE`s fueron introducidas por {cite}`chen2018neural` ...

### Ejemplo (implementación computacional)

En Myst, podemos incluir codigo!

```julia
using DifferentialEquations
using Plots
```
