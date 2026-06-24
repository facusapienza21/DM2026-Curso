---
title: No13 - Programación Diferencial Pt3
---

# Programación Diferenciable: Métodos Reverse Pt1

**Fecha:** 03/06/2026

:::{iframe} https://www.youtube.com/embed/bBnnk7YehsE
:width: 100%
:::

## Programación diferenciable utilizando métodos reverse

### Repaso de la clase anteirior: Diferenciación automática en modo forward

En la [clase 12](clase12.md), se estudió el método de *Automatic Differentiation (AD)* en el caso *Forward*, implementado mediante los números duales. Repasaremos a continuación los conceptos de este caso que reaparecerán en esta clase.

Para calcular derivadas mediante AD, representamos el cálculo de una función como un grafo computacional. En este grafo, los nodos de entrada corresponden a los $p$ parámetros del problema, mientras que el nodo final representa la cantidad que deseamos diferenciar. Esta situación se muestra en
`grafo-computacional`.
```{figure} figures/clase13/grafocomp.png
:width: 500px
:align: center
:name: grafo-computacional

Grafo computacional. Una cantidad de salida $v_{10}$, generalmente la función de pérdida, se escribe en términos de los parámetros del problema ($v_0,v_{-1},v_{-2},v_{-3}$) a partir de una serie de cantidades intermedias $v_1$ a $v_9$. Las aristas dirigidas en el grafo muestran como las cantidades iniciales se utilizan para calcular las siguientes, y así sucesivamente.
```

A partir de esta representación, la derivada de la salida respecto a algún parámetro puede interpretarse en términos de los caminos del grafo. En particular, según la fórmula de Bauer, dicha derivada se obtiene como la suma de los productos de derivadas correspondientes a los diferentes caminos del gráfo que los unen (ver [clase 11](clase11.md)).
$$\frac{\partial v_m}{\partial v_o}= \sum_{v_o \underset{\Gamma}{\to} v_m} \prod_{(w_k\to w_{k+1}) \in \Gamma} \frac{\partial w_{k+1}}{\partial w_k}.$$

donde $\Gamma$ denota un camino dirigido en el grafo computacional desde $v_o$ hasta $v_m$.



### Costo computacional en la fórmula de Bauer, métodos de evaluación.
La eficiencia en la evaluación de la fórmula de Bauer depende del método utilizado.

Un método *naïf* consiste en listar todos los caminos del grafo computacional, calcular el producto de derivadas asociado a cada camino y luego sumar los resultados. Este procedimiento es altamente ineficiente, ya que muchos caminos comparten subcaminos, lo que implica recalcular múltiples veces las mismas derivadas sin reutilizar resultados.

Por ejemplo, en {numref}`grafo-computacional-derivada` se muestran todas las aristas relevantes para calcular $\partial v_{10}/\partial v_0$. No es difícil notar que existe más de un camino que conecta los nodos relevantes y que una parte importante de ellos utiliza la arista que une $v_0$ con $v_1$, por lo que la derivada parcial $\partial v_1/\partial v_0$ se recalcula repetidamente a lo largo de distintos caminos.
```{figure} figures/clase13/grafocompderivada.png
:width: 500px
:align: center
:name: grafo-computacional-derivada

Mismo grafo computacional de ejemplo que {numref}`grafo-computacional`. Se muestran resaltadas las aristas que componen los distintos caminos que unen los nodos correspondientes a $v_0$ y $v_{10}$, y por lo tanto, que son relevantes para el cálculo de la derivada $\frac{v_{10}}{v_0}$.
```


### El método Forward
Otro método, empleado en la clase pasada, consiste en recorrer el grafo en un orden topológico, de manera que cada nodo se evalúa únicamente después de haber evaluado todos sus predecesores en el grafo. Esto permite reutilizar cálculos intermedios y evitar la enumeración explícita de todos los caminos, reduciendo significativamente el costo computacional.

Al propagar de manera ordenada los términos de la fórmula de Bauer, se van obteniendo derivadas parciales de nodos cada vez más “alejados” de las entradas. En este contexto, la estrategia *forward* para calcular la derivada $\partial v_j/\partial v_i$ consiste en expresar dicha cantidad en función de los predecesores inmediatos $\omega$ de $v_j$:
$$\frac{\partial v_j}{\partial v_i} = \sum_{\omega \to v_j} \frac{\partial v_j}{\partial \omega}\frac{\partial \omega}{\partial v_i}.$$
Al ser $\omega$ inmediatamente próximo a $v_j$, la derivada $\frac{\partial v_j}{\partial \omega}$ será sencilla de calcular, restando la dificultad en calcular $\frac{\partial \omega}{\partial v_i}$, lo que se hará de forma recursiva reaplicando la fórmula. 

En nuestro ejemplo (y como se muestra en {numref}`grafo-computacional-forward`), comenzaremos calculando $\frac{\partial v_0}{\partial v_0}$, luego $\frac{\partial v_1}{\partial v_0}$ y $\frac{\partial v_2}{\partial v_0}$, luego $\frac{\partial v_4}{\partial v_0}$, $\frac{\partial v_5}{\partial v_0}$ y $\frac{\partial v_6}{\partial v_0}$, luego $\frac{\partial v_8}{\partial v_0}$ y $\frac{\partial v_9}{\partial v_0}$ y por último el resultado buscado $\frac{\partial v_{10}}{\partial v_0}$.

Dado que las derivadas se calculan propagando información desde las entradas hacia la salida siguiendo el orden topológico del grafo, este procedimiento se denomina método *Forward*.
```{figure} figures/clase13/grafocompforward.png
:width: 500px
:align: center
:name: grafo-computacional-forward

Mismo caso que {numref}`grafo-computacional-derivada`. Se muestra ahora también un orden posible en que un método forward logrará calcular las distintas derivadas parciales $\frac{v_{i}}{v_0}$, con $i\in\{1,...,10\}$.
```

## El método Reverse

Análogamente al método forward, que consistía en explorar los nodos del grafo de forma ordenada partiendo de los parámetros y llegando al valor de salida, el método *reverse* consistirá en explorar los nodos del grafo en un orden similar pero partiendo del valor de salida y dirigiéndose hacia los parámetros, calculando en el proceso derivadas de la salida respectoa valores de indices menores. Así, el método reverse para calcular una derivada $\frac{\partial v_j}{\partial v_i}$ con $j>i$ también hará uso de cantidades intermedias $\omega$ que sucedan al nodo $v_i$ y se valdrá de la relación dinámica:
$$\frac{\partial v_j}{\partial v_i} = \sum_{v_i\to \omega} \frac{\partial v_j}{\partial \omega}\frac{\partial \omega}{\partial v_i}.$$
De las dos derivadas a la derecha, esta vez será $\frac{\partial \omega}{\partial v_i}$ la que es fácil de calcular mientras que para $\frac{\partial v_j}{\partial \omega}$ se tendrá que aplicar la relación recursivamente. A este método también se le da el nombre *back-propagation*.

Reusando el ejemplo anterior (como se muestra en {numref}`grafo-computacional-reverse`), se calcularán primero $\frac{\partial v_{10}}{\partial v_8}$ y $\frac{\partial v_{10}}{\partial v_9}$, luego $\frac{\partial v_{10}}{\partial v_4}$, $\frac{\partial v_{10}}{\partial v_5}$ y $\frac{\partial v_{10}}{\partial v_6}$, luego $\frac{\partial v_{10}}{\partial v_1}$ y $\frac{\partial v_{10}}{\partial v_2}$ y por último $\frac{\partial v_{10}}{\partial v_0}$.
```{figure} figures/clase13/grafocompreverse.png
:width: 500px
:align: center
:name: grafo-computacional-reverse

Mismo caso que {numref}`grafo-computacional-derivada`. Se muestra ahora también un orden posible en que un método reverse logrará calcular las distintas derivadas parciales $\frac{v_{10}}{v_i}$, con $i\in\{0,...,9\}$.
```

### Uso de memoria y _checkpointing_


Una desventaja del modo reverse respecto al modo forwards, es que requiere precomputar los valores asignados a cada nodo para poder evaluar las derivadas. 

Por ejemplo, para calcular un valor para $\frac{\partial v_{j+1}}{\partial v_{j}}$, uno necesita del valor de $v_j$ o $v_{j+1}$ en que evaluar la derivada. 

En el modo forward, estos valores se generan y utilizan de manera secuencial a medida que avanza el cálculo. En cambio, en el modo reverse, los valores intermedios obtenidos durante el pase forward deben permanecer disponibles durante el recorrido inverso del grafo, como se esquematiza en  {numref}`grafo_doscasos`. Esto implica un costo adicional en memoria para guardar todas las variables intermedias del programa y un costo en tiempo por el acceso de la memoria.
```{figure} figures/clase13/doscasos.png
:width: 500px
:align: center
:name: grafo_doscasos

Ventaja en memoria de los métodos modo forward con respecto a los modo reverse. Se esquematiza como en el caso forward, a medida que el programa calcula las diferentes cantidades intermedias del sistema (por ejemplo, resolviendo una ODE), estas ya se pueden utilizar para calcular los gradientes. En cambio, en el caso backwards, uno debe esperar a calcular la variable final antes de ir sucesivamente utilizando variables intermedias anteriores para propagar el gradiente para atrás. 
```
Para aliviar esta demanda de memoria, se suele usar la estrategia de *Checkpointing* ilustrada en {numref}`grafo-checkpointing`. Esta estrategia consiste en guardar en memoria el estado del sistema únicamente en ciertos puntos espaciados del programa. Si luego para el cálculo de derivadas mediante algún metodo backwards se requieren valores intermedios a los puntos disponibles, estos deberán poder obtenerse corriendo el programa de vuelta a partir del punto guardado más cercano. Con esta estrategia se balancea el úso de memoria y cómputo al decidir sobre que puntos guardar el estado del sistema.

```{figure} figures/clase13/checkpointing.png
:width: 300px
:align: center
:name: grafo-checkpointing

Mismo caso que {numref}`grafo_doscasos`. Se muestra la estrategia de checkpointing, en que el estado del sistema (en este caso valor de la ODE) se guarda únicamente en los puntos marcados por una cruz azul. Al calcular gradientes mediante backpropagation, el valor de la ODE en el punto a considerar debe calcularse a partir del valor en memoria más cercano.
```

Esta estrategia permite equilibrar el uso de memoria y el costo computacional, intercambiando almacenamiento por costo de calculo computacional: cuanto menor es la cantidad de estados guardados, mayor es el trabajo necesario para reconstruir los valores intermedios durante el barrido inverso.


### Cuando conviene usar reverse? VJP vs JVP


Vale la pena pensar en que casos las ventajas de los métodos *reverse*  los vuelven convenientes aún considerando las penalidades en memoria que estos implican. Para ver un ejemplo de estas ventajas, podemos considerar una versión vectorizada de un grafo computacional (ilustrado en {numref}`grafo-computacional-vectorizado`):
```{figure} figures/clase13/grafovectorizado.png
:width: 500px
:align: center
:name: grafo-computacional-vectorizado

Grafo computacional en una versión vectorizada, en que cada cantidad $h_j$ representa un vector de $d_j$ cantidades intermedias que se calculan de las cantidades del nodo anterior mediante una función $g_j$.
```


Aquí, los $p$ parámetros $\theta$ son la entrada del grafo computacional, mientras que la salida es la función de costo $\mathcal L$, un escalar real. En las capas intermedias, las variables intermedias $h_j\in \mathbb R^{d_j}$ se calculan como 

$$h_j=g_j(h_{j-1})$$



donde $g_j:\mathbb R^{d_{j-1}}\to \mathbb R^{d_j}$ es la función que relaciona las variables en cada nodo con el anterior. Así, la función de costo será 

$$ \mathcal L: \mathbb{R}^{p} \rightarrow \mathbb{R}
$$

$$\mathcal L(\theta)=g_m\circ g_{m-1}\circ...\circ g_1(\theta).$$

Empleando la regla de la cadena, el Jacobiano de la función de costo será el producto de los Jacobianos de cada $g_m$:
$$D\mathcal L=Dg_m\;Dg_{m-1}\;...\;Dg_1 (\theta)$$

con
$$D\mathcal{L} \in \mathbb{R}^{1 \times p}$$
$$D\mathcal{g_m} \in \mathbb{R}^{1 \times d_{m-1}}$$
$$D\mathcal{g_{m-1}} \in \mathbb{R}^{d_{m-1}\times d_{m-2}}$$ $$D\mathcal{g_1} \in \mathbb{R}^{d_1\times p}$$


Podemos observar que en este caso el cálculo del lado derecho de la ecuación mediante un método *backward* corresponde a la propagación de derivadas en sentido inverso sobre el grafo computacional, lo cual es equivalente a la evaluación de productos vector–Jacobian (VJP). Este enfoque resulta especialmente eficiente cuando la dimensión de la salida es pequeña, en particular cuando $d_m = 1$.

En cambio, el método *forward* corresponde a la evaluación de productos Jacobian–vector (JVP), en los cuales se propagan direcciones desde las entradas hacia la salida. Este enfoque es más eficiente cuando la dimensión de entrada es pequeña, es decir, cuando $d_0 = p$ no es grande.

Para calcular $AB$ con $A\in \mathbb R^{n\times r}$ y $B\in \mathbb R^{r\times m}$ se requieren $n\times m$ productos internos, cada uno de $r$ términos, el costo de calcular $AB$ será de orden $\mathcal O(mnr)$. En general, en un caso en que tenemos una función de costo $\mathcal L:\mathbb R^{p}\to \mathbb R^q$ con $m$ pasos intermedios, el número de operaciones a realizar en modo forward va como $\mathcal O(pm+q)$, mientras que en el modo backward irá como $\mathcal O(mq+p)$.

Esto, junto al mayor costo en memoria de aplicar metodos backwards, sugieren que esto es únicamente conveniente si $p>q$.


| Caso          | Forward (JVP) | Reverse (VJP) |
|---------------|---------------| ------------- |
| $p \ll q$     | ✓             | X             |
| $p \approx q$ | ✓             | X             |
| $p \gg q$     | X             | ✓             |


Para dar números más concretos. Dada una ODE de $n$ variables, con $q=1$, para $n+p\gtrsim 50$, convene utilizar métodos reverse, mientras que en el caso contrario conviene utilizar métodos forward.




## Método del adjunto

Este es un metodo estándar de programación diferencial utilizando ecuaciones diferenciales.
Supongamos que tenemos que resolver una ecuación diferencial $\frac{du}{dt}=f(u,\theta,t)$. 
Podemos discretizar la solución en el tiempo para buscar únicamente un conjunto $\{u_j\}$, tales que $u_j=u(t_j)$. 
Podemos concatenar los valores a cada tiempo en un *super-vector*, $u=(u_1,...,u_m)\in \mathbb R^{nm}$, tal que la expresión para la ecuación diferencial discreta es $G(u,\theta)=0$.

Podemos pensar, por ejemplo, en el caso de una ecuación diferencial lineal 
$$\frac{du}{dt}=A(\theta)u+b(\theta).$$Una discretización posible, si utilizamos el método de euler explicito, es $$u_{j+1}=(\mathbb I +\Delta t A(\theta)) u_j+\Delta t b(\theta).$$
Esta ecuación diferencial a su vez se puede reescribir como $$g_j(u,\theta)=u_{j+1}-(\mathbb I +\Delta t A(\theta)) u_j-\Delta t b(\theta)=0.$$En forma vectorial, esta resultará
$$G(u,\theta)=(g_1,...,g_m)=0.$$
Podemos considerar que además tendremos una función de costo para contrastar con datos observacionales, $L(u,\theta)$. Por ejemplo, podría ser de la forma $L(u,\theta)=\sum^m_{i=1}\omega_i||u_i^{OBS}-u_i||^2$, siendo $\omega_i$ el peso asignado a cada observación.

El objetivo del método del adjunto será calcular el gradiente de la función de costo respecto a los parámetros, $\frac{dL}{d\theta}(u,\theta)$. 
Utilizando la regla de la cadena podemos relacionar este gradiente con la sensibilidad:
$$\frac{dL}{d\theta}(u,\theta)=\frac{\partial L}{\partial u}\frac{\partial u}{\partial \theta}+\frac{\partial L}{\partial \theta}.$$
Para obtener la sensibilidad, que es el término difícil de esta ecuación resultante, se puede, por ejemplo, derivar la ecuación diferencial discreta:
$$0=\frac{dG}{d\theta}=\frac{\partial G}{\partial u}\frac{\partial u}{\partial \theta}+\frac{\partial G}{\partial \theta}\;\Rightarrow\;\frac{\partial u}{\partial \theta}= -\left(\frac{\partial G}{\partial \theta}\right)^{-1}\frac{\partial G}{\partial \theta}.$$
Esta expresión para la sensibilidad, al remplazarse en la expresión del gradiente de la función de costo otorga
$$\frac{dL}{d\theta}(u,\theta)=-\frac{\partial L}{\partial u}\left(\frac{\partial G}{\partial \theta}\right)^{-1}\frac{\partial G}{\partial \theta}+\frac{\partial L}{\partial \theta}.$$
El método del adjunto define una cantidad $\lambda^t=\frac{\partial L}{\partial u}\left(\frac{\partial G}{\partial \theta}\right)^{-1}$, llamada la variable adjunta, tal que
$$\frac{dL}{d\theta}(u,\theta)=-\lambda^t\frac{\partial G}{\partial \theta}+\frac{\partial L}{\partial \theta} \; \text{ con }\; \left(\frac{\partial G}{\partial \theta}\right)^T\lambda=\left(\frac{\partial L}{\partial u}\right)^T.$$

La [clase siguiente](clase14.md) continuaremos viendo el método del adjunto, continuando con el caso discreto y luego incluyendo también el caso continuo.