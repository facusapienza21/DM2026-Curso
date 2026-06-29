---
title: No11 - Programación Diferencial Pt1
---

# Programación Diferenciable: Métodos Forward

**Fecha:** 27/05/2026

:::{iframe} https://www.youtube.com/embed/SoJYG0C1kjU
:width: 100%
:::
# Programación diferenciable
Se tiene una función de costo 

$$
\mathcal{L}(\theta), \quad \theta \in \mathbb{R}^p, \quad p \gg 1
$$

donde $p$ representa los parámetros de un modelo de muchas dimensiones. 
### Optimización de Parámetros
Para resolver el problema de minimización de la función de costo:

$$
\min_{\theta \in \mathbb{R}^p} \mathcal{L}(\theta)
$$

se requiere, por ejemplo, realizar actualizaciones iterativas del parámetro $\theta$ mediante el algoritmo de descenso de gradiente. La regla de actualización en el paso $m$ se define como:

$$
\theta^{m+1} = \theta^m - \alpha^m \frac{\partial \mathcal{L}}{\partial \theta}(\theta^m)
$$

donde $\alpha^m$ representa la tasa de aprendizaje (*learning rate*) en la iteración $m$. 
Este esquema de optimización basado en gradientes es fundamental en el entrenamiento de modelos como Redes Neuronales Informadas por la Física (**PINNs**) y Ecuaciones Diferenciales Universales (**UDEs**). 

En este marco, desde una perspectiva **frecuentista**, el proceso busca obtener un estimador puntual óptimo $\theta^*$. Sin embargo, la optimización y el cálculo de estos gradientes son herramientas necesarias en ambos paradigmas:
**Frecuentista:** Para converger directamente al óptimo global o local $\theta^*$.
**Bayesiano:** Aunque el objetivo principal es obtener la distribución de probabilidad *a posteriori*, la optimización basada en gradientes sigue siendo indispensable.
Hay varios métodos que permiten calcular $\mathcal{L}(\theta)$, a nosotros nos interesan los métodos de programación diferenciable para ecuaciones diferenciales.

## Métodos de PD para ecs. diferenciables
Para UDEs y para NODEs se puede definir la función de costo asociada utilizando la norma $L_2$ al cuadrado:

$$
\mathcal{L}(\theta) = \frac{1}{N} \sum_{i=1}^{N} \| y_i - x(t_i; \theta) \|_2^2
$$

donde $x(t_i; \theta)$ representa la solución de una ecuación diferencial evaluada en el instante $t_i$. 

Dado que esta función de costo tiene dentro una ecuación diferencial, su evaluación no es analítica y requiere el uso de un *solver* numérico. A nivel de implementación, no importa estrictamente cómo se evalúa $\mathcal{L}(\theta)$ en su totalidad, ya que la computadora resuelve el problema descomponiéndolo en una secuencia de operaciones atómicas e iterativas.

Para calcular los gradientes y optimizar este sistema, los enfoques se pueden dividir según dos ejes principales, generando cuatro categorías conceptuales:

### 1. Eje del momento de discretización: Continuo vs. Discreto
* **Método Discreto (Discretizar y luego Optimizar):** El algoritmo primero discretiza la función de costo y las ecuaciones diferenciales involucradas, y luego calcula los gradientes sobre ese sistema ya discretizado.
* **Método Continuo (Optimizar y luego Discretizar):** Primero se calculan las ecuaciones diferenciales exactas que definen las sensibilidades o gradientes en el dominio continuo (diferenciación analítica) y, en un segundo paso, estas nuevas ecuaciones se resuelven numéricamente (se discretizan).

### 2. Eje de propagación de derivadas: Forward vs. Reverse
* **Modo Forward (Hacia adelante):** Propaga las derivadas desde los parámetros de entrada hacia la salida de la función. El problema de este enfoque es que **escala mal respecto a la cantidad de parámetros ($p$)**. Como nuestra función tiene $p$ entradas (los parámetros $\theta$) y 1 sola salida (el valor de $\mathcal{L}$), el modo *Forward* se ve obligado a propagar y calcular $p$ derivadas distintas en todo el grafo computacional, lo cual es ineficiente si $p \gg 1$.
* **Modo Reverse (Hacia atrás / Adjunto):** Propaga las derivadas desde la salida de la función (el error o pérdida) hacia las entradas. Al tener 1 sola salida y $p$ entradas, computacionalmente es mucho más económico propagar el gradiente "hacia atrás" en una sola pasada. Por este motivo, el modo *Reverse* (base del algoritmo de *Backpropagation*) es el estándar utilizado para entrenar redes neuronales, y representa el punto de partida óptimo para trabajar con arquitecturas de gran escala.

### Diferencias finitas

### Diferenciacion compleja
Proponemos una función de costo simplificada de un solo parámetro:

$$
\mathcal{L}(\theta), \quad \theta \in \mathbb{R}^p, \quad p = 1
$$

Si la función es **localmente analítica** (una condición que no siempre se puede garantizar teóricamente, pero que en la práctica general se cumple para la mayoría de los modelos), es posible realizar una extensión analítica de $\mathcal{L}$ hacia el plano complejo. De esta forma, el dominio de la función pasa a aceptar parámetros complejos:

$$
\mathcal{L}(\theta), \quad \theta \in \mathbb{C}^p
$$

A nivel operativo, si partimos de una variable real $x \in \mathbb{R}$ y la extendemos a una variable compleja $z = x + iy$ (con $z \in \mathbb{C}$), las transformaciones de las funciones elementales se comportan de la siguiente manera:

* **Funciones exponenciales:** $e^x \longrightarrow e^z$
* **Funciones trigonométricas:** $\cos(x), \sin(x) \longrightarrow \cos(z), \sin(z)$
* **Funciones polinómicas:** $P(x) \longrightarrow P(z)$

 **Nota computacional:** Extender las funciones al plano complejo es la base conceptual del método de diferenciación numérica por **Paso Complejo** (*Complex-Step Derivative*). Esta técnica permite calcular gradientes perturbando el sistema en el eje imaginario, lo que evita por completo los errores por cancelación catastrófica (resta de números muy parecidos) que sufren los métodos discretos tradicionales como las diferencias finitas.

Definiendo a nuestro parámetro en el dominio complejo como $z = x + iy$ (donde $z = \theta$ cuando $y=0$), podemos expresar nuestra función descompuesta en su parte real e imaginaria:

$$
\mathcal{L}(z) = u(x,y) + i v(x,y)
$$

Si $\mathcal{L}$ es localmente analítica (diferenciable en el sentido complejo) alrededor de $z$, podemos aplicar:

>### **Teorema de Cauchy-Riemann**
>**Hipótesis:** Sea una función de variable compleja $\mathcal{L}(z) = u(x,y) + i v(x,y)$ que es analítica en un entorno del punto $z$.
> 
>**Resultado:** Las derivadas parciales de primer orden de $u$ y $v$ existen, son continuas y satisfacen las ecuaciones:
> $$\frac{\partial u}{\partial x} = \frac{\partial v}{\partial y} \quad \text{y} \quad \frac{\partial u}{\partial y} = -\frac{\partial v}{\partial x}$$

De este modo, sabemos que la derivada de la función respecto al parámetro real $\theta$ equivale a la derivada parcial respecto a $x$. Y por las ecuaciones de Cauchy-Riemann, esto es igual a la derivada de la parte imaginaria respecto a $y$:

$$
\frac{d\mathcal{L}}{d\theta} = \frac{\partial u}{\partial x} = \frac{\partial v}{\partial y}
$$

Por lo tanto

$$
\frac{\partial v}{\partial y} \approx \frac{v(x, y + \epsilon) - v(x, y)}{\epsilon}
$$

Dado que partimos de un parámetro estrictamente real $\theta$, estamos evaluando en $y=0$. Además, como $\mathcal{L}(\theta)$ devuelve una función de costo real, su parte imaginaria inicial es cero ($v(x,0) = 0$). Reemplazando esto en la expresión:

$$
\frac{d\mathcal{L}}{d\theta} \approx \frac{\text{Im}[\mathcal{L}(\theta + i \epsilon)] - 0}{\epsilon}
$$

En conclusión, si $\mathcal{L}$ admite extensión compleja, podemos calcular su derivada de la siguiente manera:

$$
\frac{d\mathcal{L}}{d\theta} \approx \frac{\text{Im}[\mathcal{L}(\theta + i \epsilon)]}{\epsilon}, \quad \text{con } \epsilon \ll 1
$$


**Pros:**
* **Precisión extrema:** A diferencia de las diferencias finitas convencionales, no sufre errores de cancelación catastrófica al restar números muy cercanos, porque no hay resta en el numerador. Permite usar valores de $\epsilon$ tan chicos como la precisión de la máquina lo permita.
* **Fácil implementación:** La mayoría de los lenguajes de programación modernos tienen soporte nativo para aritmética de números complejos.

**Contras:**
* **Restricción de dominio:** Depende de que podamos evaluar y extender toda la función de costo $\mathcal{L}$ a variables complejas. No todas las funciones admiten esto (por ejemplo, funciones que usan valores absolutos o condiciones lógicas muy rígidas).
* **Sobrecarga computacional:** Matemáticamente, este enfoque termina siendo equivalente a la Diferenciación Automática (*Forward AutoDiff*), pero requiere operar algebraicamente con números complejos en cada paso, lo que consume más memoria y tiempo de cómputo.
* **Error residual:** Aunque elimina el error de cancelación, sigue siendo una aproximación numérica sujeta al error de truncamiento del orden de $\mathcal{O}(\epsilon^2)$.

Tanto el método de diferencias finitas como el de diferenciación compleja presentan un error de aproximación (de orden $\epsilon$).
Por el contrario, Forward AD nos permite obtener la derivada numérica exacta, sin error de aproximación.

### 3. Diferenciación Automática Forward (Forward AD)

#### Concepto: Grafo Computacional

Para entender AD, primero debemos modelar nuestra función como un Grafo Dirigido Acíclico (DAG).
En este esquema, definimos un conjunto de variables de entrada, variables intermedias y una variable de salida.
Las variables de entrada se denotan como $v_{-p+1}, v_{-p+2}, \dots, v_0$.
Estas variables corresponden a los parámetros del modelo, por lo que forman el vector de entrada $\theta \in \mathbb{R}^p$.
Luego, tenemos las variables intermedias, que se denotan desde $v_1$ hasta $v_{m-1}$.
Finalmente, la variable $v_m$ representa el output de nuestra función.
En este grafo computacional, una variable $v_j$ depende de $v_i$ siempre que $i < j$.
Las aristas del DAG representan las operaciones elementales que conectan a las variables de una capa con la siguiente.

**Ejemplo: $f(x) = \sin(x^2)$**

Podemos representar esta función mediante un DAG de tres capas.
En la primera capa, definimos la variable de entrada $v_0 = x$.
La arista hacia la segunda capa aplica la operación $t \to t^2$, generando la variable intermedia $v_1 = v_0^2$.
La arista hacia la tercera capa aplica la operación $t \to \sin(t)$, generando el output $v_2 = \sin(v_1)$.
De esta manera, obtenemos el resultado final $v_2 = \sin(x^2)$.

En el DAG general, nos interesa calcular la derivada del output con respecto a una de las entradas, es decir, $\frac{\partial v_m}{\partial v_{-p+1}}$.
Para ello, utilizamos la Fórmula de Bauer.
La Fórmula de Bauer establece que $\frac{\partial v_i}{\partial v_j}$ (con $i > j$) es igual a la sumatoria, sobre todos los caminos posibles $w_0 \to w_1 \to \dots \to w_k$ (donde $w_0 = v_j$ y $w_k = v_i$), del producto de las derivadas locales $\frac{\partial w_{k+1}}{\partial w_k}$.

$$\frac{\partial v_i}{\partial v_j} = \sum_{\text{caminos}} \prod_{k=1}^{K-1} \frac{\partial w_{k+1}}{\partial w_k}$$


#### Implementación de Forward AD: Números Duales

Una manera de implementar Forward AD es mediante el uso de Números Duales.
Un número dual extiende los números reales introduciendo una componente abstracta $\epsilon$.
Esta componente cumple con la propiedad de que $\epsilon^2 = 0$, con $\epsilon \neq 0$.
Un número dual se escribe de la forma $x_\epsilon = x_1 + \epsilon x_2$.
En este $x_1$ es el valor real de la variable y $x_2$ representa la variable derivada.
Ambos coeficientes, $x_1$ y $x_2$, pertenecen a los números reales.

**Propiedades de los Números Duales:**

Si tenemos dos números duales $x_\epsilon = x_1 + \epsilon x_2$ e $y_\epsilon = y_1 + \epsilon y_2$, se cumplen las siguientes propiedades:
- **Suma:** $x_\epsilon + y_\epsilon = (x_1 + y_1) + \epsilon (x_2 + y_2)$.
- **Producto:** $x_\epsilon \cdot y_\epsilon = (x_1 \cdot y_1) + \epsilon (x_1 \cdot y_2 + x_2 \cdot y_1)$.
Notemos que la componente $\epsilon$ del producto es estructuralmente idéntica a la regla de la derivada del producto.

En este contexto, $x_1$ almacena el valor de la variable original y $x_2$ almacena la derivada de esa variable con respecto a un parámetro.
Por simplicidad, asumiendo $p=1$, tenemos que $x_2 = \frac{\partial x_1}{\partial \theta}$ e $y_2 = \frac{\partial y_1}{\partial \theta}$.
Reemplazando en la regla del producto, obtenemos $x_\epsilon \cdot y_\epsilon = x_1 y_1 + \epsilon \left( x_1 \frac{\partial y_1}{\partial \theta} + y_1 \frac{\partial x_1}{\partial \theta} \right)$.
Esto equivale directamente a $x_\epsilon \cdot y_\epsilon = x_1 y_1 + \epsilon \frac{\partial (x_1 y_1)}{\partial \theta}$.

**Implementación Computacional**

Para implementar esto, extendemos el concepto de "número" en nuestro código al de "número dual".
De esta forma, cada operación matemática atómica sabe cómo multiplicar números duales y, en consecuencia, propaga la derivada automáticamente.
En el **Modo Forward** (a diferencia del modo Reverse), la evaluación avanza desde las entradas hacia la salida.
Para ello, inicializamos nuestras variables de entrada emparejando su valor con su derivada direccional.
Por ejemplo, el número dual inicializado para la entrada $v_{-p+1}$ sería:
$v_{-p+1} + \epsilon \frac{\partial v_{-p+1}}{\partial \theta_1}$.

