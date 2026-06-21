---
title: No14 - Programación Diferencial Pt4
---

# Programación Diferenciable: Métodos Reverse Pt2

**Fecha:** 08/06/2026

:::{iframe} https://www.youtube.com/embed/OYmbtWjY6HE
:width: 100%
:::
**Autores:** Felipe Cignoli (`@fcignoli`), Martin Sinnona, Noe Hsueh



## Método del Adjunto Continuo

Consideremos una ODE de primer orden dada por

$$
\frac{du}{dt} = f(u,\theta,t)
$$ (eq-ode)

sujeta a la condición inicial $u(t_0)=u_0$, donde $u \in \mathbb{R}^n$ es el vector solución desconocido de la ODE, $f:\mathbb{R}^n \times \mathbb{R}^p \times \mathbb{R} \to \mathbb{R}^n$ es una función que depende del estado $u$, $\theta \in \mathbb{R}^p$ es un vector de parámetros, y $t \in [t_0,t_1]$ se refiere al tiempo. Aquí, $n$ denota el tamaño de la ODE y $p$ el número de parámetros. Resolver la ODE implica obtener $u(t)$, que depende de $\theta$. En general no va a ser posible obtener la solución explicita de $u$. 

Recordemos que queremos obtener $\theta$ (donde $\theta$ puede ser parametros de una red o en problemas inversos coeficientes de ecuaciones diferenciales). De esta forma nos interesa generalmente $\frac{dL}{d\theta}$. 

Para ello, veamos que es $L$. Se puede escribir el término de la loss más general de la forma integral:

$$
L(u(\cdot, \theta);\theta)=\int_{t_0}^{t_1} h(u(\tau;\theta),\theta)\,d\tau
$$ (eq-loss)

donde $h$ es una funcion de costo puntual. 

- ¿Por qué podemos escribirlo como una integral?
    
    Veamos el ejemplo del caso discreto, donde la Loss es  $w_i \|u(t, \theta)-u^{\text{obs}}(t)\|_2^2.$ Podemos escribir la loss como 
    
    $$
    h(u, \theta, t) = \sum_{i} \|u(t; \theta) - u_i^{\text{obs}}\|^2 \delta(t - t_i),
    $$
    
    donde $\delta$ es la funcion delta de dirac. Si tomamos la integral de $h$ de $t_0$ a $t_1$ recuperamos la suma. 
    

Ahora, derivemos {eq}`eq-loss` con respecto de $\theta$ usando regla de la cadena:

$$
\frac{dL}{d\theta} = \int_{t_0}^{t_1} \left( \frac{\partial h}{\partial \theta} + \frac{\partial h}{\partial u} \underbrace{{\frac{\partial u}{\partial\theta}}}_{s(t)} \right) dt.
$$ (eq-dloss)

En {eq}`eq-dloss` notamos que tenemos un VJP $\frac{\partial h}{\partial u} {s(t)}$, con sensibilidad $s(t)\in \mathbb{R}^{n\times p}$ (matriz por vector, notar que $\partial h/\partial u$ es $1\times n$ (h es una función escalar, u es un vector, luego el gradiente es un vector de n num) cuyo resultado es un vector de $1\times p$. Como en el método discreto, la idea *del método adjunto consiste en* aprovechar esto para introducir una nueva variable (*adjunto)* $\lambda$ que nos permita evitar calcular el jacobiano $s(t).$

:::::{margin}
:::{dropdown} ¿Qué era sensibilidad?

$s(t)$ define qué tanto cambia mi solución $u(t)\in\mathbb{R}^n$ con respecto de $\theta$: $\frac{\partial u }{\partial \theta}$. Notemos que tiene su ecuacion diferencial asociada. Diferenciemos {eq}`eq-ode` con respecto de $\theta$

$$
\begin{align*}
\frac{d}{d\theta} \left[ \frac{d}{dt}u(t; \theta) \right] = \frac{d}{d\theta} \left[ f(u(t; \theta), \theta, t) \right] \tag{i}
\end{align*}
$$

Ahora, intercambiamos orden de derivadas parciales (asumimos que vale, solucion suave, …)

$$
\begin{align*}
\frac{d}{dt} \left[ \frac{du}{d\theta} \right] = \frac{d}{d\theta} f(u(t; \theta), \theta, t) \tag{ii}
\end{align*}
$$

En el lado deracho desarrollamos la derivada:

$$
\begin{align*}
\frac{d}{dt} \left[ \frac{du}{d\theta} \right] = \frac{\partial f}{\partial u} \frac{\partial u}{\partial \theta} + \frac{\partial f}{\partial \theta} \tag{iii}
\end{align*}
$$

Definamos ahora la matriz de sensibilidad como: $s(t) = \frac{\partial u}{\partial \theta}$ y reemplazamos, obteniendo asi la ecuacion de sensibilidad.

$$
\begin{align*}
\frac{ds}{dt} = \frac{\partial f}{\partial u} s + \frac{\partial f}{\partial \theta} \tag{iv}
\end{align*}
$$

Su condicion inicial esta dada por $u_0$:

$$
\begin{align*}
s(t_0) = \frac{du(t_0)}{d\theta} \tag{v}
\end{align*}
$$

y si $u_0$ no depende de $\theta$, $s_0=0$

De esta forma, la ecuacion de sensibilidad me dice como un cambio de los parametros afecta a mi solucion del sistema en el tiempo.
:::
:::::

Notemos que la sensibilidad tiene su ecuación diferencial asociada: 

$$
\frac{ds}{dt} = \frac{\partial f}{\partial u} s + \frac{\partial f}{\partial \theta} 
$$

Reordenamos los términos

$$
\left[ \frac{ds}{dt} = \frac{\partial f}{\partial u} s + \frac{\partial f}{\partial \theta} \right]
\Rightarrow \frac{ds}{dt} - \frac{\partial f}{\partial u} s - \frac{\partial f}{\partial \theta} = 0
$$

Esta expresion es cero para cualquier instante de tiempo, podemos multiplicar por $\lambda(t)^\top$ e integrarlo y va a seguir siendo cero. 

$$
\begin{align*}
\Rightarrow \int_{t_0}^{t_1} \lambda(\tau)^\top \left[ \frac{ds}{dt} - \frac{\partial f}{\partial u} s - \frac{\partial f}{\partial \theta} \right] d\tau &= 0\quad \forall \lambda(t):[t_0, t_1] \mapsto \mathbb{R}^n
\end{align*}
$$ (eq-integral-constraint)

Recordamos que el objetivo de ahora es sacarnos la sensibilidad $s(t)$. Primero vamos a usar integracion por partes sobre $\lambda^\top \frac{ds}{dt}$ para despejar esta derivada con respecto del tiempo, luego lo reemplazamos en {eq}`eq-integral-constraint`. 

$$
\begin{align*}
0 &= \int_{t_0}^{t_1}\left[
\lambda^\top \frac{ds}{dt}
-
\lambda^\top \frac{\partial f}{\partial u}s
-
\lambda^\top \frac{\partial f}{\partial \theta}
\right]d\tau \\
&=
-\int
{\color{red}\frac{\partial\lambda^\top}{\partial t}s}\,d\tau
+
\left. \lambda s \right|_{t_0}^{t_1}
+
\int_{t_0}^{t_1}
{\color{red}-\lambda^\top \frac{\partial f}{\partial u}s}
-
\lambda^\top \frac{\partial f}{\partial \theta}
d\tau \\
&= \lambda(t_1) s(t_1) +
\int_{t_0}^{t_1}
{\color{red}\left( -\frac{\partial \lambda^\top}{\partial t} -\lambda ^\top \frac{\partial f}{\partial u} \right)}
s\,d\tau +
\int_{t_0}^{t_1}
-\lambda^\top \frac{\partial f}{\partial \theta} d\tau
= 0
\,\,\,\,\,
\forall \lambda(t)
\end{align*}
$$ (eq-ibp)

Notar que  $\left. \lambda s \right|_{t_0}^{t_1} = \lambda(t_1)s(t_1)-\lambda(t_0)s(t_0)$. En general, la condición inicial no depende de $\theta$, por lo que: $\frac{du_0}{d\theta}=0.$ Luego, $s(t_0)=\frac{d u_0}{d \theta} = 0 \Rightarrow \left. \lambda s \right|_{t_0}^{t_1} = \lambda(t_1)s(t_1)$

Y recordemos que en la otra ecuación tenemos

$$
\frac{dL}{d\theta}=
\int_{t_0}^{t_1}
\frac{\partial h}{\partial u}s +
\frac{\partial h}{\partial \theta} \, d\tau
$$ (eq-loss-expanded)

Como vale $\forall \lambda(t)$ seleciono $\lambda(t)$ tal que

$$
{\color{red}
 -\frac{\partial \lambda^\top}{\partial t} -\lambda ^\top \frac{\partial f}{\partial u}
}
=
\frac{\partial h}{\partial u}
$$ (eq-adjoint-ode)

Notemos que {eq}`eq-adjoint-ode` es una ecuacion diferencial para $\lambda(\tau)$, necesitamos condiciones iniciales. Tomemos ${\color{blue}\lambda(t_1)=0}$ como condicion final. Sustituyendo en {eq}`eq-ibp`:

$$
\begin{align*}
0 &=
{\color{blue}0} \cdot s(t_1) +
\int_{t_0}^{t_1}
\frac{\partial h}{\partial u}s \, d\tau +
\int_{t_0}^{t_1}
-\lambda^\top \frac{\partial f}{\partial \theta} \, d\tau
\\
&\Rightarrow
\int_{t_0}^{t_1}
\frac{\partial h}{\partial u}s \, d\tau
=
\int_{t_0}^{t_1}
\lambda^\top \frac{\partial f}{\partial \theta} \, d\tau
\end{align*}
$$

Luego, reemplazo en {eq}`eq-loss-expanded`, y obtengo el gradiente:

$$
\frac{dL}{d\theta}=
\int_{t_0}^{t^1}
\lambda^\top \frac{\partial f}{\partial \theta} +
\frac{\partial h}{\partial \theta} 
\,\, d\tau
$$

:::{dropdown} Observacion: otra forma de verlo.

Consideremos el producto interno entre funciones $\langle f, g \rangle = \int_{t_0}^{t_1} f(t)g(t) dt$.

Definamos:

$$
r(t) = \frac{ds}{dt} \frac{\partial f}{\partial u}s - \frac{\partial f}{\partial\theta}
$$

Luego {eq}`eq-integral-constraint` es simplemente un producto interno:

$$
\langle \lambda, r\rangle = 0 \quad \forall\lambda.
$$

Podemos ver todo esto como una manipulación de productos internos: queremos evitar calcular: $\int \frac{\partial h}{\partial u} {s(t)}$ que es un producto interno.

$$
\langle g, s\rangle = \int_{t_0}^{t_1}\frac{\partial h}{\partial u}s \,d\tau.
$$

Definamos primero

- $g=\frac{\partial h}{ \partial u}$
- $\mathcal{A} = \tfrac{d}{dt} - \tfrac{\partial f}{\partial u}$ un operador lineal (o sea dado  una funcion devuelve funcion y ademas es lineal). Por la encuación de sensibilidad $\mathcal{A}s = b$, con $b=\partial f/\partial\theta$
- Se puede derivar que el operador adjunto es $\mathcal{A}^*\lambda = -\tfrac{d\lambda}{d\tau} - \big(\tfrac{\partial f}{\partial u}\big)^\top\lambda$ ( el operador adjunto se define como $\mathcal A^*$ tal que cumpla: $\int_{t_0}^{t_1} (\mathcal{A}v) \cdot w \, dt = \int_{t_0}^{t_1} v \cdot (\mathcal{A}^*w) \, dt$, puede verse como una generalización de transpuesta $\langle Au, v \rangle = \langle u, A^T v \rangle$ con $A$ matriz, $u,v$ vectores).

Por definición de operador adjunto

$$
\langle \mathcal{A}^*\lambda, s\rangle = \langle \lambda, \mathcal{A}s\rangle
$$

Ahora, si consideramos $\lambda$ como solución de $\mathcal{A}^*\lambda = g$,

$$
\langle g, s\rangle = \langle \mathcal{A}^*\lambda, s\rangle
$$

Usamos que sabemos que $\mathcal{A}s = b$, con $b = \partial f/\partial\theta$

$$
\langle \lambda, \mathcal{A}s\rangle = \langle \lambda, b\rangle.
$$

Luego solo basta computar este producto interno:

$$
\langle \lambda, b\rangle = \int_{t_0}^{t_1}\lambda^\top\frac{\partial f}{\partial\theta}\,d\tau.
$$
:::

### Pasos del método del Adjunto Continuo

De esta forma, obtenemos el siguiente método para computar el gradiente de $dL/d\theta$:

1. Resolver la ODE original (forward): $\frac{du}{dt} = f(u, \theta, t), \quad u(t_0) = u_0$ Se guardan los valores de $u(t)$ o se usan técnicas como checkpointing.
2. Resolver la ecuación adjunta (backward): $\frac{d\lambda}{dt} = - \left(\frac{\partial f}{\partial u}\right)^T \lambda - \left(\frac{\partial h}{\partial u}\right)^T, \quad \lambda(t_1) = 0$ La condición final $\lambda(t_1)=0$ significa que la ODE adjunta se resuelve hacia atrás en el tiempo.
3. Calcular el gradiente: $\frac{dL}{d\theta} = \int_{t_0}^{t_1} \left( \lambda^T \frac{\partial f}{\partial \theta} + \frac{\partial h}{\partial \theta} \right) dt$

#### Checkpointing

Es una técnica para balancear el uso de memoria y el tiempo de cómputo en métodos que requieren almacenar activaciones intermedias (como Reverse AD y el método del adjunto). Consiste en guardar solo algunos puntos intermedios en la memoria y recomputar los demás según sea necesario.

#### Backsolve

Una alternativa al **checkpointing** es resolver la ecuación diferencial al revés. Primero, invierto la variable *temporal* y defino un estado final, en vez de un estado inicial

$$
\begin{align*}
\frac{du}{dt} &= f(u,\theta,t) \\
\overset{t\to-t}{\Rightarrow}\quad
\frac{du}{dt} &= -f(u,\theta,t) \quad \quad u(t_1) = u_1
\end{align*}
$$

Luego, podríamos resolver la ecuación en conjunto con la del adjunto, en modo **reverse** 

$$
\left\{
\begin{aligned}
\frac{du}{dt}
&=
-f(u,\theta,t),
\qquad
u(t_1)=u_1
\\[0.8em]
\frac{d\lambda}{dt}
&=
-\left(\frac{\partial f}{\partial u}\right)^\top \lambda
+
\frac{\partial h}{\partial u}
\end{aligned}
\right.
$$

Esto se denomina **backsolve.** En ciertos casos puede presentar error numérico, pero puede implementarse junto al **checkpointing** para ir corrigiendo dichos errores.