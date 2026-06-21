---
title: No14 - Programación Diferencial Pt4
---

# Programación Diferenciable: Métodos Reverse Pt2

**Fecha:** 08/06/2026

:::{iframe} https://www.youtube.com/embed/OYmbtWjY6HE
:width: 100%
:::
**Autores:** Felipe Cignoli (`@fcignoli`), Martin Sinnona [`@martinsinnona`](https://github.com/martinsinnona), Noé Hsueh [@noehsueh](https://github.com/noehsueh)


## Repaso de la clase anterior: Método del adjunto discreto

Supongamos que la solución discreta del problema se escribe como

```{math}
U=(u_1,\ldots,u_M)\in \mathbb{R}^{nM}.
```

La trayectoria discreta no es una variable libre.
Está determinada por un conjunto de ecuaciones algebraicas que escribimos como

```{math}
G(U,\theta)=0.
```

Aquí $\theta$ representa los parámetros del modelo.
La función de costo que queremos derivar es

```{math}
L=L(U,\theta).
```

Como $U$ depende implícitamente de $\theta$, la derivada total de $L$ es

```{math}
\frac{dL}{d\theta}
=
\frac{\partial L}{\partial U}\frac{\partial U}{\partial \theta}
+
\frac{\partial L}{\partial \theta}.
```

El término difícil es $\partial U/\partial\theta$.
Ese término mide cómo cambia toda la trayectoria numérica cuando cambiamos el parámetro $\theta$.

Para eliminarlo, derivamos la restricción $G(U,\theta)=0$ respecto de $\theta$.
Se obtiene

```{math}
0=
\frac{dG}{d\theta}
=
\frac{\partial G}{\partial U}\frac{\partial U}{\partial \theta}
+
\frac{\partial G}{\partial \theta}.
```

Si $\partial G/\partial U$ es invertible, entonces

```{math}
\frac{\partial U}{\partial \theta}
=
-
\left(\frac{\partial G}{\partial U}\right)^{-1}
\frac{\partial G}{\partial \theta}.
```

Reemplazando en la derivada de $L$, queda

```{math}
\frac{dL}{d\theta}
=
-
\frac{\partial L}{\partial U}
\left(\frac{\partial G}{\partial U}\right)^{-1}
\frac{\partial G}{\partial \theta}
+
\frac{\partial L}{\partial \theta}.
```

El método del adjunto consiste en definir una variable auxiliar $\lambda$ tal que

```{math}
\left(\frac{\partial G}{\partial U}\right)^T\lambda
=
\left(\frac{\partial L}{\partial U}\right)^T.
```

Con esta definición, la derivada total queda

```{math}
\boxed{
\frac{dL}{d\theta}
=
-
\lambda^T\frac{\partial G}{\partial \theta}
+
\frac{\partial L}{\partial \theta}
}
```

Esta expresión permite calcular el gradiente sin construir explícitamente $\partial U/\partial\theta$.

## Ejemplo: solver lineal explícito

Consideremos una ecuación diferencial de la forma

```{math}
\frac{du}{dt}=f(u,\theta,t).
```

En el caso lineal, podemos escribir

```{math}
f(u,\theta,t)=A(\theta,t)u.
```

Usando Euler explícito,

```{math}
\frac{u_{j+1}-u_j}{\Delta t_j}=f(u_j,\theta,t_j).
```

Por lo tanto,

```{math}
u_{j+1}=u_j+\Delta t_j f(u_j,\theta,t_j).
```

En un problema lineal, esto puede escribirse como

```{math}
u_{j+1}=A_j(\theta)u_j+b_j(\theta).
```

En esta notación, $A_j$ representa la matriz de avance del paso temporal $j$.
Si no hay término inhomogéneo, entonces $b_j=0$.

El residuo discreto de cada paso es

```{math}
g_j(U,\theta)=u_{j+1}-A_j(\theta)u_j-b_j(\theta)=0.
```

Apilando todos los pasos temporales, el sistema completo puede escribirse como

```{math}
G(U,\theta)\equiv \mathcal{A}(\theta)U-B(\theta)=0.
```

La matriz $\mathcal{A}$ tiene estructura triangular por bloques.
Esquemáticamente,

```{math}
\begin{pmatrix}
I & 0 & 0 & \cdots & 0 \\
-A_0 & I & 0 & \cdots & 0 \\
0 & -A_1 & I & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & -A_{M-1} & I
\end{pmatrix}
\begin{pmatrix}
u_0 \\
u_1 \\
u_2 \\
\vdots \\
u_M
\end{pmatrix}
=
\begin{pmatrix}
u_0^{\mathrm{ini}} \\
b_0 \\
b_1 \\
\vdots \\
b_{M-1}
\end{pmatrix}.
```

Con esta notación,

```{math}
\frac{\partial G}{\partial \theta}
=
\frac{\partial \mathcal{A}}{\partial\theta}U
-
\frac{\partial B}{\partial\theta}.
```

Entonces el gradiente se calcula como

```{math}
\frac{dL}{d\theta}
=
-
\lambda^T
\left(
\frac{\partial \mathcal{A}}{\partial\theta}U
-
\frac{\partial B}{\partial\theta}
\right)
+
\frac{\partial L}{\partial\theta}.
```

La matriz que aparece en la ecuación adjunta es $\mathcal{A}^T$.
Por eso, aunque el problema directo se resuelve hacia adelante en el tiempo, el problema adjunto se resuelve hacia atrás.

## Ejemplo de función de costo

Una función de costo típica para comparar la simulación con datos observados es

```{math}
L(U,\theta)=
\sum_{j=1}^{M}
w_j
\left\|u_j-u_j^{\mathrm{obs}}\right\|_2^2.
```

Si se usa un factor $1/2$ delante de la norma cuadrática, la derivada queda sin el factor $2$.
Equivalentemente, ese factor puede absorberse en los pesos $w_j$.

Con esa convención, el término fuente de la ecuación adjunta es

```{math}
\frac{\partial L}{\partial u_j}
=
w_j\left(u_j-u_j^{\mathrm{obs}}\right).
```

La condición final para el adjunto es

```{math}
\lambda_M
=
w_M\left(u_M-u_M^{\mathrm{obs}}\right).
```

La recurrencia hacia atrás es

```{math}
\lambda_j
=
A_j^T\lambda_{j+1}
+
w_j\left(u_j-u_j^{\mathrm{obs}}\right),
\qquad
j=M-1,\ldots,1.
```

Esta ecuación se resuelve en modo reverso.
Primero se calcula y se guarda la trayectoria directa $u_0,u_1,\ldots,u_M$.
Luego se calcula $\lambda_M$.
Finalmente se propaga $\lambda_j$ hacia atrás.

## Algoritmo práctico

El procedimiento práctico es el siguiente.

1. Resolver el problema directo hacia adelante en el tiempo.
2. Guardar la trayectoria $u_j$.
3. Calcular la condición final del adjunto a partir de la función de costo.
4. Resolver la ecuación adjunta hacia atrás en el tiempo.
5. Usar $\lambda$ para calcular el gradiente respecto de los parámetros.

## Comentarios conceptuales

El adjunto discreto evita calcular una sensibilidad distinta para cada parámetro.
Esto es especialmente útil cuando hay muchos parámetros y una única función de costo escalar.

Aunque la ecuación diferencial original sea no lineal, la ecuación para $\lambda$ es lineal en $\lambda$.
La linealidad aparece porque el adjunto se obtiene al linearizar alrededor de la trayectoria directa ya calculada.

En muchos casos, el método del adjunto discreto es equivalente a hacer backpropagation sobre el solver numérico.
Por eso se dice que el adjunto se resuelve en modo reverso.




## Método del Adjunto Continuo

Consideremos una ODE de primer orden dada por

$$
\frac{du}{dt} = f(u,\theta,t)
$$ (eq-ode)

sujeta a la condición inicial $u(t_0)=u_0$, donde $u \in \mathbb{R}^n$ es el vector solución desconocido de la ODE, $f:\mathbb{R}^n \times \mathbb{R}^p \times \mathbb{R} \to \mathbb{R}^n$ es una función que depende de: el estado $u$, $\theta \in \mathbb{R}^p$ es un vector de parámetros, y $t \in [t_0,t_1]$ se refiere al tiempo. Aquí, $n$ denota el tamaño de la ODE y $p$ el número de parámetros. Resolver la ODE implica obtener $u(t)$, que depende de $\theta$. En general no es posible obtener una solución explícita de $u$ (salvo en casos lineales o muy particulares), por lo que debemos resolverla numéricamente.

Recordemos que queremos obtener $\theta$ (donde $\theta$ pueden ser los parámetros de una red o, en problemas inversos, coeficientes de ecuaciones diferenciales). De esta forma, nos interesa generalmente $\frac{dL}{d\theta}$.

Para ello, veamos qué es $L$. Podemos escribir el término de la *loss* de forma general como una integral[^loss-as-integral]:
[^loss-as-integral]:¿Por qué podemos escribirlo como una integral?

    Veamos el ejemplo del caso discreto, donde la *loss* es una suma ponderada $\sum_i w_i \|u(t_i, \theta)-u^{\text{obs}}_i\|_2^2$ sobre instantes de observación $t_i$. Podemos escribir el integrando como

    $$
    h(u, \theta, t) = \sum_{i} w_i\,\|u(t; \theta) - u_i^{\text{obs}}\|^2 \,\delta(t - t_i),
    $$

    donde $\delta$ es la función delta de Dirac. Si tomamos la integral de $h$ de $t_0$ a $t_1$, recupera exactamente la suma.


$$
L(u(\cdot, \theta);\theta)=\int_{t_0}^{t_1} h(u(\tau;\theta),\theta,\tau)\,d\tau
$$ (eq-loss)

donde $h$ es una función de costo puntual (evaluada en cada instante $\tau$).




Ahora, derivemos {eq}`eq-loss` con respecto de $\theta$ usando la regla de la cadena (notar que $h$ depende de $\theta$ tanto de forma directa como a través de $u(\tau;\theta)$):

$$
\frac{dL}{d\theta} = \int_{t_0}^{t_1} \left( \frac{\partial h}{\partial \theta} + \frac{\partial h}{\partial u} \underbrace{{\color{red}\frac{\partial u}{\partial\theta}}}_{{\color{red}s(t)}} \right) dt.
$$ (eq-dloss)

En {eq}`eq-dloss` notamos que aparece un VJP[^shapes] $\frac{\partial h}{\partial u}\,{\color{red}s(t)}$, con **sensibilidad** ${\color{red}s(t)} = \frac{\partial u}{\partial \theta} \in \mathbb{R}^{n\times p}$[^sensibilidad]. Como en el método discreto, la idea *del método adjunto consiste en* aprovechar esta estructura para introducir una nueva variable (*adjunto* ${\color{blue}\lambda}$) que nos permita **evitar calcular** el jacobiano ${\color{red}s(t)}$.

[^shapes]: Notar que $\partial h/\partial u$ es de tamaño $1\times n$ (ya que $h$ es una función escalar y $u\in\mathbb{R}^n$, su gradiente es un vector fila de $n$ componentes), de modo que el producto $\frac{\partial h}{\partial u}\,{\color{red}s(t)}$ es un vector de $1\times p$.


¿Por qué es costoso ${\color{red}s(t)}$? Es una matriz de $n\times p$ y, como veremos, satisface su propia ODE. Cuando $p$ es grande (por ejemplo, los parámetros de una red), esto se vuelve prohibitivo, y de ahí la motivación del método adjunto.

[^sensibilidad]: ${\color{red}s(t)}$ define qué tanto cambia mi solución $u(t)\in\mathbb{R}^n$ con respecto de $\theta$: $\frac{\partial u }{\partial \theta}$. Notemos que tiene su ecuación diferencial asociada. Diferenciemos {eq}`eq-ode` con respecto de $\theta$:

    $$
    \frac{d}{d\theta} \left[ \frac{d}{dt}u(t; \theta) \right] = \frac{d}{d\theta} \left[ f(u(t; \theta), \theta, t) \right] \tag{i}
    $$

    Ahora intercambiamos el orden de las derivadas parciales (asumimos que vale: solución suave, etc.):

    $$
    \frac{d}{dt} \left[ \frac{du}{d\theta} \right] = \frac{d}{d\theta} f(u(t; \theta), \theta, t) \tag{ii}
    $$

    En el lado derecho desarrollamos la derivada con la regla de la cadena ($f$ depende de $\theta$ directamente y vía $u$):

    $$
    \frac{d}{dt} \left[ \frac{du}{d\theta} \right] = \frac{\partial f}{\partial u} \frac{\partial u}{\partial \theta} + \frac{\partial f}{\partial \theta} \tag{iii}
    $$

    Definamos ahora la matriz de sensibilidad como $s(t) = \frac{\partial u}{\partial \theta}$ y reemplazamos, obteniendo así la **ecuación de sensibilidad**:

    $$
    \frac{ds}{dt} = \frac{\partial f}{\partial u} s + \frac{\partial f}{\partial \theta} \tag{iv}
    $$

    Su condición inicial está dada por la derivada de $u_0$:

    $$
    s(t_0) = \frac{du(t_0)}{d\theta} \tag{v}
    $$

    y si $u_0$ no depende de $\theta$, entonces $s(t_0)=0$. De esta forma, la ecuación de sensibilidad me dice cómo un cambio de los parámetros afecta a mi solución del sistema en el tiempo.

Notemos que la sensibilidad tiene su ecuación diferencial asociada:

$$
\frac{ds}{dt} = \frac{\partial f}{\partial u} s + \frac{\partial f}{\partial \theta}
$$

Reordenamos los términos para dejarla igualada a cero:

$$
\left[ \frac{ds}{dt} = \frac{\partial f}{\partial u} s + \frac{\partial f}{\partial \theta} \right]
\Rightarrow \frac{ds}{dt} - \frac{\partial f}{\partial u} s - \frac{\partial f}{\partial \theta} = 0
$$

Esta expresión es cero para cualquier instante de tiempo, así que podemos multiplicarla por ${\color{blue}\lambda(t)^\top}$ e integrarla, y va a seguir siendo cero:

$$
\begin{align*}
\Rightarrow \int_{t_0}^{t_1} {\color{blue}\lambda(\tau)^\top} \left[ \frac{d{\color{red}s}}{dt} - \frac{\partial f}{\partial u} {\color{red}s} - \frac{\partial f}{\partial \theta} \right] d\tau &= 0\quad \forall\, {\color{blue}\lambda(t)}:[t_0, t_1] \mapsto \mathbb{R}^n
\end{align*}
$$ (eq-integral-constraint)

Recordemos que el objetivo ahora es **eliminar** la sensibilidad ${\color{red}s(t)}$. Primero vamos a usar integración por partes sobre ${\color{blue}\lambda^\top} \frac{d{\color{red}s}}{dt}$ para trasladar la derivada temporal de ${\color{red}s}$ hacia ${\color{blue}\lambda}$, y luego reemplazamos en {eq}`eq-integral-constraint`.

:::{margin}
Recordemos la **integración por partes**:

$$
\int_{t_0}^{t_1} {\color{blue}\lambda^\top} \frac{d{\color{red}s}}{dt}\,d\tau = \left.{\color{blue}\lambda^\top}{\color{red}s}\right|_{t_0}^{t_1} - \int_{t_0}^{t_1} \frac{d{\color{blue}\lambda^\top}}{dt}{\color{red}s}\,d\tau.
$$
:::

$$
\begin{align*}
0 &= \int_{t_0}^{t_1}\Bigg[\;
\overbrace{{\color{olive}\lambda^\top\,\frac{d{s}}{dt}}}^{\textstyle\text{aplicamos partes}}
\;-\;{\color{black}\lambda^\top}\,\frac{\partial f}{\partial u}\,{\color{black}s}
\;-\;{\color{black}\lambda^\top}\,\frac{\partial f}{\partial \theta}\;\Bigg]\,d\tau \\[1.4em]
&= \underbrace{{\color{olive}\left.{\color{black}\lambda^\top}{s}\,\right|_{t_0}^{t_1}
\;-\int_{t_0}^{t_1}\frac{d{\lambda^\top}}{dt}\,{s}\,d\tau}}_{\text{partes}}
\;-\;\int_{t_0}^{t_1}{\lambda^\top}\,\frac{\partial f}{\partial u}\,{s}\,d\tau
\;-\;\int_{t_0}^{t_1}{\lambda^\top}\,\frac{\partial f}{\partial \theta}\,d\tau \\[1.4em]
&= \underbrace{{\lambda(t_1)^\top}\,{s(t_1)}}_{\substack{\text{frontera con }{s(t_0)}=0}}
\;+\;\int_{t_0}^{t_1}\underbrace{{\color{orange}\boxed{-\dfrac{d{\lambda^\top}}{dt}-{\lambda^\top}\,\dfrac{\partial f}{\partial u}}}}_{\text{coeficiente de }{\color{red}s}}\;{\color{red}s}\,d\tau
\;-\;\int_{t_0}^{t_1}{\lambda^\top}\,\frac{\partial f}{\partial \theta}\,d\tau
\qquad \forall\,{\color{blue}\lambda(t)}
\end{align*}
$$ (eq-ibp)

Notar que $\left. {\lambda^\top}{s} \right|_{t_0}^{t_1} = {\lambda(t_1)^\top}{s(t_1)}-{\lambda(t_0)^\top}{s(t_0)}$. En general, la condición inicial no depende de $\theta$, por lo que $\frac{du_0}{d\theta}=0$. Luego, ${s(t_0)}=\frac{d u_0}{d \theta} = 0 \Rightarrow \left. {\lambda^\top}{s} \right|_{t_0}^{t_1} = {\lambda(t_1)^\top}{s(t_1)}$.

Y recordemos que en la ecuación {eq}`eq-dloss` teníamos:

$$
\frac{dL}{d\theta}=
\int_{t_0}^{t_1}
{\color{orange}
\frac{\partial h}{\partial u}{\color{red}s} }+
\frac{\partial h}{\partial \theta} \, d\tau
$$ (eq-loss-expanded)

Como {eq}`eq-ibp` vale $\forall\, {\color{blue}\lambda(t)}$, podemos **elegir** ${\color{blue}\lambda(t)}$ de forma inteligente. Seleccionamos ${\color{blue}\lambda(t)}$ tal que el coeficiente de ${\color{red}s}$ coincida con $\frac{\partial h}{\partial u}$:

$$
{
 -\frac{d \lambda^\top}{dt} -\lambda ^\top \frac{\partial f}{\partial u}
}
=
\frac{\partial h}{\partial u}
$$ (eq-adjoint-ode)

Notemos que {eq}`eq-adjoint-ode` es una ecuación diferencial para ${\color{blue}\lambda(\tau)}$; necesitamos una condición para resolverla. Tomemos ${\color{blue}\lambda(t_1)=0}$ como **condición final**. Sustituyendo ambas elecciones en {eq}`eq-ibp`:

$$
\begin{align*}
0 &=
{\color{blue}\underbrace{\lambda(t_1)^\top}_{=\,0}}\, {\color{red}s(t_1)} +
\int_{t_0}^{t_1}
\frac{\partial h}{\partial u}{\color{red}s} \, d\tau -
\int_{t_0}^{t_1}
{\color{blue}\lambda^\top} \frac{\partial f}{\partial \theta} \, d\tau
\\
&\Rightarrow
\int_{t_0}^{t_1}
\frac{\partial h}{\partial u}{\color{red}s} \, d\tau
=
\int_{t_0}^{t_1}
{\color{blue}\lambda^\top} \frac{\partial f}{\partial \theta} \, d\tau
\end{align*}
$$

Luego, reemplazamos en {eq}`eq-loss-expanded` y obtenemos el gradiente, ya **sin** la sensibilidad ${\color{red}s}$:

$$
\frac{dL}{d\theta}=
\int_{t_0}^{t_1}
{\color{blue}\lambda^\top} \frac{\partial f}{\partial \theta} +
\frac{\partial h}{\partial \theta}
\,\, d\tau
$$

Finalmente, transponiendo {eq}`eq-adjoint-ode` obtenemos la forma estándar de la **ecuación adjunta** 
:::{margin}
Recordando que $(\lambda^\top A)^\top = A^\top \lambda$:
:::

$$
\frac{d{\color{blue}\lambda}}{dt} = - \left(\frac{\partial f}{\partial u}\right)^\top {\color{blue}\lambda} - \left(\frac{\partial h}{\partial u}\right)^\top, \qquad {\color{blue}\lambda(t_1)=0}.
$$

:::{dropdown} Observación: otra forma de verlo (operadores adjuntos)

Consideremos el producto interno entre funciones $\langle f, g \rangle = \int_{t_0}^{t_1} f(t)^\top g(t)\, dt$.

Definamos el residuo de la ecuación de sensibilidad:

$$
r(t) = \frac{ds}{dt} - \frac{\partial f}{\partial u}s - \frac{\partial f}{\partial\theta}
$$

Luego {eq}`eq-integral-constraint` es simplemente un producto interno nulo:

$$
\langle \lambda, r\rangle = 0 \quad \forall\lambda.
$$

Podemos ver todo esto como una manipulación de productos internos: queremos evitar calcular $\int \frac{\partial h}{\partial u}\, s(t)$, que también es un producto interno:

$$
\langle g, s\rangle = \int_{t_0}^{t_1}\frac{\partial h}{\partial u}\,s \,d\tau.
$$

Definamos primero:

- $g=\left(\frac{\partial h}{ \partial u}\right)^\top$
- $\mathcal{A} = \tfrac{d}{dt} - \tfrac{\partial f}{\partial u}$, un operador lineal (dada una función devuelve una función, y además es lineal). Por la ecuación de sensibilidad, $\mathcal{A}s = b$, con $b=\partial f/\partial\theta$.
- Se puede derivar que el operador adjunto es $\mathcal{A}^*\lambda = -\tfrac{d\lambda}{d\tau} - \big(\tfrac{\partial f}{\partial u}\big)^\top\lambda$. El operador adjunto $\mathcal{A}^*$ se define como aquel que cumple $\int_{t_0}^{t_1} (\mathcal{A}v)^\top w \, dt = \int_{t_0}^{t_1} v^\top (\mathcal{A}^*w) \, dt$; puede verse como una generalización de la transpuesta, $\langle Au, v \rangle = \langle u, A^\top v \rangle$, con $A$ matriz y $u,v$ vectores.

Por definición de operador adjunto:

$$
\langle \mathcal{A}^*\lambda, s\rangle = \langle \lambda, \mathcal{A}s\rangle.
$$

Ahora, si consideramos $\lambda$ como solución de $\mathcal{A}^*\lambda = g$,

$$
\langle g, s\rangle = \langle \mathcal{A}^*\lambda, s\rangle.
$$

Usamos que $\mathcal{A}s = b$, con $b = \partial f/\partial\theta$:

$$
\langle \lambda, \mathcal{A}s\rangle = \langle \lambda, b\rangle.
$$

Luego, solo basta computar este producto interno:

$$
\langle \lambda, b\rangle = \int_{t_0}^{t_1}\lambda^\top\frac{\partial f}{\partial\theta}\,d\tau.
$$
:::

### Pasos del método del Adjunto Continuo

De esta forma, obtenemos el siguiente método para computar el gradiente $dL/d\theta$:

1. **Resolver la ODE original (forward):** $\dfrac{du}{dt} = f(u, \theta, t), \quad u(t_0) = u_0$. Se guardan los valores de $u(t)$ o se usan técnicas como *checkpointing*.
2. **Resolver la ecuación adjunta (backward):** $\dfrac{d\lambda}{dt} = - \left(\dfrac{\partial f}{\partial u}\right)^\top \lambda - \left(\dfrac{\partial h}{\partial u}\right)^\top, \quad \lambda(t_1) = 0$. La condición final $\lambda(t_1)=0$ significa que la ODE adjunta se resuelve **hacia atrás** en el tiempo (de $t_1$ a $t_0$).
3. **Calcular el gradiente:** $\dfrac{dL}{d\theta} = \displaystyle\int_{t_0}^{t_1} \left( \lambda^\top \dfrac{\partial f}{\partial \theta} + \dfrac{\partial h}{\partial \theta} \right) dt$.

#### Checkpointing

Es una técnica para balancear el uso de memoria y el tiempo de cómputo en métodos que requieren almacenar activaciones intermedias (como Reverse AD y el método del adjunto). Consiste en guardar solo algunos puntos intermedios en memoria y recomputar los demás según sea necesario, intercambiando memoria por cómputo.

#### Backsolve

Una alternativa al **checkpointing** es no almacenar la trayectoria $u(t)$, sino reconstruirla resolviendo la ODE hacia atrás junto con la del adjunto. Primero invertimos la variable *temporal* y definimos un estado final, en vez de un estado inicial:

$$
\begin{align*}
\frac{du}{dt} &= f(u,\theta,t) \\
\overset{t\to-t}{\Rightarrow}\quad
\frac{du}{dt} &= -f(u,\theta,t) \quad \quad u(t_1) = u_1
\end{align*}
$$

Bajo el cambio $t\to -t$, todo lado derecho cambia de signo. Aplicándolo también a la ecuación adjunta (cuya forma estándar es $\frac{d\lambda}{dt} = -(\partial f/\partial u)^\top \lambda - (\partial h/\partial u)^\top$), podemos resolver el sistema acoplado hacia atrás, en modo **reverse**:

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
\left(\frac{\partial f}{\partial u}\right)^\top \lambda
+
\left(\frac{\partial h}{\partial u}\right)^\top,
\qquad
\lambda(t_1)=0
\end{aligned}
\right.
$$

Esto se denomina **backsolve**. Su ventaja es que evita almacenar la trayectoria completa (poca memoria); su desventaja es que, en ciertos casos, reconstruir $u$ hacia atrás puede acumular error numérico. Por eso puede combinarse con **checkpointing** para reanclar la solución en puntos guardados e ir corrigiendo dichos errores.
