---
title: No14 - Programación Diferencial Pt4
---

# Programación Diferenciable: Métodos Reverse Pt2

**Fecha:** 08/06/2026

:::{iframe} https://www.youtube.com/embed/OYmbtWjY6HE
:width: 100%
:::
**Autores:** Felipe Cignoli (`@fcignoli`)


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


