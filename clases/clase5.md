---
title: No5 - Optimización
---

**Fecha:** 22/04/2026


:::{iframe} https://www.youtube.com/embed/SKOspQoy9N4
:width: 100%
:::

Estas notas introducen el problema general de optimización que aparece al ajustar modelos dinámicos a datos. En particular, veremos cómo se define una función de costo, qué significa buscar parámetros óptimos, qué diferencias hay entre búsqueda global y búsqueda local, y cómo se conectan los métodos de optimización con el ajuste de trayectorias de ODEs.

# Problema general de optimización

En muchos problemas de modelado queremos encontrar parámetros que hagan que un modelo describa lo mejor posible los datos observados. Para eso definimos una **función de costo** o **función de pérdida** $L$, que mide qué tan mal se ajusta el modelo para un cierto valor de los parámetros.

De manera general, buscamos resolver el problema

$$
\min_{\theta} L(\theta; X, Y), \qquad \theta \in \Omega \subseteq \mathbb{R}^p.
$$

Donde:

* $\theta$ es el vector de parámetros del modelo.
* $p$ es la cantidad de parámetros.
* $\Omega$ es el espacio de búsqueda o conjunto de valores admisibles para $\theta$.
* $X$ e $Y$ representan los datos disponibles.
* $L(\theta; X,Y)$ mide el desacuerdo entre el modelo y los datos.

El objetivo es encontrar un parámetro $\theta^*$ que minimice la función de costo:

$$
\theta^* = \arg\min_{\theta \in \Omega} L(\theta).
$$

:::{note} Estimación paramétrica
En este contexto, la optimización aparece como una herramienta para realizar **estimación paramétrica**: buscamos parámetros que no conocemos, pero que queremos inferir a partir de datos observados.
:::

En la práctica, incluso si escribimos el problema como una minimización exacta, no siempre podemos encontrar el mínimo global de forma analítica. Por eso necesitamos algoritmos numéricos de optimización.

# ¿De dónde sale la función de costo?

Una pregunta central es cómo construir la función $L(\theta)$. En general, la función de costo se elige para cuantificar el error entre las observaciones y las predicciones del modelo.

Por ejemplo, si tenemos observaciones $Y_i$ y un modelo que predice valores $\hat{Y}_i(\theta)$, una elección clásica es la suma de errores cuadráticos:

$$
L(\theta) = \sum_{i=1}^N \|Y_i - \hat{Y}_i(\theta)\|_2^2.
$$

Esta elección penaliza más fuertemente errores grandes y conduce al método de **cuadrados mínimos**. En clases anteriores apareció en el contexto de ajuste de trayectorias: se comparan datos observados contra trayectorias generadas por una ODE.

:::{tip} Interpretación estadística(próxima clase visto más a detalle)
Cuando se minimiza una suma de cuadrados, muchas veces se está asumiendo implícitamente un modelo de ruido Gaussiano independiente con varianza constante. Bajo ese supuesto, minimizar cuadrados mínimos coincide con maximizar la verosimilitud.
:::

# Tipos de algoritmos de optimización

Podemos clasificar los algoritmos de optimización según la información que usan de la función de costo y según la región del espacio de parámetros que exploran.

Una primera distinción importante es entre:

* **Búsqueda global:** intenta explorar ampliamente el espacio $\Omega$ para encontrar el mejor mínimo posible.
* **Búsqueda local:** parte de un valor inicial $\theta_0$ y mejora iterativamente la solución en una vecindad.

## Búsqueda global

La búsqueda global intenta encontrar el mínimo global de $L(\theta)$ en todo el espacio de parámetros $\Omega$.

Una estrategia simple consiste en evaluar la función de costo en muchos puntos del espacio de parámetros, por ejemplo mediante una grilla. En dimensión baja esto puede ser razonable, pero en dimensión alta se vuelve rápidamente inviable.

:::{warning} Maldición de la dimensionalidad
Si cada parámetro se evalúa en $k$ posibles valores y tenemos $p$ parámetros, una grilla completa requiere $k^p$ evaluaciones de la función de costo. Por eso la búsqueda exhaustiva escala muy mal con la dimensión.
:::

La búsqueda global puede ser útil cuando el espacio de parámetros es pequeño, cuando la función tiene muchos mínimos locales o cuando no tenemos una buena inicialización. Sin embargo, suele ser costosa.
Para tener una idea mas concreta de cómo escala este método con la dimensión, en clase vimos que este método de es útil hasta 4 o 5 dimensiones. Luego rápidamente se vuelve incalculable




## Búsqueda local

La búsqueda local comienza desde un punto inicial $\theta_0$ y genera una sucesión de parámetros

$$
\theta_0, \theta_1, \theta_2, \ldots
$$

buscando que la función de costo disminuya en cada paso.

La idea general es actualizar los parámetros de la forma

$$
\theta^{m+1} = \theta^m + \Delta \theta^m,
$$

para alguna dirección de actualización $\Delta \theta^m$.

A diferencia de la búsqueda global, estos métodos no garantizan necesariamente encontrar el mínimo global. Pueden converger a un mínimo local, a un punto silla, o incluso fallar si la función está mal condicionada o si la inicialización es mala.

:::{note}
En problemas de ajuste de modelos dinámicos, muchas veces se combinan estrategias: primero se usa una búsqueda global o una inicialización heurística, y luego se refina la solución con un método local.
:::

# Métodos locales según la información que usan

Dentro de los métodos locales, podemos distinguir tres grandes familias:

1. **Métodos de orden cero:** sólo evalúan $L(\theta)$.
2. **Métodos de primer orden:** usan $L(\theta)$ y su gradiente $\nabla_\theta L(\theta)$.
3. **Métodos de segundo orden:** usan $L(\theta)$, el gradiente y la información de curvatura, usualmente mediante la matriz Hessiana.

## Métodos de orden cero

Los métodos de orden cero no usan derivadas. Sólo necesitan evaluar la función de costo en distintos valores de $\theta$.

Esto puede ser útil cuando:

* La función no es diferenciable.
* Calcular derivadas es demasiado costoso.
* El modelo funciona como una caja negra.

Sin embargo, al no usar información geométrica de la función, suelen necesitar muchas evaluaciones de $L$.

## Métodos de primer orden

Los métodos de primer orden usan el gradiente de la función de costo:

$$
\nabla_\theta L(\theta) = \frac{\partial L}{\partial \theta}.
$$

El gradiente indica la dirección de mayor crecimiento local de la función. Por eso, para minimizar, nos movemos en la dirección opuesta al gradiente.

# Descenso por gradiente

El método más básico de primer orden es el **descenso por gradiente**. La actualización es

$$
\theta^{m+1} = \theta^m - \alpha^m \nabla_\theta L(\theta^m), \qquad m=0,1,2,\ldots
$$

Donde $\alpha^m > 0$ es el tamaño de paso o *learning rate* en la iteración $m$.

La intuición es simple: si el gradiente apunta hacia donde $L$ aumenta más rápido, entonces $-\nabla_\theta L$ apunta hacia donde $L$ disminuye más rápido localmente.

:::{warning} Tamaño de paso
Si $\alpha^m$ es demasiado grande, el método puede oscilar o divergir. Si es demasiado chico, puede converger muy lentamente. Elegir el tamaño de paso es una parte central del problema de optimización.
:::

## Momentum

Una mejora común del descenso por gradiente es agregar **momentum**. La idea es que la dirección de actualización no dependa solamente del gradiente actual, sino también de las direcciones tomadas en pasos anteriores.

Una forma de escribirlo es:

$$
\theta^{m+1} = \theta^m - \alpha^m g^m,
$$

con

$$
g^m = \eta g^{m-1} + (1-\eta) \nabla_\theta L(\theta^m),
$$

para $0 \leq \eta < 1$.

El término $g^m$ actúa como una media móvil de gradientes. Esto puede suavizar oscilaciones y acelerar el avance en direcciones consistentes.



## Forma general de los métodos de primer orden

Una forma más general de escribir estos algoritmos es

$$
\theta^{m+1} = \theta^m - \alpha^m g^m,
$$

con

$$
g^m = g\left(\nabla_\theta L(\theta^m), \nabla_\theta L(\theta^{m-1}), \ldots, \theta^m, \theta^{m-1}, \ldots\right).
$$

Es decir, la dirección de actualización puede depender del gradiente actual, de gradientes anteriores y de parámetros anteriores.

:::{note} Adam
El optimizador Adam entra dentro de esta familia. Combina ideas de momentum con estimaciones adaptativas de escala para cada coordenada del gradiente. Por eso es muy usado en entrenamiento de redes neuronales.
:::

# Métodos de segundo orden

Los métodos de segundo orden usan información sobre la curvatura de la función de costo. Esta información está contenida en la matriz Hessiana:

$$
H_L(\theta) = \nabla^2_\theta L(\theta).
$$

Mientras que el gradiente indica una dirección de descenso local, el Hessiano describe cómo cambia ese gradiente alrededor del punto actual.

## Método de Newton

El método de Newton usa una aproximación cuadrática local de la función de costo. Su actualización puede escribirse como

$$
\theta^{m+1}
= \theta^m - \alpha^m \left(H_L(\theta^m)\right)^{-1} \nabla_\theta L(\theta^m).
$$

La diferencia principal con descenso por gradiente es que la dirección de descenso se corrige usando la curvatura local de la función.

:::{important}
Cuando el Hessiano está bien definido y es positivo definido cerca del óptimo, Newton puede converger muy rápido. Sin embargo, calcular, almacenar e invertir el Hessiano puede ser muy costoso en problemas de alta dimensión.
:::

## Métodos quasi-Newton

Los métodos quasi-Newton buscan aproximar la información de segundo orden sin calcular explícitamente el Hessiano exacto.

En lugar de usar $H_L(\theta)$ directamente, construyen una aproximación a partir de gradientes y pasos anteriores. Un ejemplo muy usado es **BFGS**.

:::{note} BFGS
BFGS estima progresivamente una matriz que aproxima la curvatura de la función de costo. Esto permite capturar parte de la ventaja de Newton sin pagar el costo completo de calcular el Hessiano exacto.
:::

# Optimización en ajuste de trayectorias de ODEs

Volvamos ahora al caso que aparece en el curso: queremos ajustar parámetros de una ODE a partir de datos observados.

Supongamos que tenemos un sistema dinámico

$$
\frac{dx}{dt} = f(x,t,\theta),
\qquad
x(t_0) = x_0.
$$

Para cada valor de $\theta$, resolver la ODE produce una trayectoria $x(t;\theta)$. Si observamos datos $Y_i$ en tiempos $t_i$, una función de pérdida natural es

$$
L(\theta) = \sum_{i=1}^N \|Y_i - x(t_i;\theta)\|_2^2.
$$

Como $x(t_i;\theta)$ se obtiene resolviendo numéricamente una ODE, también podemos escribir esta pérdida como

$$
L(\theta) = \sum_{i=1}^N \|Y_i - \operatorname{ODESolve}(t_i,\theta)\|_2^2.
$$

:::{note}
En este caso, la función de costo no es una fórmula cerrada simple en $\theta$. Cada evaluación de $L(\theta)$ requiere resolver la ODE para esos parámetros.
:::

Desde el punto de vista de la optimización, el problema sigue siendo

$$
\min_{\theta} L(\theta),
$$

pero la dependencia de $L$ respecto de $\theta$ está mediada por la solución de una ecuación diferencial.

# Gradiente de la pérdida en ODEs

Para usar métodos de primer orden necesitamos calcular

$$
\frac{dL}{d\theta}.
$$

Si

$$
L(\theta) = \sum_{i=1}^N \|Y_i - x(t_i;\theta)\|_2^2,
$$

entonces, aplicando regla de la cadena, obtenemos una expresión del tipo

$$
\frac{dL}{d\theta}
= -\sum_{i=1}^N 2\left(Y_i - x(t_i;\theta)\right)^\top
\frac{\partial x(t_i;\theta)}{\partial \theta}.
$$

El término clave es

$$
S(t_i,\theta) = \frac{\partial x(t_i;\theta)}{\partial \theta},
$$

que se conoce como **sensibilidad** de la solución respecto de los parámetros.

:::{important} Sensibilidades
Las sensibilidades miden cuánto cambia la trayectoria $x(t;\theta)$ cuando perturbamos los parámetros $\theta$. Son necesarias para calcular gradientes de la pérdida respecto de los parámetros.
:::

Si $x(t;\theta) \in \mathbb{R}^n$ y $\theta \in \mathbb{R}^p$, entonces

$$
S(t,\theta) \in \mathbb{R}^{n \times p}.
$$

Cada columna de $S$ representa cómo cambia el estado del sistema cuando se modifica uno de los parámetros.

En la practica NO es necesario calcular esta matriz S(el profe no explicó bien cómo), ya que calcular esa matriz es muy costoso.


# Cierre de la clase

La optimización es el puente entre el modelo dinámico y los datos. En el contexto de ODEs, cada valor de los parámetros define una trayectoria; la función de costo mide qué tan lejos está esa trayectoria de las observaciones; y el algoritmo de optimización actualiza los parámetros para reducir ese desacuerdo.

En resumen:

* Definimos una función de costo $L(\theta)$.
* Buscamos un parámetro óptimo $\theta^*$ que minimice esa función.
* Podemos usar búsqueda global o búsqueda local.
* Los métodos locales pueden usar sólo evaluaciones de $L$, gradientes o información de segundo orden.
* En ODEs, evaluar la pérdida implica resolver numéricamente el sistema.
* Para optimizar eficientemente, necesitamos entender cómo cambia la solución de la ODE respecto de los parámetros.

---

# Implementación computacional en Julia: optimización de una UDE

A continuación presentamos un ejemplo computacional cuyo objetivo principal es mostrar cómo aparece la **optimización** cuando ajustamos una dinámica parcialmente desconocida. 

El código completo de este ejemplo se encuentra disponible en [`05_LV_inverse_UDE`](https://github.com/facusapienza21/DM2026-Curso/tree/main/code/05_LV_inverse_UDE).

Para una explicación previa del modelo Lotka-Volterra y del problema dinámico que se busca resolver, se puede consultar la [`Clase_2`](https://github.com/facusapienza21/DM2026-Curso/blob/main/clases/clase2.md).
En este ejemplo, partimos de datos generados por un sistema Lotka-Volterra conocido y entrenamos una **Universal Differential Equation** (UDE) para reproducir esa trayectoria. La parte central del problema es estimar los parámetros de una red neuronal que aparece dentro de la ecuación diferencial.

Desde el punto de vista de la optimización, queremos resolver un problema de la forma

$$
\min_\theta L(\theta),
$$

donde $\theta$ representa los pesos y sesgos de la red neuronal, y $L(\theta)$ mide qué tan lejos está la trayectoria predicha por el modelo de la trayectoria observada.

---

## ¿Qué se está optimizando?

En una ODE clásica, los parámetros suelen ser coeficientes explícitos del modelo. Por ejemplo, en Lotka-Volterra podríamos querer estimar parámetros como $\alpha$, $\beta$, $\delta$ o $\gamma$.

En este ejemplo, en cambio, parte de la dinámica se reemplaza por una red neuronal. La UDE tiene la forma general:

$$
\frac{du}{dt}
=
f_{\text{conocida}}(u)
+
f_{\text{NN}}(u;\theta),
$$

donde:

* $u(t) = (x(t), y(t))$ es el estado del sistema.
* $f_{\text{conocida}}$ representa la parte de la dinámica que dejamos fija.
* $f_{\text{NN}}(u;\theta)$ es la parte aprendida por la red neuronal.
* $\theta$ son los parámetros entrenables.

En el código, la dinámica aprendible se define como:

```julia
function lotka_volterra_ude!(du, u, p, t, p_true)
    x, y = u
    interaction, _ = nn([x, y], p, nn_st)
    du[1] =  p_true[1] * x + interaction[1]
    du[2] = -p_true[4] * y + interaction[2]
end
```

La red neuronal recibe el estado actual $(x,y)$ y devuelve dos términos de interacción. Por lo tanto, los parámetros que se optimizan no son directamente los coeficientes clásicos del sistema, sino los pesos internos de la red que modela la parte faltante de la dinámica.



---

## Predicción: evaluar un punto del espacio de parámetros

Para evaluar la función de costo en un punto $\theta$, el código debe resolver la UDE con esos parámetros. Esto aparece en la función `predict`:

```julia
function predict(ps)
    _prob = remake(prob_ude, u0=u0, tspan=tspan, p=ps)
    Array(solve(_prob, Vern7(), saveat=t_obs,
                abstol=1e-6, reltol=1e-6,
                sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))))
end
```

Conceptualmente, esta función realiza el siguiente procedimiento:

$$
\theta
\longrightarrow
\text{resolver ODE}
\longrightarrow
\hat{u}(t_1;\theta), \ldots, \hat{u}(t_N;\theta).
$$

Es decir, `predict(ps)` no hace una predicción directa como en una regresión lineal simple. Primero construye un problema de ODE usando los parámetros actuales de la red neuronal y luego lo resuelve numéricamente.

Esto hace que el problema de optimización sea más costoso: cada evaluación de la pérdida requiere integrar la dinámica.

La función `remake` permite reutilizar la estructura del problema `prob_ude`, pero reemplazando los parámetros actuales por `ps`. Esto es útil porque durante el entrenamiento los parámetros cambian en cada iteración.
---

## Función de costo

La función de costo está definida en el código como:

```julia
function loss(ps)
    ŷ = predict(ps)
    size(ŷ, 2) != length(t_obs) && return Inf
    return mean((ŷ .- target) .^ 2)
end
```

Esta función toma los parámetros actuales `ps`, resuelve la UDE y compara la trayectoria predicha con los datos de referencia `target`.

Matemáticamente, implementa un error cuadrático medio:

$$
L(\theta)
=
\frac{1}{N}
\sum_{i=1}^{N}
\left\|
\hat{u}(t_i;\theta) - u_i^{obs}
\right\|^2.
$$

Donde:

* $\hat{u}(t_i;\theta)$ es el estado predicho por la UDE en el tiempo $t_i$.
* $u_i^{obs}$ es el dato observado o de referencia.
* $\theta$ son los parámetros de la red neuronal.

Entonces el problema completo es:

$$
\theta^*
=
\arg\min_{\theta}
\frac{1}{N}
\sum_{i=1}^{N}
\left\|
\hat{u}(t_i;\theta) - u_i^{obs}
\right\|^2.
$$

Esta es la conexión principal con la clase: la función de costo mide el desacuerdo entre modelo y datos, y el optimizador busca parámetros que reduzcan ese desacuerdo.
## Entrenamiento en dos fases

El código entrena el modelo en dos etapas:

```julia
adam_iters = 1500
bfgs_iters = 1000
```

Primero usa **Adam** y luego **BFGS**. La idea no es cambiar la función de costo, sino cambiar el método que se usa para minimizarla.

El esquema general es:

$$
\text{inicialización aleatoria}
\longrightarrow
\text{Adam}
\longrightarrow
\text{BFGS}
\longrightarrow
\theta^*.
$$

Ambas fases minimizan la misma pérdida:

$$
L(\theta)
=
\frac{1}{N}
\sum_i
\left\|
\hat{u}(t_i;\theta) - u_i^{obs}
\right\|^2.
$$

Lo que cambia es la estrategia de actualización de los parámetros.

---

## Primera fase: Adam

La primera fase se define como:

```julia
optf = OptimizationFunction((ps, _) -> loss(ps), AutoZygote())

optprob1 = OptimizationProblem(optf, nn_ps)

res1 = Optimization.solve(optprob1, OptimizationOptimisers.Adam(),
                           callback = (ps, l) -> begin
                               push!(loss_history, l)
                               length(loss_history) % 50 == 0 &&
                                   println("  [Adam] iter $(length(loss_history))  —  Loss: $(round(l, digits=6))")
                               false
                           end,
                           maxiters = adam_iters)
```

Adam es un método local de primer orden. Esto significa que usa información del gradiente de la función de costo, pero no usa explícitamente el Hessiano.

En términos generales, Adam actualiza los parámetros usando una dirección basada en gradientes recientes y una escala adaptativa para cada coordenada. Por eso suele ser robusto en problemas con redes neuronales, donde la función de costo puede ser irregular, no convexa y de alta dimensión.


## Segunda fase: BFGS

Una vez terminada la etapa con Adam, el código usa la solución obtenida como punto inicial para BFGS:

```julia
optprob2 = OptimizationProblem(optf, res1.u)

res2 = Optimization.solve(optprob2,
                           OptimizationOptimJL.BFGS(linesearch=BackTracking()),
                           callback = (ps, l) -> begin
                               push!(loss_history, l)
                               n = length(loss_history) - adam_iters
                               n % 10 == 0 &&
                                   println("  [BFGS] iter $n  —  Loss: $(round(l, digits=6))")
                               false
                           end,
                           maxiters = bfgs_iters)
```

La parte clave es:

```julia
optprob2 = OptimizationProblem(optf, res1.u)
```

BFGS es un método quasi-Newton. En lugar de usar solamente el gradiente actual, intenta construir una aproximación de la curvatura de la función de costo a partir de la historia de pasos y gradientes.


En este ejemplo, Adam cumple el rol de acercar los parámetros a una región razonable del espacio de búsqueda. Luego BFGS refina la solución aprovechando información aproximada de curvatura.
## Resultado final del entrenamiento

Al finalizar BFGS, el código toma los parámetros finales:

```julia
nn_ps = res2.u
println("Entrenamiento finalizado. Loss final: $(round(loss_history[end], digits=6))")
```

Estos parámetros son la solución aproximada del problema de optimización. No se garantiza que sean el mínimo global, porque el problema es no convexo. Sin embargo, si la pérdida final es baja y las trayectorias ajustan bien, podemos interpretar que el entrenamiento encontró una solución útil.

Luego se resuelve nuevamente la UDE con los parámetros entrenados:

```julia
prob_trained = ODEProblem(nn_dynamics!, u0, tspan, nn_ps)
sol_trained  = solve(prob_trained, Vern7(), saveat=0.1, abstol=1e-6, reltol=1e-6)
```

Esto permite construir la trayectoria final aprendida por el modelo.

---

## Interpretación de la curva de pérdida

La primera salida importante para analizar la optimización es la curva de pérdida:

```julia
p2 = plot(1:length(loss_history), loss_history, label="Loss (MSE)",
          xlabel="Época", ylabel="MSE", title="Curva de entrenamiento",
          linewidth=2, color=:purple, yscale=:log10, legend=false)
```

Esta curva muestra el valor de la función de costo a lo largo de las iteraciones.


Conceptualmente, esperamos observar dos comportamientos:

1. Una primera caída importante durante Adam.
2. Una etapa de refinamiento durante BFGS.

La curva de pérdida es el gráfico más directamente relacionado con la clase de optimización: muestra cómo el algoritmo va reduciendo la función objetivo.
---

## Interpretación del gráfico de trayectorias

La figura principal se guarda como:

```julia
savefig(fig, "lotka_volterra_inverse_ude.png")
```

En ella se comparan:

* los datos de referencia,
* la trayectoria producida por el modelo entrenado,
* la curva de pérdida del entrenamiento.


Si la UDE entrenada se superpone bien con los puntos observados, significa que los parámetros encontrados por el optimizador generan una dinámica compatible con los datos.


## Interpretación de los términos de interacción aprendidos

El código también genera una figura para comparar las salidas de la red neuronal con los términos verdaderos de interacción:

```julia
true_f1 = @. -p_true[2] * X_traj * Y_traj
true_f2 = @.  p_true[3] * X_traj * Y_traj
```

y luego:

```julia
savefig(p3, "lotka_volterra_inverse_ude_interactions.png")
```

En este ejemplo sintético sabemos que la red debería aprender:

$$
NN_1(x,y) \approx -\beta xy,
$$

$$
NN_2(x,y) \approx \delta xy.
$$

Este gráfico permite analizar si la optimización recuperó algo parecido a la función faltante verdadera, no sólo si logró ajustar la trayectoria.


Desde el punto de vista de optimización, este gráfico ayuda a distinguir dos niveles de éxito:

* **Ajuste predictivo:** la trayectoria producida por el modelo se parece a los datos.
* **Recuperación de dinámica:** la red aprendió una función interna parecida a la verdadera.

El primer objetivo puede alcanzarse sin que el segundo sea perfecto, especialmente si hay pocos datos o si sólo observamos una parte limitada del espacio de estados.

---

## Interpretación de los heatmaps

La última figura se guarda como:

```julia
savefig(fig_nn, "lotka_volterra_inverse_ude_nn_heatmap.png")
```

Los heatmaps muestran los valores que devuelve la red neuronal en distintas regiones del espacio de estados $(x,y)$.

Desde el punto de vista de entrenamiento, esto permite inspeccionar cómo se comporta la función aprendida fuera de los puntos exactos usados durante la optimización.

Sin embargo, el código aplica una máscara:

```julia
sqrt(d) < mask_threshold ? nn([x, y], nn_ps, nn_st)[1][output_idx] : NaN
```

Esto hace que sólo se grafiquen valores cerca de la trayectoria observada. La razón es importante: la red fue entrenada principalmente con información en esa región del espacio de estados. Lejos de esa región, la salida de la red puede ser una extrapolación poco confiable.


---

## Qué nos enseña este ejemplo sobre optimización

Este ejemplo resume varias ideas centrales de la clase.

Primero, la función de costo no es simplemente una fórmula cerrada en los parámetros. Para evaluarla, hay que resolver una ODE. 

Segundo, los gradientes se calculan a través del solver de la ODE.

Tercero, se combinan dos métodos locales:

* Adam, como primera etapa robusta basada en gradientes.
* BFGS, como etapa de refinamiento quasi-Newton.

Cuarto, los resultados deben interpretarse en varios niveles:

* disminución de la pérdida,
* ajuste visual de trayectorias,
* recuperación de los términos de interacción,
* comportamiento de la función aprendida en el espacio de estados.

En resumen, el ejemplo muestra que ajustar una UDE no consiste sólo en resolver una ecuación diferencial, sino en resolver repetidamente muchas ecuaciones diferenciales dentro de un proceso de optimización.
