---
title: No4 - Métodos numéricos
---

# Metodos Numericos

**Fecha:** 20/04/2026

---

## Video de la clase

:::{iframe} https://www.youtube.com/embed/Rdvz8KRA1JQ 
:width: 100%
:::

Hasta ahora hablamos de NODEs y ODEs. Sabemos que hay una ecuación diferencial, pero no discutimos cómo se calculan sus soluciones en la práctica. A lo largo de esta clase abordaremos cómo se resuelven numéricamente las trayectorias $\boldsymbol{u}(t)$.

Las ideas que veremos se extienden naturalmente a ecuaciones más complejas, como ecuaciones en derivadas parciales, ecuaciones estocásticas, etc., donde los métodos numéricos son más sofisticados.


---

## Solvers numéricos para ODEs

Estamos pensando que queremos resolver,

$$
\dfrac{d\boldsymbol{u}}{dt} = \boldsymbol{f}_\theta(\boldsymbol{u},t)
$$

donde $\boldsymbol{f}_\theta : \mathbb{R}^N \times \R \to \mathbb{R}^N$ es una función parametrizada por un conjunto de parametros $\theta$ fijo, y $\boldsymbol{u}(t) \in \mathbb{R}^N$ es la trayectoria buscada.


Existen dos grandes familias de métodos que engloban al resto:
- [Multi-step](#metodos-multi-step)
- [Runge Kutta](#metodos-runge-kutta)

:::{note}
En cualquiera de los casos, tendremos algún tipo de discretización del eje temporal. De modo que la trayectoria $\boldsymbol{u}(t)$ no estará evaluada de manera continua a cada tiempo $t$, sino que será evaluada en una grilla temporal $t_0, t_1, \ldots, t_m$. Esto da lugar a $\boldsymbol{u}_0, \boldsymbol{u}_1, \ldots, \boldsymbol{u}_m$, donde $\boldsymbol{u}(t_0) = \boldsymbol{u}_0$ es la condición inicial y $\boldsymbol{u}(t_m) = \boldsymbol{u}_m$ es la solución dada por el solver en los puntos de la grilla.

El paso temporal $\Delta t$ se puede fijar o determinar dinámicamente por un algoritmo, pero los métodos descritos a continuación no dependen de una forma específica de discretización temporal (pueden aplicarse tanto a pasos fijos como variables).
:::


### Métodos Multi-Step

Los métodos Multi-Step se posicionan en un tiempo $t^m$ y calculan la solución del paso siguiente basándose en una historia de pasos anteriores, siguiendo una relación de recurrencia general:

$$
\sum_{i=1}^{d_1} \alpha_{i} \boldsymbol{u}^{m+i} = \sum_{j=1}^{d_2} \beta_{j} \boldsymbol{f}(\boldsymbol{u}^{m-j}, t^{m-j})
$$


con $d_1, d_2 \in \mathbb{N}$ hiperparámetros del solver que determinarán la precisión del método con respecto a la resolución temporal $\Delta t$, es decir, indexan el orden del solver. El lado izquierdo de la ecuación representa la evolución hacia el futuro, mientras que el lado derecho utiliza la pendiente calculada en puntos del pasado.

:::{note} Ejemplo: Adams-Bashforth

Un caso muy simple y conocido de este tipo de métodos es el descrito por Adams-Bashforth de segundo orden:

$$
\boldsymbol{u}^{m+1} = \boldsymbol{u}^m + \dfrac{\Delta t}{2} \left(3 \boldsymbol{f}(\boldsymbol{u}^m, t^m) - \boldsymbol{f}(\boldsymbol{u}^{m-1}, t^{m-1}) \right)
$$

Este algoritmo determina el siguiente paso $\boldsymbol{u}^{m+1}$ utilizando la información de los dos pasos anteriores. Las ventajas de este método son su bajo costo computacional, ya que solo requiere una nueva evaluación de la función $\boldsymbol{f}$ por cada paso temporal (reutilizando las evaluaciones ya calculadas de pasos previos).
:::

### Métodos Runge-Kutta

Los métodos de Runge-Kutta se basan en aproximar la solución evaluando la función $\boldsymbol{f}$ en distintos puntos intermedios dentro del paso temporal.

$$
\boldsymbol{u}^{m+1} = \boldsymbol{u}^m + \sum_{i=1}^{s} b_i \boldsymbol{k}_i 
$$
donde $b_i$ son coeficientes, $s \in \mathbb{N}$ son las *etapas*, y
$$
\boldsymbol{k}_i = \boldsymbol{f}\left(\boldsymbol{u}^m + \sum_{j=1}^{s} a_{ij} \boldsymbol{k}_j, t^m + c_i \Delta t \right)
$$
con $a_{ij}$ y $c_i$ hiperparámetros del solver. La idea de este algoritmo es no evaluar $\boldsymbol{f}$ únicamente en $\boldsymbol{u}^m$, sino aproximar mejor la dinámica evaluándola en puntos intermedios dentro del intervalo temporal $[t^m, t^{m+1}]$. Esto se logra construyendo una combinación de pendientes $\boldsymbol{k}_i$ que capturan mejor la curvatura de la trayectoria.


:::{note} Ejemplo: Método del punto medio

Especificando para $s=2$ se recupera el conocido método de punto medio:

$$
\boldsymbol{k}_1 = \boldsymbol{f}(\boldsymbol{u}^m, t^m)
$$

$$
\boldsymbol{k}_2 = \boldsymbol{f}\left(\boldsymbol{u}^m + \dfrac{\Delta t}{2} \boldsymbol{k}_1, t^m + \dfrac{\Delta t}{2} \right)
$$

$$
\boldsymbol{u}^{m+1} = \boldsymbol{u}^m + \Delta t \boldsymbol{k}_2
$$
:::

:::{important} Runge-Kutta de cuarto orden (RK4)
El método de Runge-Kutta de cuarto orden (RK4) es uno de los más populares debido a su equilibrio entre precisión y costo computacional. Para $s=4$, los coeficientes específicos son:
$$
\begin{aligned}
\boldsymbol{k}_1 &= \boldsymbol{f}(\boldsymbol{u}^m, t^m) \\
\boldsymbol{k}_2 &= \boldsymbol{f}\left(\boldsymbol{u}^m + \dfrac{\Delta t}{2} \boldsymbol{k}_1, t^m + \dfrac{\Delta t}{2} \right) \\
\boldsymbol{k}_3 &= \boldsymbol{f}\left(\boldsymbol{u}^m + \dfrac{\Delta t}{2} \bold{k}_2, t^m + \dfrac{\Delta t}{2} \right) \\
\boldsymbol{k}_4 &= \boldsymbol{f}\left(\boldsymbol{u}^m + \Delta t \boldsymbol{k}_3, t^m + \Delta t \right) \\
\boldsymbol{u}^{m+1} &= \boldsymbol{u}^m + \dfrac{\Delta t}{6} (\boldsymbol{k}_1 + 2\boldsymbol{k}_2 + 2\boldsymbol{k}_3 + \boldsymbol{k}_4)
\end{aligned}
$$
Este método es de orden 4, lo que significa que el error local por paso es proporcional a $\Delta t^5$, y el error global es proporcional a $\Delta t^4$. RK4 es ampliamente utilizado en la práctica dada su convergencia en $O(\Delta t^4)$.
:::

### Euler explícito

El método numérico más conocido para la resolución de ODEs es el Euler explícito. Este método es el de orden 1 (el más simple dentro de los métodos consistentes) y puede reobtenerse como un método de [Runge-Kutta](#metodos-runge-kutta) de primer orden ($s=1$), o como un método [Multi-Step](#metodos-multi-step) con $d_1=1$ y $d_2=0$:

$$
\boldsymbol{u}^{m+1} = \boldsymbol{u}^m + \boldsymbol{f}(\boldsymbol{u}^m, t^m)\Delta t
$$

> **Obs 1:** Los puntos donde el solver devuelve la solución $\boldsymbol{u}^m$ no necesariamente coinciden con los tiempos donde el solver evalúa la función. Una opción es devolver un interpolador (solución densa) y evaluarlo en cualquier tiempo deseado.

> **Obs 2:** Estos métodos son conocidos como *forward*: dada la ecuación, devuelven la trayectoria $\boldsymbol{u}^m$. Los métodos inversos (inverse problems) buscan, dada la solución, determinar los parámetros $\theta$. Esto es parte de lo que nos interesará resolver más adelante.

## Segunda Motivación de las NODEs

En {cite}`chen2018neural`, se presenta la idea de usar una red neuronal para parametrizar la función $f$ de una ODE dada su similitud con el método de [Euler explícito](#euler-explicito). Dada una red neuronal con capas ocultas $h_i \in \mathbb{R}^{n_i}$, se puede definir la acción de una capa cualquiera a partir de la anterior como 

$$
h_{i+1} = h_i + g(h_i)
$$

donde $g(h_i)$ es una función no lineal, por ejemplo $g(h_i) = \sigma(W_i h_i + b_i)$, con parámetros $\theta_i = (W_i, b_i)$. Este tipo de arquitectura se conoce como una red neuronal **residual** (ResNet).

Este esquema suele ser mejor para modelar dinámicas o series temporales, ya que introduce una noción de continuidad entre capas:

$$
h_{i+1} = h_i + \tilde{g}(h_i)\Delta t_i
$$

De aquí se observa una discretizacion de un sistema dinamico (Ver: [Euler explícito](#euler-explicito)), lo que motiva la idea de "pasar al continuo" y definir una NODE como la solución de la siguiente ODE:

$$
\dfrac{dh}{dt} = \tilde{g}(h(t))
$$

Dado que desconocemos quién es $\tilde{g}$, entonces debemos aprenderla a partir de la "trayectoria", es decir, resolver esta ecuación utilizando un aproximador universal:

$$
\dfrac{dh}{dt} = \mathrm{NN}_\theta(h(t))
$$

## Pre-Entrenamiento

A veces resulta necesario pre-entrenar la red neuronal que reemplaza la función $f$ de la ODE, para aproximar los parámetros de la red a una región del espacio donde la solución física tenga sentido. Esto funciona como una "inicialización": le damos a la red neuronal una idea previa de cómo debería comportarse la dinámica.

Si tenemos una idea de la física subyacente dada por $\tilde{F}_\theta$, entonces podemos usar esta información para guiar el entrenamiento de la red neuronal optimizando primero la siguiente función de pérdida:

$$
\mathcal{L}_{pre}(\theta) = \sum_{i=1}^{N} \|\tilde{F}_\theta(X_i) - \mathrm{NN}_\theta(X_i)\|^2
$$

Esto es mucho menos costoso que optimizar la función de pérdida que se obtiene al comparar la solución de la ODE con datos observados, ya que no requiere resolver la ecuación diferencial en cada iteración del entrenamiento.

Observemos que en esta etapa no utilizamos datos experimentales, sino datos sintéticos generados a partir del modelo físico $\tilde{F}_\theta$. Es decir, usamos nuestro conocimiento previo del sistema para "enseñarle" a la red una primera aproximación de la dinámica.


## Ecuaciones Diferenciales Universales (UDE) — Sistema de Lotka-Volterra - Implementación en Julia


> El siguiente código [03_LV_forward_UDE](https://github.com/facusapienza21/DM2026-Curso/tree/main/code/03_LV_forward_UDE) presenta un ejemplo concreto de pre-entrenamiento aplicado al modelo de Lotka-Volterra.

Una **Ecuación Diferencial Universal** o **UDE** (*Universal Differential Equation*) es una técnica que combina modelos físicos mecánicos con redes neuronales. La idea central es conservar las leyes físicas conocidas y usar una red neuronal para aprender los componentes del sistema que son desconocidos o difíciles de modelar.

En este ejemplo se aplica una UDE al sistema de **Lotka-Volterra** (presa-depredador). El sistema verdadero es:
$$
\frac{dx}{dt} = \alpha x - \beta x y \qquad \frac{dy}{dt} = \delta x y - \gamma y
$$

El enfoque UDE aquí consiste en suponer que conocemos las tasas de nacimiento/muerte lineal ($\alpha x$ y $-\gamma y$), pero desconocemos cómo interactúan las especies. Por lo tanto, reemplazamos los términos de interacción por una **Red Neuronal** $\mathrm{NN}(x, y)$:

$$
\frac{dx}{dt} = \alpha x + \mathrm{NN}(x,y)_1 \qquad \frac{dy}{dt} = -\gamma y + \mathrm{NN}(x,y)_2
$$
El objetivo es que la red neuronal aprenda los términos originales.

Utilizaremos la librería **Lux.jl**, que se destaca en el ecosistema de Julia por su enfoque funcional, separando la estructura de los datos (parámetros).

```{code-cell} julia
using Lux       
```
En esta sección definimos la arquitectura de la red neuronal que actuará como el componente "universal" de nuestra ecuación diferencial. 

```{code-cell} julia
begin
    # Definición de la arquitectura de la red 
    nn = Chain(
        Dense(2 => n_hidden, tanh),
        Dense(n_hidden => n_hidden, tanh),
        Dense(n_hidden => 2)
    )

    # Inicialización de parámetros y estado de forma reproducible
    rng_nn       = MersenneTwister(seed)
    nn_ps, nn_st = Lux.setup(rng_nn, nn)
end;
```

Esto define una red neuronal dada por los siguientes puntos:

1. **Arquitectura de la Red (`Chain`):**
   Se define una red neuronal densa secuencial. La entrada son 2 valores (población de presas $x$ y depredadores $y$). La red tiene capas ocultas cuyo tamaño está definido por la variable `n_hidden` y una capa de salida de 2 valores, que representan los términos de interacción del sistema físico.

2. **Funciones de Activación (`tanh`):**
   Se utiliza la tangente hiperbólica en lugar de funciones como ReLU. Es fundamental usar funciones de activación "suaves" (continuamente diferenciables) para que el *solver* no encuentre problemas numéricos al calcular las trayectorias y sus gradientes.

3. **Funcionamiento de Lux:** 
   *   **`nn`**: La estructura lógica de la red (inmutable).
   *   **`nn_ps` (Parameters)**: Los pesos y sesgos que se van a entrenar. 
   *   **`nn_st` (State)**: Variables que la red necesita pero que no son entrenables.

   Este esquema es lo que permite "inyectar" la red neuronal dentro de una ecuación diferencial y optimizarla como si fuera un parámetro físico más.

4. **Reproducibilidad:**
   Mediante `MersenneTwister(seed)`, nos aseguramos de que los pesos iniciales de la red sean siempre los mismos en cada ejecución del *notebook*, lo cual es vital para el proceso de depuración y para compartir resultados científicos.

Antes de integrar la red neuronal dentro de la ecuación diferencial, realizamos la etapa de [Pre-Entrenamiento](#pre-entrenamiento). El objetivo es que la red aprenda algo sobre los términos de interacción antes de enfrentarse al problema de optimización completo. Esto puede ser útil, por ejemplo, si se tiene un sistema de presas-depredadores que no sigue estrictamente la dinámica de Lotka-Volterra y queremos que una red neuronal aprenda este término de interacción desconocido.

Para no enfrentarnos directamente al problema de optimización completo, el cual es altamente complicado, una opción es decirle a la red neuronal lo que sabemos: "la dinámica se parece a Lotka-Volterra", y entrenarla para situar los parámetros en una región adecuada del espacio de parámetros.

Este paso es opcional (controlado por la variable `do_pretrain`), pero mejora significativamente la convergencia del entrenamiento posterior. Si uno intenta resolver el problema sin entrenar la red, el solver puede no "explotar", pero el resultado probablemente no tenga sentido físico (la red neuronal está "descalibrada").

```{code-cell} julia
trained_ps, loss_history = if do_pretrain
    let nn_ps = nn_ps
        X_mat = reduce(hcat, [[x, y] for x in range(0.0, 60.0, length=30)
                                      for y in range(0.0, 60.0, length=30)])
        Y_mat = reduce(hcat, [[-β * x * y, δ * x * y]
                               for x in range(0.0, 60.0, length=30)
                               for y in range(0.0, 60.0, length=30)])

        pretrain_loss(ps) = mean((nn(X_mat, ps, nn_st)[1] .- Y_mat) .^ 2)
        opt_state = Optimisers.setup(Adam(1e-3), nn_ps)
        history = Float64[]

        for epoch in 1:pretrain_epochs
            grads = Zygote.gradient(pretrain_loss, nn_ps)[1]
            opt_state, nn_ps = Optimisers.update(opt_state, nn_ps, grads)
            push!(history, pretrain_loss(nn_ps))
        end
        nn_ps, history
    end
else
    nn_ps, Float64[]
end;
```
El propósito de este bloque es generar un dataset sintético, dado por las variables:

- `X_mat`: contiene los pares $(x, y)$ como columnas, que son las entradas de la red.
- `Y_mat`: contiene los valores "correctos" de los términos de interacción, calculados analíticamente a partir del modelo de Lotka-Volterra.

A partir de estos datos sintéticos, se define una función de costo (`pretrain_loss`). Se utiliza el Error Cuadrático Medio (MSE) entre la salida de la red y los valores analíticos:

$$
\mathcal{L}_{pre}(\theta) = \frac{1}{N} \sum_{i=1}^{N} \left\| \mathrm{NN}_\theta(x_i, y_i) - \begin{pmatrix} -\beta x_i y_i \\ \delta x_i y_i \end{pmatrix} \right\|^2
$$

Luego, comienza la etapa de entrenamiento. Se usa el optimizador Adam con una tasa de aprendizaje de $10^{-3}$. `Optimisers.setup` inicializa el estado interno del optimizador (momentos de primer y segundo orden) asociado a los parámetros `nn_ps`. En cada época:

- `Zygote.gradient` calcula el gradiente de la función de pérdida respecto a `nn_ps` mediante diferenciación automática.
- `Optimisers.update` aplica el paso de Adam, devolviendo el nuevo estado del optimizador y los parámetros actualizados.
- Se registra el valor de la pérdida en `history` para visualizar la convergencia.

:::{note}  
Este bloque `let` crea una copia local mutable de `nn_ps`. En Pluto, las variables son reactivas y globalmente inmutables; el bloque `let` permite modificar `nn_ps` localmente dentro del loop sin afectar el estado global del notebook.
:::


:::{bibliography}
:::