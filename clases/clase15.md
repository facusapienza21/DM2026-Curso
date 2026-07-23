---
title: No15 - Generalizaciones
---

# Generalizaciones

**Fecha:** 10/06/2026

:::{iframe} https://www.youtube.com/embed/joZAfB1-1Ns
:width: 100%
:::
**Autores:** Juan Ignacio Tollo ([@JuaniTollo](https://github.com/JuaniTollo))

Esta clase funciona como una **recapitulación y consolidación** de los temas del cuatrimestre: un pantallazo general de la metodología, más dos o tres ideas nuevas (gradient matching, SINDy, generalización a PDEs y SDEs) ahora que ya están las herramientas para entenderlas.

## 1. Metodologías científicas

La metodología científica puede clasificarse, de forma simplificada, en dos grandes grupos:

*   **Basados en Teorías Físicas:** ej. leyes físicas, ecuaciones diferenciales, leyes de conservación.
*   **Basados en Datos:** estadística y Machine Learning tradicional (deep learning, random forests, métodos estadísticos, AI).
    Son enfoques *data-driven*.
    La filosofía consiste en detectar patrones a partir de datos.
*   **Enfoque Híbrido:** el punto medio que integra modelos físicos y componentes aprendidas de los datos.
    Aquí se ubicaría *Physics Informed Machine Learning*.
    Se usan **ecuaciones diferenciales con componentes que se aprenden de los datos**.

La dicotomía no describe el funcionamiento de la ciencia en sí (las teorías físicas siempre usaron datos), sino la **metodología** que usamos para aprender algo de la naturaleza o construir conocimiento científico.

## 2. Clasificación de Modelado y Métodos

### Tipos de Modelado

*   **Modelado Directo (Forward):** conocido el modelo completo (ecuaciones y parámetros), calcular qué produce — ir de causas a efectos (ej. resolver Navier-Stokes para obtener el flujo de un fluido, o la ecuación de Schrödinger para obtener la función de onda).
*   **Modelado Inverso:** conocidas las observaciones, inferir lo que falta dentro de un modelo ya prescrito (parámetros, campos, funciones) — ir de efectos a causas.
    Hoy mucho de esto cae dentro de machine learning, pero es lo que tradicionalmente se conoce como *Data Assimilation*.
    Es el problema computacionalmente **difícil**: cada paso de la inferencia requiere resolver el forward, a menudo millones de veces ("correr una ecuación la puede correr cualquiera; el problema es correrla un millón de veces").

### Familias de Métodos Híbridos

1.  **Emuladores:** redes neuronales (típicamente convolucionales cuando los datos son tipo imagen) entrenadas para imitar un solver numérico de forma más **rápida** (no "más eficiente": es sobre velocidad, no eficiencia).
    Se entrenan generando muchos pares input–output con el solver original.
    *   *Ejemplo:* input = dominio, condiciones de borde, condición inicial, viscosidad → red convolucional → output = campo de velocidades del fluido.
    *   *Motivación:* uno no quiere resolver Schrödinger una vez, la quiere resolver un millón de veces (ej. probar un millón de moléculas para diseño de drogas).
        Es decir: el emulador acelera el forward porque en el fondo uno está pensando en el problema inverso.
    *   Metodológicamente usan poca física: son esencialmente machine learning entrenado con soluciones físicas.
2.  **Restricción Suave (PINNs):** una red neuronal $u_\theta$ juega el papel de la solución, y la ecuación entra como penalización en la función de costo (por eso "suave": se cumple solo aproximadamente).
    Pueden usarse en modo directo, pero un solver clásico tiende a ser mucho mejor para forward; el gran interés es el **modelado inverso**: se optimizan a la vez la red $u_\theta$ (la solución) y los parámetros desconocidos de la ecuación (un coeficiente, o un campo entero como la viscosidad $\nu(x)$).
    Los datos anclan a $u_\theta$ y el residuo de la ecuación obliga a que solución y parámetros sean consistentes — así se recupera lo que un solver no puede dar, porque el solver necesita conocer esos parámetros de antemano.
3.  **Restricción Fuerte (NODEs):** la ecuación diferencial está embebida en la arquitectura y se resuelve numéricamente.
    Lo interesante es el problema inverso: **aprender la física desde los datos**.

En resumen (el cuadro del pizarrón): las tres familias sirven para ambos tipos de problema, pero el interés principal siempre está en el inverso.

| Familia de Métodos | Modelado Directo / Forward | Modelado Inverso (*Data Assimilation* — **difícil**) |
|---|:---:|:---:|
| 1) Emuladores | ✓ | ✓ |
| 2) Restricción Suave (PINNs) | ✓ (Modo 1) | ✓ (Modo 2) |
| 3) Restricción Fuerte (NODEs, UDEs) |  | ✓ |

*Modo 1 / Modo 2:* los dos modos de una PINN — Modo 1 = usarla como solver de la ecuación (forward); Modo 2 = modelado inverso, aprender también los parámetros subyacentes (el modo interesante).

## 3. Problema general: optimización

En todos los casos, el problema se reduce a uno de optimización, que modelamos siempre como un problema de minimización.
Cualquier problema de estadística, machine learning, deep learning (y también física: mínima acción, mínima energía, mínimo potencial) puede pensarse como un problema de minimización.
Es en el fondo un **problema de minimización de una función de costo** $L(\theta)$.

*   **Frecuentistas:** buscan un estimador puntual $\theta^* = \arg\min_\theta L(\theta)$ (minimizar error / maximizar verosimilitud).
*   **Bayesianos:** buscan una distribución (*posterior*) sobre el espacio de parámetros, combinando *prior* (el análogo de la regularización) con la verosimilitud.
    Se samplean puntos aleatorios alrededor del mínimo del likelihood con un término de ruido que depende de la regularización.
    *   *Intuición térmica:* para ciertos priors, el posterior puede pensarse como un movimiento aleatorio sobre el mínimo de la función de costo con un término de ruido.
        Si se "enfría", converge al mínimo (caso frecuentista); si se "calienta", se dispersa y da una distribución.
*   Ambos enfoques tienen **más en común que diferencias**: en los dos hay que definir costo/verosimilitud y regularización/prior, en los dos el problema central es de optimización, y en los dos aparece el mismo problema de dimensión del espacio de parámetros.
    El bayesiano es más costoso porque devuelve una distribución (más información) en vez de un punto.

Para optimizar: métodos libres de gradiente sirven con $\dim(\theta)$ chica (~1–20, depende del problema); con dimensión grande se necesitan **gradientes**, y ahí aparece el problema técnico central del curso: cómo calcularlos (diferenciación automática, adjuntos).

*   Libre de gradiente
*   Gradiente
*   Mayor orden

## 4. Estrategias de Ajuste (Matching)

### Trajectory Matching

El método usado en todo el curso (PINNs, NODEs): minimizar la distancia entre la trayectoria del modelo y las observaciones:

$$\min_{\theta} \sum_{i} (y_i - u(t_i, \theta))^2$$

### Gradient Matching

Metodología alternativa en dos pasos, para la ecuación $\frac{du}{dt} = f(u, \theta, t)$ con observaciones ruidosas $y_i$ en tiempos $t_i \in [t_0, t_1]$:

1.  **Suavizado:** ajustar primero los puntos observados con un suavizador (ej. *cubic splines*), sin usar nada de la ecuación diferencial, para obtener $\hat{u}(t)$.
    Su derivada $\frac{d\hat{u}}{dt}$ sale analíticamente del suavizador.
2.  **Matching de gradientes:** con $\hat{u}$ ya fija, buscar $\theta$ tal que la derivada del suavizado coincida con lo que la dinámica predice (si la ecuación se satisface, el integrando es cero):

$$\min_{\theta} \int_{t_0}^{t_1} \left\| \frac{d\hat{u}}{dt}(t) - f(\hat{u}(t), \theta, t) \right\|^2 dt$$

En la práctica la integral se aproxima con una suma sobre los $t_i$.
A diferencia de trajectory matching, acá nunca se integra la ODE: se compara la derivada del suavizado contra $f$, no la trayectoria del modelo contra los datos.

**Ventajas:**

*   **No hay que resolver la ecuación diferencial** en cada paso de la optimización.
*   La función de costo es "más cuadrática" que en trajectory matching: si $f$ es **lineal en los coeficientes** $\theta$, el costo es cuadrático y el problema tiene solución (casi) analítica.
    En trajectory matching el costo sigue siendo no cuadrático aun con $f$ lineal.

**Desventajas:**

*   **Hay que poder "fitear" $u$:** el método implica tener datos de la trayectoria y poder suavizarlos por sí solos — esa es la gran limitación.
    En trajectory matching la ecuación diferencial genera la trayectoria y los datos solo la corrigen; acá el spline tiene que reconstruirla usando únicamente las observaciones, así que estas tienen que alcanzar para eso.
*   **Y no basta con fitear la trayectoria: hay que estimar bien su derivada**, que es mucho más sensible al ruido (aproximar bien $u$ no implica aproximar bien $u'$; un cociente incremental con ruido y $\Delta t$ chico amplifica el error).

**¿Cuándo conviene?** En **sistemas caóticos**: no tiene sentido fitear la trayectoria fuera del horizonte de Lyapunov, porque una diferencia infinitesimal en la condición inicial (incluso error de máquina) destruye la trayectoria a largo plazo.
Ahí matchear gradientes en vez de trayectorias es lo natural.
La idea es que el caos destruye la predictibilidad a largo plazo, no la validez local de la ecuación.

## 5. SINDy (Sparse Identification of Nonlinear Dynamics)

Método (familia gradient matching) para descubrir ecuaciones dinámicas simples a partir de datos.
Muy usado para sistemas caóticos (mínimo 3 ODEs para tener caos).

*   **Idea:** en los sistemas de ODEs típicos, los términos de $f$ son simples (lineales, productos $u_i u_j$, senos...) y la dinámica es **esparsa**: no aparecen todos los términos posibles, solo unos pocos.
    Se propone $f$ como combinación lineal de un **diccionario de $M$ funciones base** ($1,\ u_1,\ u_2,\ u_3,\ u_1 u_2,\ u_1 u_3,\ \ldots$), con coeficientes $\theta \in \mathbb{R}^{n \times M}$ ($n$ = número de ODEs) mayormente nulos:

$$\frac{du}{dt} = f(u, \theta) = \begin{bmatrix} \theta_{11} + \theta_{12}\,u_1 + \theta_{13}\,u_2 + \theta_{14}\,u_3 + \theta_{15}\,u_1 u_2 + \theta_{16}\,u_1 u_3 + \cdots \\ \theta_{21} + \theta_{22}\,u_1 + \theta_{23}\,u_2 + \cdots \\ \vdots \end{bmatrix}$$

    La fila $k$ de $\theta$ son los coeficientes de la ecuación $k$; todas las filas usan las mismas funciones base.
    La no-linealidad está toda en el diccionario: en $\theta$ el problema es **lineal**.

*   **Optimización con regularización $L_1$ (Lasso)** para forzar esparsidad:

$$\min_{\theta}\ \sum_i \left\| \frac{d\hat{u}}{dt}(t_i) - f(\hat{u}(t_i), \theta) \right\|^2 \;+\; \lambda \sum_{k,j} |\theta_{kj}|$$

    Cuanto más grande $\lambda$, más esparsidad → menos coeficientes activos → ecuaciones más simples.
    El primer término es cuadrático en $\theta$; el término $L_1$ no es diferenciable en cero pero se resuelve bien numéricamente.

*   **Pre-procesamiento — normalizar:** los términos del diccionario tienen **unidades diferentes** (adimensional, metros, metros cuadrados...), así que sus coeficientes no son comparables y penalizarlos juntos está mal.
    Hay que escalar cada función base (dividir por su desvío).
    Sutileza: los términos no son independientes ($u_1 u_2$ vs. $u_1$ y $u_2$), lo que genera además mucha colinealidad en el diccionario.

*   **Ventajas de que el problema sea Lasso clásico:**
    *   La **selección de $\lambda$** tiene teoría: los grados de libertad efectivos son el número de coeficientes no nulos, con lo cual se puede estimar el error de generalización y elegir $\lambda$ de manera fundada (a diferencia de una NODE regularizada, donde uno está más a ciegas).
    *   La **cuantificación de incertidumbre viene gratis**: al ser cuadrático en los coeficientes activos, se puede dar la incertidumbre de cada $\theta_{kj}$ analíticamente.
    *   Agrandar el diccionario no daña mucho: la esparsidad filtra los términos importantes, incluso con más features que observaciones.

*   **El gran desafío:** todo lo anterior asume el paso 1 (suavizado) resuelto, pero las derivadas $\frac{d\hat{u}}{dt}$ son muy sensibles al ruido de observación.
    Con datos sintéticos sin ruido las diferencias finitas funcionan; con ruido real lo amplifican y el método sufre.
    Si el suavizador es analíticamente diferenciable (ej. cubic splines), conviene derivar el suavizado en vez de los datos.

## 6. Generalizaciones a PDEs y SDEs

### Ecuaciones en Derivadas Parciales (PDEs)

Es el caso de interés en gran parte de la física y la geofísica (atmósfera, océanos, glaciares, manto terrestre): por lo general relacionan variables independientes como tiempo y espacio, y son más difíciles de resolver que las ODEs.
Los métodos usados pueden ser **diferencias finitas** y **elementos finitos**, y en general aparece el problema adicional de cómo mallar el dominio espacial.

Lo importante: **los métodos del curso se generalizan naturalmente** a PDEs.

*   **El método del adjunto se traslada directo.**
    Recordar el caso ODE lineal: si la dinámica es $\frac{du}{dt} = A u$, la variable adjunta satisface $\frac{d\lambda}{dt} = A^\top \lambda$ con condición **final** $\lambda(t_1) = \lambda_1$ — la misma dinámica pero con la matriz transpuesta, integrada hacia atrás.
    En PDEs el rol de "transponer la matriz" lo cumple tomar el **adjunto del operador diferencial**, que se obtiene por integración por partes.
    *Ejemplo (ecuación del calor):* si $\frac{\partial u}{\partial t} = \nabla \cdot (D \nabla u)$, el operador $\nabla \cdot (D \nabla\,\cdot)$ es **autoadjunto** (es su propia "transpuesta"), así que la ecuación adjunta es otra ecuación de calor: $\frac{\partial \lambda}{\partial t} = \nabla \cdot (D \nabla \lambda)$, con condición final $\lambda(t_1, x) = \lambda_1(x)$.

*   **Resolver hacia atrás = anti-difusión.**
    Como la ecuación adjunta tiene condición final, hay que integrarla hacia atrás en el tiempo.
    Haciendo el cambio de variable $t \to -t$ para escribirla como ecuación hacia adelante, el signo se invierte: $\frac{\partial \lambda}{\partial t} = -\nabla \cdot (D \nabla \lambda)$.
    Eso es una ecuación de **anti-difusión**: en vez de desparramar, concentra.
    Numéricamente es mucho más delicada que la difusión — el análogo en ODEs de integrar una exponencial creciente en vez de una decreciente.

:::{note}
Ojo con la convención de signos: en otras referencias la ecuación adjunta aparece como $\frac{\partial \lambda}{\partial t} = -\nabla\cdot(D\nabla\lambda)$ con condición final, que es la misma ecuación del pizarrón tras el cambio $t \to -t$.
:::

*   **Por qué el adjunto es crítico en PDEs:** los parámetros suelen ser **campos** ($D(x)$, la condición inicial $u_0(x)$, forzantes), o sea millones de parámetros tras discretizar la malla.
    Calcular el gradiente perturbando costaría una resolución de la PDE *por parámetro*; con el adjunto cuesta **dos resoluciones en total**, independiente del número de parámetros.
    Es la diferencia entre imposible y factible (así funciona la asimilación de datos en meteorología, e.g. 4D-Var).

*   **Costo de memoria:** almacenar la trayectoria de un campo con $x, y, z, t$ es mucho más pesado que la de una ODE; checkpointing y adjuntos se vuelven imprescindibles en alta dimensión.

### Ecuaciones Diferenciales Estocásticas (SDEs)

Se agrega ruido a la dinámica.
Partiendo de la discretización de una ODE ($du = f(u)\,dt$), se suma un término de ruido:

$$du = f(u, t)\,dt + g(u, t)\,\sqrt{dt}\; \xi, \qquad \xi \sim \mathcal{N}(0,1)$$

que en notación continua se escribe con el **proceso de Wiener**, cuyo incremento distribuye $dW_t \sim \mathcal{N}(0, dt)$:

$$du = f(u, t)\,dt + g(u, t)\,dW_t$$

*   El escaleo $\sqrt{dt}$ es el único que funciona: con potencias mayores o menores de $dt$ la solución degenera o explota.
*   $f$ es el **drift** (hacia dónde va la partícula sin ruido) y $g$ la **difusión** (el empuje aleatorio constante).
    La trayectoria ya no es única: cada realización del ruido da una distinta, y no tiene por qué parecerse a la solución suave.
*   *Usos:* sistemas intrínsecamente aleatorios (el ejemplo canónico: precios de mercado; también variables físicas como temperatura o presión), y también **aproximar sistemas deterministas** complejos, especialmente al modelar parametrizaciones (una parte del sistema en vez de la física completa).

**Estimación de parámetros.** En el diagrama del curso había un *modelo de estado* ($\theta \to u(t)$, determinista) y un *modelo observacional* ($y_i = u(t_i) + \varepsilon_i$, aleatorio).
Con SDEs el modelo de estado también es aleatorio (con ruido independiente del observacional).
El principio de **máxima verosimilitud no cambia**, pero ahora hay que promediar sobre las trayectorias posibles:

$$\mathcal{L}(\theta) = \mathbb{E}_{u} \left[ \prod_i \exp\left(-\frac{(y_i - u(t_i))^2}{2\sigma_i^2}\right) \right]$$

El orden importa: primero el producto sobre observaciones dada una trayectoria, después la esperanza sobre trayectorias (invertirlo asumiría observaciones desacopladas de la historia, lo cual es incorrecto).

*   **Problema:** samplear trayectorias a fuerza bruta no funciona — la probabilidad de que una trayectoria aleatoria pase cerca de *todos* los puntos observados es ínfima (como sacar cara muchas veces seguidas): casi todas contribuyen $\approx 0$ a la esperanza.
*   **Particle Filtering** (basado en *importance sampling / resampling*): se samplean trayectorias, se descartan las que se alejan de los datos, se conservan las cercanas y se re-samplean nuevas partículas a partir de ellas.
    Así se estima la verosimilitud de forma inteligente y eficiente.
    (Filosóficamente emparentado con los algoritmos genéticos: matar las partículas que no funcionan y reproducir las que sí.)
*   También existen métodos para calcular **gradientes** de esta cantidad (una esperanza estimada por Monte Carlo) respecto de $\theta$.
