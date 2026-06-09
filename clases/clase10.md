---
title: No10 - PINNs Pt2
---

# Tips para entrenar PINNs

**Fecha:** 13/05/2026

:::{iframe} https://www.youtube.com/embed/5-SWMQoxbKs
:width: 100%
:::
---

# Visto en clase anterior

En la clase anterior se introdujeron las PINNs ({terms}`PINNs`) y se analizó la dualidad entre las restricciones suaves y fuertes estableciendo que entrenar una red de estas características es un problema de optimización de alta dificultad dado que consta de minimizar una función de costo que no sólo depende de los datos empíricos, sino también de que se satisfaga una ecuación diferencial sobre el dominio de la misma.


# Desafíos en el tratamiento de PINNs

Resulta conveniente mencionar y catalogar las posibles dificultades a las que se enfrentan las PINNs sea en la construcción de estas y/o en su entrenamiento.


## 1. Mal condicionamiento del problema de optimización

El problema de optimización sobre la función de costo puede ser altamente sensible a variaciones en los parámetros, por ejemplo, por poseer parámetros ajustados sobre escalas ampliamente distintas.

Por ello resulta de interés, encontrar algoritmos de búsqueda eficientes en función a la complejidad del problema y la capacidad disponible. Utilizando métodos de, al menos, primer orden y, preferentemente, cuasi-newton; véase [Clase No5](https://facusapienza.org/DM2026-Curso/clase5/#Métodos-de-búsqueda-local).


## 2. Balanceo de la función de costo

La función de costo total $\mathcal{L}_{TOT}$ resulta de una combinación lineal

$$\mathcal{L}_{TOT}(\theta) = \mathcal{L}_{EMP}(\theta) + \sum_i \lambda_i \mathcal{L}_{FIS}^{(i)}(\theta;x)$$

de la función de costo empírica $\mathcal{\mathcal{L}_{EMP}(\theta)}$ y las dadas por términos físicos $\mathcal{L}_{FIS}^{(i)}(\theta;x)$ (ecuación diferencial, condiciones iniciales, condiciones de borde, etcétera).
ads
Es fácil ver que el problema radica en **cómo elegir los hiperparámetros $\lambda_i$**.

* **Hiperparámetros fijos:** Establezco de antemano los hiperparámetros y los dejo constantes entre iteraciones. 

  * **L-Curve** Metodología empírica donde se grafica la pérdida empírica contra la pérdida física para diferentes valores de $\lambda$ y se toma por adecuado aquel que se encuentra en el punto máximo de curvatura como se ve en [fig. 1](#fig:seleccion_hiperparametro) a partir de escalar las magnitudes logarítmicamente. Se recomienda ver {cite}`Hansen_2001`.
    :::{figure} ./figures/no10_l_curve.svg
    :width: 100%
    :align: center
    :label: fig:seleccion_hiperparametro

    Selección de hiperparámetro $\lambda$ óptimo en función de la L-curve descrita por logaritmos de funciones de costo.
    :::

* **Hiperparámetros Adaptativos:** La mejor práctica en PINNs es recurrir a la meta-optimización de estos parámetros actualizando los $\lambda_i$ dinámicamente durante el entrenamiento buscando equiparar peso de las funciones de costo (o su efecto sobre los parámetros). 
  Algunos criterios para la adaptación incluyen:
  * Forzar a que todos los términos de la función de costo contribuyan de manera similar en magnitud.

    $$\mathcal{L}_{Emp} \simeq \lambda_{i} \mathcal{L}_{FIS}^{(i)}(\theta;x)$$

  * Asegurar que los gradientes de la pérdida física y empírica tengan magnitudes similares
   
    $$\lVert\nabla\mathcal{L}_{Emp}\rVert_2 \simeq \lambda_{i} \lVert\mathcal{L}_{FIS}^{(i)}(\theta;x)\rVert_2$$

    De esta manera, se evita que una componente domine la actualización del gradiente y, por ende, el siguente paso en el espacio de parámetros; de forma similar a la vista en  [fig. 2](#fig:desplazamiento_parametros) .

  :::{figure} ./figures/no10_desplazamiento_loss.svg
  :width: 90%
  :align: center
  :label: fig:desplazamiento_parametros

  Trayectoria de elección a través de igualar gradientes de las funciones de costo sobre el espacio de parámetros $\theta$.
  :::


## 3. Puntos de colocación 

Para imponer la restricción de la física, la ecuación diferencial se evalúa en un conjunto discreto de puntos en el dominio; de ahora en más denominado $\Omega$. Dado que evaluar la red (y calcular derivadas) tiene un alto costo computacional, la elección de estos puntos resulta crítica para el óptimo rendimiento y eficiencia del sistema propuesto.

Para asegurar esa elección existen distintas estrategias de *sampleo*:

1. **Grilla Uniforme (Latin Hypercube):** Distribución uniforme determinística de los puntos a evaluar sobre $\Omega$ como se ve en [fig. 3](#fig:latin_hypercube). 
    :::{figure} ./figures/no10_latin_hypercube.svg
    :width: 90%
    :align: center
    :label: fig:latin_hypercube

    Representación de la grilla de sampleo uniforme sobre dominio arbitrario $\Omega$.
    :::

    :::{caution} 
    Éste método corre riesgo de acoplarse erróneamente con las frecuencias características de la ecuación que se intenta resolver como se puede observar en [fig. 4](fig:sincronizacion)

    :::{figure} ./figures/no10_sincronizacion.svg
    :width: 90%
    :align: center
    :label: fig:sincronizacion

    Representación de inconsistencias por sampleo sobre grilla uniforme aplicado a función periódica. La vista en el dominio (derecha) ilustra conceptualmente la grilla uniforme subyacente que genera la distribución en el dominio físico $\Omega$ (izquierda).
    :::
    :::    
   
2. **Sampleo Uniforme Aleatorio:** Colección aleatoria y uniforme de puntos sobre $\Omega$ que pretende no poseer un claro patrón de sampleo como se puede identificar en [fig. 5](fig:sampleo_aleatorio).
    :::{figure} ./figures/no10_unif_random.svg
    :width: 90%
    :align: center
    :label: fig:sampleo_aleatorio

    Representación de sampleo uniformemente aleatorio sobre dominio arbitrario $\Omega$.La vista en el dominio (derecha) ilustra conceptualmente la grilla uniformemente aleatoria subyacente que genera la distribución en el dominio físico $\Omega$ (izquierda)
    :::
    
    En general, elude posibles sinterizaciones con frecuencias características de la ecuación analizada dada la aleatoriedad de la muestra.
    
    :::{caution} 
    Estos métodos de sampleo pueden generar *overfitting* sobre los puntos evaluados pues busca "simplemente" que la función de costo física se cumpla sobre los puntos específicos de la matriz de muestreo generando residuos nulos allí como se ve en [fig. 6](#fig:overfitting).
      :::{figure} ./figures/no10_overfitting.svg
      :width: 90%
      :align: center
      :label: fig:overfitting

      Representación de fenómeno de overfitting basado en una función de costo a partir del sampleo sobre una grilla fija en el dominio $\Omega$.
      :::
    :::
   
3. **Resampleo uniforme por época:** Sampleo uniforme de $\Omega$ de manera, o no aleatoria, en la iteración $k-ésima$ y, luego, se genera un nuevo muestreo  para la siguiente iteración $k+1$ como se puede apreciar en [fig.7](fig:resampleo). Buscando, de esta manera, no sobrepesar los puntos muestreados evitando comportamientos similares al visto en [fig. 6](#fig:overfitting).

    :::{figure} ./figures/no10_resampleo_epoch.svg
    :label: fig:resampleo
    :width: 90%
    :align: center

    Representación de fenómeno de overfitting visto en residuos de la función de costo.
    :::
4. **Importance Sampling (Sampleo por Importancia):** Estrategia más compleja donde no se samplea uniformemente, sino que se da mayor probabilidad de muestreo a las regiones de $\Omega$ donde la función de costo parece mayor buscando de esta manera que
   $$P(x_{i}^{k+1}) \propto e^{\alpha \mathcal{L}_{FIS}(x_{i}^{k})}$$
  donde $x_{i}^{k}$ es el la posición de sampleo $i$ en la iteración *$k$-ésima*.
  Concentrando, así, el poder de computo y el coste del mismo en las regiones de "mayor interés" o donde las pesa más tener un amplio muestreo.



## 4. Sesgo Espectral (Spectral Bias)

Las redes neuronales clásicas tienen la propiedad intrínseca de aprender **bajas frecuencias** mejor y más rápido que las altas frecuencias.

  :::{important} Comentario
    Por el comportamiento sobre las frecuencias, las redes neuronales clásicas, suelen ser tratadas como filtros pasa bajos dada esta facilidad en la implementación de frecuencias bajas. Se puede ver, por ejemplo, en [fig. 8](fig:sesgo_y_ff) donde parece ajustarse a la oscilación de menor frecuencia de la función buscada, un fenómeno formalmente teorizado por {cite}`rahaman2019spectral`.
  :::


Esto es muy problemático en ecuaciones diferenciales (como Navier-Stokes), donde las **altas frecuencias tienen un significado físico importante** y no son sólo ruido despreciable (por ejemplo, vórtices pequeños, turbulencias o singularidades). Si la red no puede capturar altas frecuencias, no podrá aproximar estas soluciones.


### **Soluciones al Sesgo Espectral:**

* **Fourier Features:** Se introduce una capa inicial no entrenable que pre-procesa las entradas espaciales/temporales utilizando funciones sinusoidales con desfases $\phi_k$ (o combinación lineal de senos y coseno) con diferentes frecuencias $\omega_k$ de manera que

  $$a^{(0)}(x) =  \begin{bmatrix}
    a^{(0)}_1(x)\\\
    a^{(0)}_2(x)\\
    \vdots 
    \end{bmatrix} =
    \begin{bmatrix}
    \sin(\omega_1 x + \phi_1)\\\
    \sin(\omega_2 x + \phi_2)\\
    \vdots 
    \end{bmatrix}$$
    
  siendo $a^{(0)}_{k}(x)$ la neurona $i-ésima$ de la capa "0", es decir, se mapea la variable de entrada a la red neuronal sobre una base de Fourier utilizando un diccionario de funciones, o en mejor de los casos, Transformada de Fourier (*FFT*).

  :::{figure} ./figures/no10_sesgo_espectral_y_ff.svg
  :label: fig:sesgo_y_ff
  :width: 90%
  :align: center

  Representación de predicciones de una red neuronal clásica con función de activación $\tanh(v)$ y otra haciendo uso *Fourier Features* en 200 iteraciones sobre una ecuación base con dos frecuencias características.
  :::

* **Múltiples Redes Concurrentes (Filosofía Bagging):** Otro método es aumentar los valores de entrada para mapear y capturar de mejor manera las frecuencias bajas que tenga la señal.
El desempeño de esta metodología se puede ver en [fig. 9](fig:bagging)
  :::{figure} ./figures/no10_bagging.svg
  :label: fig:bagging
  :width: 90%
  :align: center

  Representación de predicciones de una red neuronal aplicando *Bagging* en 200 iteraciones sobre una ecuación base con dos frecuencias características.
  :::


  * **Esquema de la Arquitectura**
    En este diseño, el input escalar $x$ se multiplica por factores enteros sucesivos en la primera capa modificada.
    Cada  una de estas entradas escaladas alimenta a una red neuronal independiente (arquitectura tipo *Ensemble/Bagging*), cuyas salidas parciales se combinan finalmente en una neurona de salida global $u$:

    $$x \longrightarrow \begin{bmatrix} 1x \\ 2x \\ 3x \\ \vdots \\ nx \end{bmatrix} \longrightarrow \begin{aligned} &\boxed{\text{Red}_1(\theta_1)} \longrightarrow u_1 \\ &\boxed{\text{Red}_2(\theta_2)} \longrightarrow u_2 \\ &\boxed{\text{Red}_3(\theta_3)} \longrightarrow u_3 \\ &\ \vdots \\ &\boxed{\text{Red}_n(\theta_n)} \longrightarrow u_n \end{aligned} \longrightarrow \boxed{\sum} \longrightarrow u_{\text{final}}$$

  * **Expresión Formal del Modelo**

    La salida final del sistema $u(x)$ se define como la combinación (frecuentemente un promedio o una suma ponderada) de $n$ redes neuronales "vainilla" independientes, donde cada red $i$ está parametrizada por sus propios pesos $\theta_i$:

    $$u_{\text{final}}(x) = \sum_{i=1}^{n} \text{Red}_i(i \cdot x ; \theta_i)$$

  * **Intuición del modelo de Bagging**

    Si la señal original contiene una frecuencia extremadamente baja o "chiquita" del orden de:

    $$f_{\text{baja}} = \frac{1}{n}$$

    Una red neuronal estándar (vainilla) tendría severas dificultades para aprenderla debido al sesgo hacia las altas frecuencias.
    Al forzar la multiplicación de $x$ por el factor $n$ en la $n$-ésima neurona de la primera capa, la frecuencia efectiva se transforma:

    $$f_{\text{efectiva}} = n \cdot \left(\frac{1}{n}\right) = 1$$

    Al convertirla en una frecuencia fundamental ($1$), la red asociada a ese bloque ya tiene la capacidad de aprender esa componente de la señal sin problemas.
    Finalmente, todas estas capacidades predictivas de diferentes escalas de frecuencia se combinan en la neurona de salida general $u$.

* **Multistage Networks (Filosofía Boosting):** 
  Se entrena una red neuronal que capturará principalmente las frecuencias bajas. Luego, se calcula el residuo (lo que no se pudo aprender, que son, por lo ya mencionado, frecuencias altas) y se entrena una *segunda* red neuronal para predecir exclusivamente ese residuo e iterar hasta abarcar la totalidad del espectro deseado y componiendo la predicción como la combinación lineal de estas redes, véase {cite}`Wang_Lai_Chiang_2023`. El ejemplo más común de este sistema es el algoritmo *LightGBM (LGBM)* que es utilizado para el gráfico presente en [fig. 10](lgbm).

  
  :::{figure} ./figures/no10_lgbm.svg
  :label: fig:lgbm
  :width: 90%
  :align: center

  Representación de predicciones de una red neuronal haciendo uso de *LightGBM* en 200 iteraciones sobre una ecuación base con dos frecuencias características.
  :::


   * **Esquema de la Arquitectura**

      A diferencia de la arquitectura en paralelo de la solución anterior, aquí las redes se encadenan secuencialmente en etapas (*multistage*).
      Cada red aprende la señal residual de la etapa previa:

      $$\boxed{\text{Red}^{(1)}} \longrightarrow \text{Residuo}^{(1)} \longrightarrow \boxed{\text{Red}^{(2)}} \longrightarrow \dots \longrightarrow \boxed{\sum} \longrightarrow \hat{y}$$

   * **Expresión Formal del Modelo**

      La predicción final del modelo compuesto por $K$ etapas es la suma acumulada de las salidas de cada una de las sub-redes entrenadas individualmente:

      $$\text{Red}_{\Theta}(x) = \sum_{i=1}^{K} \text{Red}^{(i)}_{\theta_i}(x)$$

    * **Mecanismo de Aprendizaje por Residuos**

      El entrenamiento se realiza de forma estrictamente secuencial bajo la siguiente lógica:

      1. **Etapa 1:** Se entrena la primera red $\text{Red}^{(1)}$ con los datos originales.
      Esta red suele capturar los patrones más fáciles y dominantes de la señal (por ejemplo, las **frecuencias bajas**).
      2. **Etapa 2:** Se calcula el primer residuo, es decir, lo que la primera red no pudo aprender:
        $$r^{(1)} = y - \text{Red}^{(1)}(x)$$
        Luego, la segunda red $\text{Red}^{(2)}$ se entrena utilizando $r^{(1)}$ como su objetivo (*target*).
        Esta red se ve forzada a capturar detalles más finos (como **frecuencias intermedias o altas**).
      3. **Etapa $i+1$:** De manera general, cada red subsiguiente aprende el residuo acumulado de todas las anteriores:
        $$r^{(i)} = y - \sum_{j=1}^{i} \text{Red}^{(j)}(x)$$

      Al finalizar el proceso, la combinación de todas las etapas permite reconstruir tanto la estructura global de los datos como sus detalles de alta frecuencia.

:::{bibliography}
:::