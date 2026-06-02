---
title: No1 - Introducción
---

## Introducción

**Fecha:** 08/04/2026

**Material:** [Slides de la clase](../slides/2026-04-08.pdf)

:::{iframe} https://www.youtube.com/embed/LDp4yZyDHdQ
:width: 100%
:::


En esta clase se introducen los temas que se van a estudiar durante la cursada.
Trataremos de entender con ecuaciones físicas y herramientas de aprendizaje automático cuál es la dinámica de distintos sistemas físicos.
Las aplicaciones de la materia van a estar motivadas por la geofísica.

## Motivación

Comprender y predecir el comportamiento de sistemas físicos complejos es uno de los grandes desafíos de la computación científica moderna.
Como hilo conductor de la clase nos enfocaremos en el estudio de **glaciares** como ejemplo de sistema dinámico. 

:::{note} ¿Qué es un glaciar?

Un glaciar es un cuerpo persistente de hielo que evoluciona dinámicamente a lo largo del tiempo.

Aunque pueda parecer una estructura estática, un glaciar fluye lentamente bajo el efecto de la gravedad y responde a cambios climáticos y ambientales.
:::

### Importancia de los glaciares

* **Reguladores del clima a nivel global**
    - Influencia sobre el nivel del mar.
    - Reflexión de radiación solar.
    - Relación con el cambio climático.

* **Recursos hídricos**
    - Reserva de agua dulce.
    - Regulación del suministro de agua.

* **Componente cultural**
    - Impacto social y económico.
    - Importancia turística y ecológica

### Dinámica de un glaciar

Los glaciares son sistemas dinámicos.
Aunque estén formados por hielo, fluyen lentamente de manera similar a un río.
Su movimiento ocurre principalmente debido a la gravedad.

Existen distintos mecanismos que explican este movimiento.

- **Deformación interna:**  
  El hielo puede modelarse como un medio viscoso que se deforma internamente.
  Las capas superiores suelen moverse más rápido que las inferiores.

- **Deslizamiento:**  
  El glaciar interactúa con la roca subyacente mediante fuerzas de fricción.
  Dependiendo de estas interacciones, el hielo puede quedar adherido o deslizarse sobre la superficie.

Dado que los glaciares son sistemas dinámicos, existen leyes físicas que describen su evolución.

### ¿Cómo se modela un glaciar?

El modelado de glaciares combina física, matemática y computación.
Las ecuaciones que describen estos sistemas suelen ser complejas y, en general, no poseen soluciones analíticas exactas.
Por esta razón, es necesario utilizar métodos numéricos y simulaciones computacionales para aproximar sus soluciones.


Para calibrar los parámetros de las ecuaciones que describen la dinámica de un glaciar se utilizan distintas fuentes de datos.

1. **Registros históricos:**  
   por ejemplo mapas, dibujos y pinturas históricas.

2. **Observaciones in-situ:**  
   mediciones directas realizadas sobre el glaciar en su entorno natural.

3. **Observaciones satelitales:**  
   han revolucionado las Ciencias de la Tierra en las últimas décadas.
   Permiten obtener grandes volúmenes de datos observacionales.

4. **Paleogeografía:**  
   permite reconstruir el comportamiento histórico de glaciares a partir de evidencia geológica.

## La revolución de la IA

En los últimos años, distintos avances tecnológicos impulsaron el uso de técnicas de Inteligencia Artificial en Ciencias de la Tierra.

- **Volumen de datos:**  
  todo algoritmo de Inteligencia Artificial requiere grandes cantidades de datos.
  Existen alrededor de 280.000 glaciares en el mundo.
  Además, las observaciones satelitales generan grandes volúmenes de información.

- **Capacidad de cómputo:**  
  actualmente se dispone de hardware con gran capacidad de procesamiento.

- **Algoritmos y software:**  
  el desarrollo de herramientas y bibliotecas de código abierto facilitó la aplicación de modelos de IA.

Todos estos factores permiten aplicar técnicas de Inteligencia Artificial al modelado de glaciares.
Estos avances permiten combinar modelos físicos tradicionales con herramientas modernas de aprendizaje automático.

## ¿Física o Datos?

Pondremos en práctica lo que se conoce como {term}`Modelado híbrido`.
El modelado híbrido combina modelos físicos tradicionales con herramientas modernas de aprendizaje automático.

```{figure} figures/ModeladoHibrido.png
---
width: 90%
name: fig-hybrid-modeling
---
Relación entre modelos físicos, modelos basados en datos y modelado híbrido.
```

:::{important} Más allá de los glaciares

Las ideas presentadas en esta clase no son exclusivas del modelado de glaciares.

La combinación de modelos físicos, métodos numéricos, datos observacionales y herramientas de aprendizaje automático aparece en numerosas disciplinas científicas e ingenieriles.

Algunos ejemplos incluyen:
- predicción climática,
- dinámica de fluidos,
- oceanografía,
- y sistemas complejos en general.

El modelado híbrido es actualmente una de las estrategias más utilizadas en computación científica moderna.
:::

## Objetivo de la materia

* Desarrollar y entender metodología que nos permita asimilar datos (heterogéneos, irregulares) a nuestros modelos físicos.
    1. Entender mejor la física del sistema que estamos tratando de modelar.
    2. Mejorar capacidad predictiva de los modelos.

* Desarrollar un marco estadístico para definir dichos problemas.

* Desarrollar software que nos permita realizar esta tarea de manera eficiente y accesible.

* Tener un marco conceptual para pensar en problemas donde leyes físicas y aprendizaje automático hablen el mismo lenguaje.

:::{important} Software Abierto

Vamos a hacer énfasis en software abierto

1. Ciencia abierta y colaborativa

2. Reproducibilidad y replicabilidad de resultados

3. Accesibilidad a recursos y metodología sin limitación de capacidad de computo
:::
---
## Introducción a Physics-Informed Machine Learning (PIML)
El {term}`PIML <Physics Informed Machine Learning (PIML)>` es una disciplina donde convergen tres conceptos principales:
* **Modelado físico**: ecuaciones diferenciales que modelan la física del problema.
* **Aprendizaje automático**: algoritmos de *machine learning*.
* **Simulaciones numéricas**: para resolver las ecuaciones diferenciales ya mencionadas.
```{figure} figures/venn-diagram-piml.png
---
width: 60%
name: fig-piml-diagram
---
Caracterización del PIML a partir de la intersección entre el modelado físico, aprendizaje automático y simulaciones numéricas.
```
## Machine Learning: ¿es una caja negra de algoritmos?
1. Tiene poca interpretabilidad, no sabemos cómo aprende.
2. No extrapola bien fuera de los datos entrenados.
3. Es difícil de implementar.

Bajo esas críticas se fundamenta la idea de apoyarnos en leyes físicas a la hora de hacer uso del aprendizaje automático.
:::{note} Ejemplo: oscilador armónico amortiguado
Se pueden ver diferencias a la hora de entrenar una red neuronal **sin** restricciones físicas ({numref}`fig-nn`) y **con** restricciones físicas ({numref}`fig-pinn`).
```{figure} figures/nn.png
---
width: 60%
name: fig-nn
---
Aproximación de la solución del problema del oscilador armónico amortiguado utilizando aprendizaje automático entrenado en los puntos naranjas de la curva.
```
```{figure} figures/pinn.png
---
width: 60%
name: fig-pinn
---
Aproximación de la solución del problema del oscilador armónico amortiguado usando PIML entrenado en los puntos naranjas de la curva.
```
:::
---

## Diferentes métodos dentro del SciML

Dentro del **Scientific Machine Learning (SciML)** se pueden distinguir tres formas principales de combinar datos, aprendizaje automático y leyes físicas.

```{figure} figures/fig-sciml-methods.png
---
width: 85%
name: fig-sciml-methods
---
Tres categorías principales dentro del SciML: métodos puramente basados en datos, restricciones físicas suaves y restricciones físicas fuertes.
```

La diferencia central entre estas estrategias está en **cómo aparece la física dentro del modelo**.

En un enfoque **data-driven**, el modelo aprende directamente desde datos. La física no aparece de manera explícita, sino que queda implícita en los ejemplos utilizados para entrenar.

En un enfoque con **restricciones suaves**, la física aparece como un término adicional en la función de pérdida. El modelo no está obligado a satisfacer exactamente la ecuación diferencial, pero es penalizado cuando se aleja de ella.

En un enfoque con **restricciones fuertes**, la física forma parte de la estructura misma del modelo. La red neuronal se integra dentro de una ecuación diferencial, un solver numérico o una simulación diferenciable.


### Problemas directos e inversos

En el modelado físico aparecen dos tipos de problemas fundamentales: los **problemas directos** y los **problemas inversos**.

```{figure} figures/forward-inverse.png
---
width: 75%
name: fig-forward-inverse
---
Comparación entre un problema directo y un problema inverso.
```

Un **problema directo** consiste en predecir el estado de un sistema o su evolución temporal a partir de condiciones iniciales, parámetros y leyes físicas conocidas.

Por ejemplo, si conocemos el espesor inicial de un glaciar, su geometría, ciertas condiciones climáticas y una ley de movimiento, podemos intentar predecir cómo evolucionará en el tiempo.

Un **problema inverso**, en cambio, parte de observaciones y busca inferir parámetros, estados ocultos o leyes del sistema.

:::{note} Ejemplo: Aplicado a glaciares
En glaciología, un problema inverso típico sería estimar parámetros relacionados con la fricción basal, el deslizamiento del hielo o la geometría de la roca subyacente a partir de observaciones satelitales de velocidad, espesor o elevación.
:::


### 1. Aprendizaje automático de sistemas físicos

El primer enfoque consiste en usar aprendizaje automático para aproximar directamente la relación entre entradas y salidas de un sistema físico.

```{figure} figures/data-driven-physical-systems.png
---
width: 80%
name: fig-data-driven-physical-systems
---
Ejemplo de aprendizaje automático aplicado a sistemas físicos.
```

La idea principal es reemplazar, aproximar o acelerar una simulación costosa mediante un modelo entrenado con datos.

En lugar de resolver una ecuación diferencial compleja cada vez que se necesita una predicción, se entrena una red neuronal para aprender una aproximación de la solución.

$$
\text{condiciones iniciales} \longmapsto \text{estado futuro del sistema}.
$$

Esto permite obtener predicciones mucho más rápidas, aunque con una limitación importante: el modelo depende fuertemente de los datos con los que fue entrenado.

:::{warning}
El riesgo principal de este enfoque es que el modelo puede fallar al extrapolar fuera del régimen observado durante el entrenamiento.
:::


### 2. Leyes físicas impuestas de manera suave

El segundo enfoque consiste en incorporar la física como una penalización dentro de la función objetivo.

```{figure} figures/soft-constraints-pinn.png
---
width: 80%
name: fig-soft-constraints-pinn
---
Ejemplo de incorporación de leyes físicas mediante restricciones suaves.
```

En este caso, el modelo se entrena minimizando una pérdida de la forma

$$
\mathcal{L}(\theta)
=
\mathcal{L}_{\text{datos}}(\theta)
+
\lambda \mathcal{L}_{\text{física}}(\theta).
$$

El primer término mide qué tan bien el modelo ajusta los datos observados. El segundo término mide qué tan bien respeta la ecuación diferencial que describe el sistema.

Por ejemplo, si una ecuación física tiene la forma

$$
F(u) = 0,
$$

entonces podemos penalizar al modelo cuando su predicción no cumple esa ecuación:

$$
\mathcal{L}_{\text{física}}(\theta)
=
\|F(u_\theta)\|^2.
$$

Este es el principio detrás de las **Physics-Informed Neural Networks (PINNs)**.

La física no se impone exactamente, sino que se introduce como una restricción suave: el modelo puede violarla, pero paga un costo por hacerlo.

:::{important}
La elección de $\lambda$ es importante: si es demasiado chico, el modelo ignora la física; si es demasiado grande, puede ajustar mal los datos.
:::



### 3. Leyes físicas impuestas de manera fuerte

El tercer enfoque consiste en imponer la física directamente dentro de la estructura del modelo.

```{figure} figures/hard-constraints-ude.png
---
width: 80%
name: fig-hard-constraints-ude
---
Ejemplo de modelo híbrido donde una red neuronal se integra dentro de un solver físico.
```

En este caso, no se agrega simplemente una penalización física a la función de pérdida. En cambio, el modelo se construye de forma tal que la evolución del sistema debe pasar por una ecuación diferencial o un solver numérico.

Una forma general de escribir este enfoque es

$$
\frac{du}{dt} = f(u,t,\theta),
$$

donde una parte de $f$ puede ser conocida por la física y otra parte puede ser aprendida mediante una red neuronal.

Por ejemplo,

$$
\frac{du}{dt}
=
f_{\text{física}}(u,t)
+
f_{\text{NN}}(u,t;\theta).
$$

La red neuronal aprende la parte desconocida del modelo, pero la evolución completa sigue estando gobernada por una ecuación diferencial.

Estos modelos suelen llamarse **Universal Differential Equations (UDEs)**.


---
---

## Ecuaciones diferenciales universales (UDE)

Las **ecuaciones diferenciales universales** introducen una red neuronal dentro de una ecuación diferencial.


La idea central es utilizar la mayor cantidad posible de conocimiento físico existente y usar modelos basados en datos sólo para las partes de la ecuación que necesitan ser aprendidas o extendidas.

En este enfoque, la física se respeta porque la ecuación diferencial se resuelve numéricamente. La red neuronal no reemplaza por completo al modelo físico, sino que se introduce dentro de la dinámica del sistema.

Una forma general de escribir este tipo de modelos es

$$
\frac{du}{dt} = f(u,t,\theta),
$$

donde una parte de $f$ puede estar definida por leyes físicas conocidas y otra parte puede ser aprendida a partir de datos.

Por ejemplo,

$$
\frac{du}{dt}
=
f_{\text{física}}(u,t)
+
f_{\text{NN}}(u,t;\theta).
$$

Rackauckas et al., *Universal Differential Equations for Scientific Machine Learning* (2020).

---

## Al final del día, todo es un problema de optimización

Todo aquello que no conocemos o no sabemos describir exactamente se puede **parametrizar** de alguna manera.

Una vez parametrizado el problema, el objetivo es encontrar los valores de los parámetros que mejor ajusten los datos, respeten las restricciones físicas o reproduzcan el comportamiento observado del sistema.

Esto lleva naturalmente a algoritmos de búsqueda de óptimos.

Salvo en casos de dimensión muy chica, no alcanza con probar valores manualmente. Es necesario calcular gradientes para saber cómo modificar los parámetros y mejorar la función objetivo.

---

## Programación diferencial

Los modelos inversos están basados en **programación diferencial**.

La programación diferencial es un paradigma de programación orientado a calcular gradientes o sensibilidades de programas de computadora.

En este contexto, entrenar o calibrar modelos con muchos parámetros requiere diferenciar a través del modelo computacional.

En particular, para modelos físicos y modelos híbridos, esto puede implicar diferenciar a través de un solver numérico.

La idea general es:

$$
\text{modelo} \longrightarrow \text{predicción} \longrightarrow \text{función objetivo} \longrightarrow \text{gradientes}.
$$

Estos gradientes permiten ajustar los parámetros del modelo mediante algoritmos de optimización.



### Métodos de programación diferencial

Existen distintos métodos para calcular derivadas y sensibilidades en modelos computacionales.

```{figure} figures/differential-programming-methods.png
---
width: 70%
name: fig-differential-programming-methods
---
Métodos de programación diferencial según si son discretos o continuos, y si son forward o reverse.
```

Los métodos pueden clasificarse según dos ejes principales.

Por un lado, pueden ser **discretos** o **continuos**. Los métodos discretos diferencian el programa o algoritmo ya discretizado. Los métodos continuos trabajan sobre la formulación continua del problema, por ejemplo sobre una ecuación diferencial.

Por otro lado, pueden ser **forward** o **reverse**. Los métodos forward propagan sensibilidades hacia adelante, mientras que los métodos reverse propagan información hacia atrás, como ocurre en backpropagation.

Algunos ejemplos son:

- diferencias finitas,
- diferenciación simbólica,
- diferenciación automática forward,
- diferenciación automática reverse,
- métodos adjuntos discretos,
- métodos adjuntos continuos,
- ecuaciones de sensibilidad forward.



### Aplicaciones

La programación diferencial y el modelado híbrido aparecen en distintas áreas de las ciencias de la Tierra.

Entre las aplicaciones principales se encuentran:

- glaciología,
- paleomagnetismo y tectónica de placas,
- geomagnetismo.

También pueden aparecer aplicaciones en otros sistemas físicos donde sea necesario combinar modelos, datos y calibración de parámetros.

---

## Qué vamos a tener que aprender

Para trabajar con Physics-Informed Machine Learning es necesario combinar herramientas de varias áreas.

Por el lado del modelado físico, necesitaremos trabajar con **ecuaciones diferenciales**.

Por el lado computacional, necesitaremos herramientas de **análisis numérico**, ya que muchas ecuaciones no pueden resolverse de forma analítica y requieren métodos numéricos.

También necesitaremos conceptos de **aprendizaje automático**, especialmente métodos de regresión como redes neuronales y procesos Gaussianos.

Además, será importante incorporar herramientas de **estadística**, tanto frecuentista como bayesiana, incluyendo el modelado del ruido.

Finalmente, aparecerán temas vinculados con **programación diferencial**, teoría de control y asimilación de datos.

