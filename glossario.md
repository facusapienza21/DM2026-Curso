---
title: Glosario
---

# Glosario

:::{glossary}
ODE
: Ecuación Diferencial Ordinaria (*Ordinary Differential Equation*). Ecuación que describe la evolución de un estado $u(t)$ en función de sus derivadas respecto al tiempo: $\frac{du}{dt} = f(u, t, \theta)$. — [Clase 2](clases/clase2.md)

PDE
: Ecuación en Derivadas Parciales (*Partial Differential Equation*). Generalización de las ODEs donde el estado depende de múltiples variables independientes (e.g., tiempo y espacio). — [Clase 2](clases/clase2.md)

NODE
: Ecuación Diferencial Ordinaria Neuronal (*Neural Ordinary Differential Equation*). Modelo donde la dinámica del estado está parametrizada por una red neuronal: $\frac{du}{dt} = f_\theta(u, t)$. Introducido por {cite}`chen2018neural`. — [Clase 2](clases/clase2.md)

UDE
: Ecuación Diferencial Universal (*Universal Differential Equation*). Generalización de las NODEs donde partes de una ecuación diferencial conocida son reemplazadas por redes neuronales. Introducido por {cite}`rackauckas2020universal`.

Parámetro
: Vector $\theta \in \mathbb{R}^p$ que caracteriza el comportamiento de un sistema dinámico. Inferir $\theta$ a partir de datos observados es uno de los problemas centrales del curso. — [Clase 2](clases/clase2.md)

Condición inicial
: Valor $u_0 = u(t_0)$ que especifica el estado del sistema en el tiempo inicial $t_0$. Junto con la ecuación diferencial ordinaria, determina unívocamente la solución (bajo condiciones de regularidad). — [Clase 2](clases/clase2.md)

Ajuste de trayectorias
: Método de inferencia estadística (*trajectory matching*) que estima los parámetros $\theta$ minimizando la discrepancia entre la solución numérica $u(t;\theta)$ y las observaciones $\{y_i\}$. Ver {cite}`ramsay2017dynamic`. — [Clase 2](clases/clase2.md)

Programación diferenciable
: Paradigma computacional (*differentiable programming*) que permite calcular gradientes de funciones definidas por programas, incluyendo simulaciones numéricas de ecuaciones diferenciales. Habilita el entrenamiento de modelos híbridos física-datos. Ver {cite}`Sapienza_2024` y {cite}`blondel2024elements`.

Sistema de Lotka-Volterra
: Modelo depredador-presa descripto por el sistema de ODEs: $\frac{dx}{dt} = \alpha x - \beta x y$, $\frac{dy}{dt} = \delta x y - \gamma y$. Es uno de los ejemplos recurrentes del curso. — [Clase 2](clases/clase2.md)

Inductive bias
: El sesgo inductivo (inductive bias) es el onjunto de supuestos, restricciones o conocimientos previos utilizados para condicionar el resultado de un algoritmo ante datos observados. — [Clase 7](clases/clase7.md)

:::