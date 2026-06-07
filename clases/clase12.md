---
title: No12 - Programación Diferencial Pt2
---

# Programación Diferenciable: Métodos Forward Pt2

**Fecha:** 01/06/2026

:::{iframe} https://www.youtube.com/embed/hlUUiCIZ_mE
:width: 100%
:::

Estas notas prosiguen la temática de {term}`Programación diferenciable` de la clase pasada, retomando la introducción a {term}`Diferenciación automática` (*Automatic Differentation (AD)*).
En particular, comenzamos con una implementación de *AD* en Julia, haciendo uso de los {term}`números duales`.

## Números duales

Se define a los números duales como una extensión de los reales, comenzando por definir a $\epsilon$ número abstracto cumpliendo

$$
\epsilon^{2} = 0\quad ;\quad \epsilon \neq 0 
$$

y escribiendo a todo número dual como 

$$
x_{\epsilon} = x_1 + \epsilon x_2
$$
en donde a $x_1$ será el *valor* de $x_{\epsilon}$ y $x_2$ será su *derivada*.

En lo que a la *AD* respecta, los números duales son una manera muy fácil y directa de implementarla en un lenguaje de programación con herramientas de POO como lo es Julia.

```julia
@ksdef struct DualNumber{F <: AbstractFloat}
            value: F
            derivative: F
end
```
Queremos ser capaces de operar sobre los números duales, sumar, multiplicar, aplicar funciones elementales.
Una de las cosas buenas de Julia es su simplicidad a la hora de extender operaciones a nuevas estructuras de datos.

```julia
#Define operations on dual numbers
function Base.:(+)(a::DualNumber, b::DualNumber)
    res_value = a.value + b.value
    res_derivative = a.derivative + b.derivative
    return DualNumber(res_value, res_derivative)
end

function Base.:(*)(a::DualNumber, b::DualNumber)
    res_value = a.value * b.value
    res_derivative = a.value * b.derivative + a.derivative * b.value 
    return DualNumber(res_value, res_derivative)
end
```
Esto nos va a permitir instanciar los números duales y operar sobre ellos, trasladando siempre en la *parte dual* la derivada correspondiente a la ejecución de las operaciones.

Si creamos 2 números duales:

```julia
a = DualNumber(1.0, 0.0)
b = DualNumber(2.0, 1.0)
```
Y hacemos la operación *a + b* deberia devolver:

```julia
a + b = DualNumber(3.0, 1.0)
```
Analogamente, si creamos 2 números duales:

```julia
a = DualNumber(0.9, 0.0)
b = DualNumber(1.4, 1.0)
```
Y hacemos la operación *a * b* deberia devolver:

```julia
a * b = DualNumber(1.4 * 0.9, 0.9)
```
De esta forma, se puede observar como la parte dual arrastra el valor de la derivada, y esto se puede hacer cuantas veces uno quiera.

Sigamos con uan función un poco más compleja.
```julia
#Define operations on dual numbers
function Base.:(sin)(a::DualNumber)
    res_value = sin(a.value)
    res_derivative = cos(a.value) * a.derivative
    return DualNumber(res_value, res_derivative)
end
```
Como estas funciones, se pueden crear tantas como operaciones tengamos, sin importar que sean unitarias, binarias, etc.

Siempre lo que uno consigue es que la primer componente tenga el *valor* y la segunda componente sea su *derivada*.

Los números duales son muy útiles a la hora de calcular derivadas parciales e incluso direccionales,
debido a que esta estructura permite flexibilizar hacia donde esta derivando uno.

*Ejemplo*<br>
Si hubiesemos querido derivar respecto a b, solo deberiamos haber modificado el input de la siguiente manera:
```julia
a = DualNumber(0.9, 1.0)
b = DualNumber(1.4, 0.0)
a * b = DualNumber(1.4 * 0.9, 1.4)
```
En la práctica uno no crea todas estas funciones desde 0, ya que existe una libreria que contiene todas estas funciones y muchas más.
Uno solo la importa y usa todas las herramientas que provee esta librería.

*Ejemplo*
```julia
using ForwardDiff

x = ForwardDiff.Dual(2.0, 1.0)

y = x^2 + 3x
```

A diferencia de lo que hacíamos en diferencias finitas el valor de la derivada usando *AD* es exacto.

Si probamos esto con diferencias finitas:

La derivada se aproxima mediante

$$
\frac{f(a+\epsilon)-f(a)}{\epsilon}
$$

```julia
f(x) = sin(x * 0.9)

epsilon = 1e-10
@show (f(a.value + ϵ ) - f(a.value)) / ϵ  
```
Mientras que con *AD* el cálculo de la derivada es exacto, con diferencias finitas comienzan a haber errores de truncación debido a
la sensibilidad del resultado con respecto a $\epsilon$.

Una contra de este método es que viene con un costo de memoria más alto, debido a que ahora estamos trabajando no solo con su valor
sino que también con su derivada.

En ecuaciones diferenciales, uno propaga el número dual en el solver númerico y consigue la solución y la derivada de esa solución
con respecto a los parámetros.

Veamos una representación de lo que sucede con cada método:

![Gráfico](../images/clase12.png)

Se observa que con la *solución exacta de diferencias finitas*, el error baja y luego vuelve a subir debido al error de truncado.

Además, la otra curva refleja la *solución exacta de diferenciacion compleja*, donde se ve que la misma baja hasta $10^{-16}$, el
error de máquina.

Por último, la curva violeta representa el error de *forward AD*. En este caso, se puede observar que el mismo se adapta totalmente
a la tolerancia ya que en ambos gráficos la curva se mantiene constante sobre la tolerancia en cada caso respectivamente.

En resumen, *diferencias fínitas* es el método menos exacto ya que contiene error de truncado, mientras que *diferenciación compleja* baja hasta error de máquina a partir de un cierto $\epsilon$. *Forward AD* no depende de $\epsilon$, por lo que, en caso de que la tolerancia fuese el error de máquina, la curva se mantendría constante en ese valor.

## Comparación matemática: Diferenciación Compleja vs Forward AD

Para entender la diferencia fundamental entre ambos métodos forward, desarrollamos paso a paso la derivada de $f(x) = \sin(x^2)$ en un punto $x$ cualquiera.

El resultado esperado es: $f'(x) = 2x\cos(x^2)$

### Método 1: Diferenciación Compleja

**Idea central**: Evaluamos la función en un punto complejo $x + i\varepsilon$ y extraemos la derivada de la parte imaginaria.

- **Parte real** del resultado → nos da el **valor** de la función $f(x)$
- **Parte imaginaria** del resultado → contiene la información de la **derivada** $f'(x)$


**Paso 1 — Definimos la Variable compleja**

Definimos:

$x = x_1 + ix_2$ con $i^2 = -1$

y para calcular la derivada, tomamos $x_2 = \varepsilon$ muy pequeño:

$x = x_1 + i\varepsilon$

Aca:
- $x_1$ es el punto donde evaluamos (parte real)
- $\varepsilon$ es el "paso" de perturbación (parte imaginaria)

**Paso 2 — Elevamos al cuadrado la variable compleja**

$$(x_1 + i\varepsilon)^2 = x_1^2 + 2i\varepsilon x_1 + (i\varepsilon)^2 = x_1^2 + 2i\varepsilon x_1 - \varepsilon^2$$

En particular,

- **Parte real**: $a = x_1^2 - \varepsilon^2$ ← contiene el valor real
- **Parte imaginaria**: $b = 2x_1\varepsilon$ ← contiene información de la derivada

**Paso 3 — Aplicamos el seno**

Para este caso, usamos la siguiente identidad 

$$\sin(a + bi) = \sin(a)\cosh(b) + i\cos(a)\sinh(b)$$

$$\sin(x^2) = \sin(x_1^2 - \varepsilon^2)\cosh(2x_1\varepsilon) + i\cos(x_1^2 - \varepsilon^2)\sinh(2x_1\varepsilon)$$

**Paso 4 — Extraemos la derivada**

La fórmula de diferenciación compleja nos dice que la derivada está en la parte imaginaria, dividida por $\varepsilon$:

$$\frac{df}{dx} = \lim_{\varepsilon \to 0} \frac{\text{Im}(f(x + i\varepsilon))}{\varepsilon} = \lim_{\varepsilon \to 0} \frac{\cos(x_1^2 - \varepsilon^2)\sinh(2x_1\varepsilon)}{\varepsilon}$$

**Paso 5 — Tomamos el límite $\varepsilon \to 0$**

Como $\varepsilon^2 = 0$

$$\frac{df}{dx} = \lim_{\varepsilon \to 0} \frac{\cos(x_1^2 - \varepsilon^2) \cdot \sinh(2x_1\varepsilon)}{\varepsilon} = \lim_{\varepsilon \to 0} \frac{\cos(x_1^2) \cdot \sinh(2x_1\varepsilon)}\varepsilon$$

Luego, usando que $\sinh(z) \approx z$ para $z \to 0$ entonces $\sinh(2x_1\varepsilon) = 2x_1\varepsilon$:

$$\frac{df}{dx} = \lim_{\varepsilon \to 0} \frac{\cos(x_1^2) \cdot 2x_1\varepsilon}{\varepsilon}$$

y dividiendo por $\varepsilon$:

$$\frac{df}{dx} = \cos(x_1^2) \cdot 2x_1 $$

### Método 2: Forward AD (Números Duales)

**Idea central**: Evaluamos la función en un "número dual" $x + \varepsilon$ y la derivada aparece directamente como coeficiente de $\varepsilon$.

- Parte sin $\varepsilon$ (término "real") → nos da el valor de la función $f(x)$
- Coeficiente de $\varepsilon$ (término "dual") → nos da directamente la derivada $f'(x)$

**Paso 1 — Definimos la Variable dual**

Definimos $x_\varepsilon = x_1 + \varepsilon$ con $\varepsilon^2 = 0, \quad \varepsilon \neq 0$

Aca:
- $x_1$ es el punto donde evaluamos (parte "valor real")
- El coeficiente de $\varepsilon$ será la derivada (parte "dual")

**Paso 2 — Elevamos al cuadrado**

$$x_\varepsilon^2 = (x_1 + \varepsilon)^2 = x_1^2 + 2x_1\varepsilon + \varepsilon^2 = x_1^2 + 2x_1\varepsilon$$

Obtenemos un número dual donde:
- **Parte valor** (sin $\varepsilon$): $x_1^2$ 
- **Parte dual** (coeficiente de $\varepsilon$): $2x_1$ 

**Paso 3 — Aplicamos Seno**

Expandimos por Taylor:

$$\sin(x_1^2 + 2x_1\varepsilon) = \sin(x_1^2) + \cos(x_1^2) \cdot (2x_1\varepsilon) - \frac{\sin(x_1^2)}{2}(2x_1\varepsilon)^2 + ...$$

Pero $(2x_1\varepsilon)^2 = 4x_1^2\varepsilon^2 = 0$, por lo que todos los términos de orden $\geq 2$ desaparecen:

$$\sin(x_\varepsilon^2) = \sin(x_1^2) + \varepsilon \cdot 2x_1\cos(x_1^2)$$

Entonces volviendo a la idea central, ya tenemos de forma explicita la derivada (tomamos la parte del coeficiente de $\varepsilon$):

$$\frac{df}{dx} = 2x_1\cos(x_1^2)$$


**Conclusión:** En el método de números duales (Forward AD) podemos obtener la derivada de forma inmediata (solamente una en una evaluación de la función tenemos el valor real y su derivada) sin necesidad de usar trucos matematicos, meternos con los limites o tener que extraer la parte imaginaria de un resultado. Es decir, si integramos esto dentro de un solver de ODEs, en cada paso tenemos la derivada correspondiente y de forma directa y automática.