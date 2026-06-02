# Información

<!-- :::{table} -->
<!-- :align: center -->
| | |
| --- | --- |
| Fecha final para asignar proyecto | 06/05/2026 |
| Fecha final de validación de proyectos | 13/05/2026 |
| Fecha de exposición de posters | 17/06/2026 12:00-14hrs |
| Fecha de entrega | 19/06/2026 |
| | |
<!-- ::: -->

## Contenido del proyecto

El proyecto final consiste en la aplicación de los contenidos del curso a un problema o tema de interés del estudiante.
El proyecto puede combinar teoría, metodología y aplicación en la proporción que cada estudiante o grupo considere más adecuada según sus intereses.
Se recomienda discutir la elección del tema durante las [horas de oficina](../general/informacion.md) o en el horario de clase.
Cada proyecto debe incluir:

- **Informe escrito** de entre 5 y 10 páginas, describiendo el problema, los datos, el método utilizado y los resultados obtenidos. Todas las fuentes utilizadas (artículos, libros, código, datasets) deben estar correctamente citadas en el informe.
- **Código** (opcional) reproducible que respalde los resultados presentados en el informe. En caso de incluir código, el link a un repositorio debe estar incluido en el informe. Cualquier lenguaje de programación es bienvenido, aunque el instructor recomienda aprovechar y darle una oportunidad a Julia.

:::{note} Proyectos grupales
¡Trabajar en grupo es altamente recomendado! Colaborar permite abordar problemas más ambiciosos y aprender de otros estudiantes, tal como ocurre en la práctica científica real. Cuando el proyecto es realizado en grupo, se espera una mayor cantidad y profundidad de material.
:::

:::{note} Proyectos teóricos
Los proyectos de naturaleza teórica son más que bienvenidos. En ese caso, se espera un informe más extenso y un análisis en mayor profundidad del tema elegido, demostrando un entendimiento sólido de los conceptos involucrados.
:::

## Cómo presentar el proyecto

La presentación de proyectos por parte de los estudiantes se realizará mediante un Issue en el [repositorio de la materia](https://github.com/facusapienza21/DM2026-Curso).
Para ver cómo abrir un issue, consultar la sección [Reportar un problema](../contribucion.md#reportar-un-problema-issue) en la página de contribución.

:::{important} Plantilla para la propuesta
La propuesta debe abrirse usando la plantilla predefinida haciendo click [aquí](https://github.com/facusapienza21/DM2026-Curso/issues/new?template=proyecto.md).
Por favor mantener la estructura general de la misma (titulo, integrantes), pero sentirse libre de agregar cuanta información sea necesaria en la descripción del proyecto.
Para ver un ejemplo de cómo debe verse una propuesta, pueden consultar este [issue de ejemplo](https://github.com/facusapienza21/DM2026-Curso/issues/9).
:::

El issue también puede usarse como espacio de conversación para discutir y refinar el proyecto antes de su aprobación.
El contenido del issue puede ser modificado por los autores en cualquier momento, por lo que no es necesario completarlo todo de una sola vez.

El proceso de creación del proyecto tiene dos fases:

**Fase 1 — Aprobación del estudiante/grupo** (`Proyecto: asignado`) — fecha límite: **06/05/2026**
Una vez abierto el issue, el instructor verificará que los integrantes del grupo estén correctamente registrados.
Al aprobar el grupo, el issue recibirá la etiqueta `Proyecto: asignado`.
Esto **no** implica que el contenido del proyecto haya sido aprobado.

**Fase 2 — Aprobación del contenido** (`Proyecto: validado`) — fecha límite: **13/05/2026**
El instructor revisará la descripción del proyecto, el dataset y el método propuesto.
En caso de ser necesario, se pedirán cambios o aclaraciones a través de comentarios en el issue.
Una vez que el contenido sea satisfactorio, el issue recibirá la etiqueta `Proyecto: validado`, lo que indica que el proyecto está oficialmente aprobado y los estudiantes pueden comenzar a trabajar.

## Entrega

El proyecto final debe ser enviado a `fsapienza@dm.uba.ar` antes de la medianoche de la fecha límite de entrega.
El informe debe tener entre 7 y 12 páginas (sin contar referencias ni apéndices) y describir el proyecto de forma completa: motivación, metodología, resultados y conclusiones.
Escribir el reporte como si fuera un artículo científico dirigido a un lector con conocimientos generales del área pero no necesariamente familiarizado con el problema específico.
El informe tiene que ser lo más autocontenido posible: incluir definiciones, notación y contexto suficientes para que el problema y la solución sean comprensibles sin necesidad de consultar fuentes externas.
Incluir todas las referencias académicas relevantes al final del documento (las mismas no cuentan como páginas del reporte).
Cuanto mayor sea la cantidad de miembros de un grupo/proyecto, más producción se espera en dicho reporte.

El informe debe seguir la siguiente estructura general:

1. **Introducción** — motivación del problema, contexto y objetivos.
2. **Métodos** — descripción del modelo, ecuaciones, datos y técnicas utilizadas.
3. **Resultados** — presentación clara de los resultados obtenidos, con figuras y tablas cuando corresponda.
4. **Discusión y conclusiones** — interpretación de los resultados y posibles trabajos futuros.
5. **Agradecimientos** *(Acknowledgments)* — reconocimiento de colaboraciones, financiamiento y uso de herramientas de IA.
6. **Datos y Software** — información sobre datos y software utilizados.
7. **Referencias** — todas las fuentes citadas en el texto.

:::{important} Datos y Software
Incluir una sección breve al final del informe (antes de las referencias) indicando cómo acceder a los datos y al código utilizados en el proyecto.

- **Datos:** indicar la fuente del dataset; si es de acceso público, incluir el link o DOI y la licencia de uso si corresponde.
- **Software:** mencionar las librerías y versiones utilizadas. Si el código desarrollado está disponible en un repositorio de GitHub público (altamente recomendado), incluir el link.


**Ejemplo de Data Availability Statement:**
> *Los datos utilizados en este trabajo corresponden al dataset público XYZ, disponible en [enlace/DOI]. El código desarrollado para este proyecto está disponible en [github.com/usuario/proyecto](https://github.com) bajo licencia MIT. Este trabajo utilizó Julia v1.10 con los paquetes `DifferentialEquations.jl` (REF) y `Lux.jl`.*
:::

:::{important} Agradecimientos y use de AI
Los Agradecimientos son el lugar apropiado para declarar el uso de herramientas de inteligencia artificial.
Pueden usar LLMs u otras herramientas de IA para realizar el trabajo, pero deben aclarar qué herramientas usaron, para qué y cómo.
La responsabilidad sobre el contenido del informe es siempre de los autores: cualquier afirmación incorrecta o referencia inexistente generada por un LLM es responsabilidad del grupo.

**Ejemplo de declaración de uso de IA en Agradecimientos:**
> *Los autores utilizaron ChatGPT (OpenAI, GPT-4o) para asistir en la redacción y revisión del texto, y GitHub Copilot para completar fragmentos de código repetitivos. Todos los resultados, ecuaciones y conclusiones fueron verificados manualmente por los autores.*
:::

## Posters

La exposición de posters se realizará el **17/06/2026** durante el horario de clase.
Cada grupo presentará su trabajo en formato poster, de manera similar a como se hace en congresos y conferencias científicas.

El poster debe incluir:

- **Título y autores**
- **Introducción y motivación** — el problema abordado y por qué es relevante.
- **Métodos** — descripción concisa del modelo o técnica utilizada, con ecuaciones clave si corresponde.
- **Resultados principales** — figuras y/o tablas que ilustren los hallazgos más importantes.
- **Conclusiones** — síntesis de los resultados y preguntas abiertas.

El poster debe poder leerse de forma autónoma y ser visualmente claro: priorizar figuras y esquemas por sobre texto extenso.
No hay restricciones de formato o herramienta; pueden usar PowerPoint, LaTeX o cualquier otra.

Durante la exposición, el instructor y los demás estudiantes podrán hacer preguntas sobre el trabajo.
Seguramente hayan profesores y estudiantes invitados recorriendo los posters.

## Evaluación

La evaluación final consiste en una nota combinada entre la presentación de posters y el proyecto final.
La nota es individual, aun cuando el proyecto haya sido realizado en grupo.
Se evaluará la claridad de la exposición, la solidez metodológica, la calidad del análisis y la capacidad de responder preguntas sobre el trabajo realizado.
