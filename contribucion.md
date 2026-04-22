# Contribución

El material de esta materia es **abierto**, lo cual quiere decir que cualquier persona puede usarlo para sus propios fines pedagógicos, y **colaborativo**, con lo cual contribuciones de alumnos son bienvenidas.
Eso significa que los estudiantes también pueden ayudar a mejorarlo.

Te invitamos a contribuir en cualquiera de los siguientes casos:
- Errores de tipeo u ortográficos
- Links que no funcionan
- Explicaciones poco claras o incorrectas
- Código que no funciona
- Material adicional que sería útil


## ¿Por qué fomentar contribuciones?

Además de mejorar el material del curso, esto también es una buena forma de practicar el workflow de colaboración con Git y GitHub, que es una herramienta clave en computación científica.
Este curso intenta acercarse a la forma en que realmente se produce conocimiento científico y software:
- Trabajo colaborativo
- Revisión por pares
- Mejora continua del material

Incluso pequeñas correcciones ayudan muchísimo a mejorar el curso para todos los estudiantes e interesados.
Además, tu nombre quedara para siempre en el repositorio del curso por haber contribuido y aportado cambios al mismo.

## Reportar un problema (Issue)

Si encontrás un problema o querés sugerir una mejora, lo más simple es abrir un issue en GitHub.
Un issue sirve para reportar errores, proponer mejoras, hacer preguntas sobre el material y discutir cambios antes de implementarlos.

:::{note} Cómo abrir un _issue_
Podes encontrar información paso a paso sobre como abrir un issue en el siguiente [link](https://docs.github.com/es/issues/tracking-your-work-with-issues/creating-an-issue
). Los issues se deben abrir en el repositorio asociado a la página web de este curso. Para abrir un issue nuevo, podes ir al siguiente [link](https://github.com/facusapienza21/DM2026-Curso/issues/new).
:::

:::{attention} Atención
Antes de abrir un issue nuevo, asegurate que dicho issue no existe. Podes chequear todos los issues en el siguiente [link](https://github.com/facusapienza21/DM2026-Curso/issues).
:::

## Proponer cambios al material (Pull Request)

Si querés corregir o mejorar directamente el contenido del curso, podés directamente enviar un Pull Request (PR).
Un Pull Request propone cambios al repositorio original y permite que el instructor revise el cambio antes de integrarlo.

El flujo típico para contribuir a un repositorio en GitHub es el siguiente:

### 1️⃣ Crear un fork del repositorio

Ir al repositorio del curso en GitHub y hacer click en fork.
Esto crea una copia del repositorio en tu propia cuenta.
Podes encontrar información crear un fork en el siguiente [link](https://docs.github.com/es/get-started/quickstart/fork-a-repo).

### 2️⃣ Clonar tu fork y trabajar localmente

Cloná tu repositorio en local y realiza nuevos cambios en local.
Recomiendo crear una nueva _branch_ en local para realizar los cambios.

:::{attention} 📚 MyST
Los contenidos de esta página están escritos usando MyST (Markedly Structured Text), una extensión de Markdown utilizada por [Jupyter Book](https://jupyterbook.org/) para crear documentación científica.
Esto permite combinar fácilmente texto, ecuaciones, código y notebooks en un mismo documento, al mismo tiempo que permite tener todos estos contenidos en una página web.
Si querés contribuir al material del curso, en muchos casos solo tendrás que editar archivos `.md` o notebooks `.ipynb`.
Podes encontrar más información sobre como editar estos archivos o qué otros atributos tiene MyST en la [documentación](https://mystmd.org/guide/quickstart).
:::

:::{tip} Cómo testear que mis cambios estan bien?
Para asegurarte que los cambios que hiciste no rompen nada y son compatibles con MyST, podes generar una versión local de la página web.
Para aprender como hacer esto, seguí las intrucciones del siguiente [link](https://mystmd.org/guide/quickstart#preview-your-myst-site-locally).
:::

### 3️⃣ Guardar los cambios y agregarlos a GitHub

Este es la secuencia clásica `git add`, `git commit` y `git push` donde se agregan los cambios nuevos a tu versión local del repositorio y luego remota (GitHub).

### 4️⃣ Abrir un Pull Request

Desde GitHub (la página web):
- Abrí un PR. Podes encontrar más información de como abrir un PR en el siguiente [link](https://docs.github.com/es/pull-requests/collaborating-with-pull-requests/creating-a-pull-request).
- Seleccioná el [repositorio original](https://github.com/facusapienza21/DM2026-Curso) como destino.
- Agregá al instructor como reviewer, asi puedo avaluar los cambios y potencialmente dar feedback sobre más cambios.

:::{tip} Buenas prácticas para contribuciones
Algunas recomendaciones:
- Hacé cambios pequeños y específicos por cada commit.
- Al momento de abrir un PR, describí claramente qué cambia tu PR con respecto a la versión original.
- Revisá que el material compile correctamente y la página web no se rompa.
:::
