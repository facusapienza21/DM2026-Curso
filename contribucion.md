# Contribución

El material de esta materia es **abierto**, lo cual quiere decir que cualquier persona puede usarlo para sus propios fines pedagógicos, y **colaborativo**, con lo cual contribuciones de alumnos son bienvenidas.
Eso significa que los estudiantes también pueden ayudar a mejorarlo.

Te invitamos a contribuir en cualquiera de los siguientes casos:
- 🔤 Errores de tipeo u ortográficos
- 🔗 Links que no funcionan
- 📖 Explicaciones poco claras o incorrectas
- 💻 Código que no funciona
- ➕ Material adicional que sería útil

:::{important}
Todos los estudiantes que estén tomando la materia por los puntos de optativa deberán contribuir a los apuntes de la materia al menos una vez.
Más información al final de esta página bajo `Requisito de contribuir a la materia`.
:::


## ¿Por qué fomentar contribuciones?

Además de mejorar el material del curso, esto también es una buena forma de practicar el workflow de colaboración con Git y GitHub, que es una herramienta clave en computación científica.
Este curso intenta acercarse a la forma en que realmente se produce conocimiento científico y software:
- 🤝 Trabajo colaborativo
- 🔍 Revisión por pares
- 🔄 Mejora continua del material

Incluso pequeñas correcciones ayudan muchísimo a mejorar el curso para todos los estudiantes e interesados.
Además, tu nombre quedará para siempre en el repositorio del curso por haber contribuido y aportado cambios al mismo.

## Reportar un problema (Issue)

Si encontrás un problema o querés sugerir una mejora, lo más simple es abrir un issue en GitHub.
Un issue sirve para reportar errores, proponer mejoras, hacer preguntas sobre el material y discutir cambios antes de implementarlos.

:::{note} Cómo abrir un _issue_
Podés encontrar información paso a paso sobre cómo abrir un issue en el siguiente [link](https://docs.github.com/es/issues/tracking-your-work-with-issues/creating-an-issue
). Los issues se deben abrir en el repositorio asociado a la página web de este curso. Para abrir un issue nuevo, podés ir al siguiente [link](https://github.com/facusapienza21/DM2026-Curso/issues/new).
:::

:::{attention} Atención
Antes de abrir un issue nuevo, asegurate de que dicho issue no existe. Podés chequear todos los issues en el siguiente [link](https://github.com/facusapienza21/DM2026-Curso/issues).
:::

## Proponer cambios al material (Pull Request)

Si querés corregir o mejorar directamente el contenido del curso, podés enviar un Pull Request (PR).
Un Pull Request propone cambios al repositorio original y permite que el instructor revise el cambio antes de integrarlo.

El flujo típico para contribuir a un repositorio en GitHub es el siguiente:

### 1️⃣ Crear un fork del repositorio

Ir al repositorio del curso en GitHub y hacer click en fork.
Esto crea una copia del repositorio en tu propia cuenta.
Podés encontrar información para crear un fork en el siguiente [link](https://docs.github.com/es/get-started/quickstart/fork-a-repo).

### 2️⃣ Clonar tu fork y trabajar localmente

Cloná tu repositorio en local y realizá los nuevos cambios.
Recomiendo crear una nueva _branch_ en local para realizar los cambios.

:::{attention} 📚 MyST
Los contenidos de esta página están escritos usando MyST (Markedly Structured Text), una extensión de Markdown utilizada por [Jupyter Book](https://jupyterbook.org/) para crear documentación científica.
Esto permite combinar fácilmente texto, ecuaciones, código y notebooks en un mismo documento, al mismo tiempo que permite tener todos estos contenidos en una página web.
Si querés contribuir al material del curso, en muchos casos solo tendrás que editar archivos `.md` o notebooks `.ipynb`.
Podés encontrar más información sobre cómo editar estos archivos o qué otros atributos tiene MyST en la [documentación](https://mystmd.org/guide/quickstart).
:::

:::{tip} ¿Cómo testear que mis cambios están bien?
Para asegurarte de que los cambios que hiciste no rompen nada y son compatibles con MyST, podés generar una versión local de la página web.
Para aprender cómo hacer esto, seguí las instrucciones del siguiente [link](https://mystmd.org/guide/quickstart#preview-your-myst-site-locally).
:::

### 3️⃣ Guardar los cambios y agregarlos a GitHub

Esta es la secuencia clásica `git add`, `git commit` y `git push` donde se agregan los cambios nuevos a tu versión local del repositorio y luego a la remota (GitHub).

### 4️⃣ Abrir un Pull Request

Desde GitHub (la página web):
- Abrí un PR. Podés encontrar más información de cómo abrir un PR en el siguiente [link](https://docs.github.com/es/pull-requests/collaborating-with-pull-requests/creating-a-pull-request).
- Seleccioná el [repositorio original](https://github.com/facusapienza21/DM2026-Curso) como destino.
- Agregá al instructor como reviewer, así puedo evaluar los cambios y potencialmente dar feedback sobre más cambios.

:::{tip} Buenas prácticas para contribuciones
Algunas recomendaciones:
- Hacé cambios pequeños y específicos por cada commit.
- Al momento de abrir un PR, describí claramente qué cambia tu PR con respecto a la versión original.
- Revisá que el material compile correctamente y la página web no se rompa.
:::

## Requisito de contribuir a la materia

Todos los estudiantes que estén tomando la materia por los puntos de optativa deberán contribuir a los apuntes de la materia al menos una vez.
Para ello, tenemos una planilla donde los estudiantes deben anotar su nombre y anotarse para estar a cargo de una de las clases (la misma está disponible por medio de las cadenas de mails de la materia).
Idealmente, dos estudiantes se deben encargar de las notas de cada clase.

:::{warning}
En caso de no estar familiarizado con cómo contribuir a un repositorio en GitHub, por favor leer cuidadosamente las notas previas.
Todo el contenido necesario y los links útiles se encuentran ahí!
:::

Al momento de abrir un PR para contribuir a la materia, por favor tener en cuenta los siguientes puntos:
- 🤝 Los estudiantes encargados de los apuntes de una clase deben coordinar sus PRs. Hay varias maneras de hacer esto, pero debe haber un único PR por clase con commits de todos los estudiantes. Luego de que el PR se abra, Facu va a encargarse de hacer una revisión y posiblemente sugerir cambios y mejoras. Luego de iterar y haber conseguido una buena versión de las notas, los cambios se van a integrar en la página web de la materia.
- 📝 En las notas, por favor incluir cuanto detalle sea posible sobre la clase. Piensen que las notas tienen que ser entendibles para alguien que no fue a la clase. Pueden usar las grabaciones de la clase para guiarse!
- ✨ No es necesario volverse un experto en myst para poder contribuir a la página web, pero recomiendo que vean la planilla de myst en la página web y cómo la misma fue generada. Myst es super elegante para escribir! Pueden usar cuantas herramientas quieran para que la página web quede lo mejor posible.
- 🖼️ Pueden incluir imágenes o figuras en las notas (de hecho, mucho mejor si lo hacen!). En tal caso, por favor incluir imágenes dentro de la carpeta `images`.
- 💻 En caso de que en la clase se haya mostrado código, por favor incluir en las notas el link al mismo y también las líneas de código mencionadas (myst deja poner código!). En vez de simplemente copiar y pegar el código, tratar de partir el mismo en bloques explicando qué está pasando en el mismo y qué se fue explicando durante la clase.
- 📚 Notar que en las notas se pueden poner palabras clave que luego pueden ser incluidas en el glosario.
- 🏅 Agregar sus nombres junto con sus usuarios de GitHub en la tabla al final del `README`. De esta manera podemos llevar un conteo de quienes ya han contribuido a la materia y también podemos dar crédito por el trabajo de los estudiantes.

¡Así es como se trabaja en GitHub! Bienvenido al mundo del software abierto. 🚀
