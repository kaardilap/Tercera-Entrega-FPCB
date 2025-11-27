README — Análisis Bioinformático y Modelización Biológica

Este repositorio contiene todos los análisis desarrollados durante el curso, incluyendo exploración de datos, análisis multivariado, modelado biológico y simulaciones dinámicas. Cada sección está documentada y acompañada de código en Python.

CONTENIDO

Análisis PCA

Clustering / DSTLM

Modelo Biológico I — Monod Extendido

Modelo Biológico II — Monod con término basal

Modelo Dinámico — Sistema C–N–O–B

ANÁLISIS PCA

Objetivo:
Reducir la dimensionalidad del dataset de ácidos grasos y explorar diferencias entre condiciones experimentales.

Qué se hizo:

Cargar y estandarizar los datos.

Calcular los componentes principales.

Graficar PC1 vs PC2.

Resultados:
El PCA muestra separación entre muestras principalmente asociada a la disponibilidad de nutrientes (Replete vs Limited) y gradientes de luz. Los ácidos grasos más influyentes en PC1 explican la mayor variabilidad del sistema.

CLUSTERING / DSTLM

Objetivo:
Evaluar si las muestras forman grupos naturales usando clustering jerárquico.

Qué se hizo:

Calcular distancias entre perfiles lipídicos.

Construir un dendrograma usando método Ward.

Comparar los grupos con las condiciones experimentales.

Resultados:
El clustering formó grupos coherentes con la condición de nutrientes y en algunos casos con la intensidad de luz. Es complementario al PCA y muestra estructuras jerárquicas en los datos.

MODELO BIOLÓGICO I — MONOD EXTENDIDO

Objetivo:
Ajustar un modelo que relacione la disponibilidad de nutrientes y luz con la producción de PUFA.

Modelo usado:
G(N, L) = Gmax * (N / (KN + N)) * (L / (KL + L))

El modelo tiene tres parámetros biológicos: Gmax, KN y KL.

Problema encontrado:
El ajuste produjo valores muy cercanos a cero para KN y KL. Esto implica saturación instantánea, lo cual no es biológicamente correcto. Además, el modelo predecía valores cercanos a cero en condiciones de baja luz/nutrientes, lo cual no concuerda con la fisiología real del organismo.

MODELO BIOLÓGICO II — MONOD REFORMULADO

Objetivo:
Incorporar restricciones biológicas ausentes en el modelo I y obtener un ajuste fisiológicamente realista.

Cambios implementados:

Se añadió un término basal Gmin:
G(N, L) = Gmin + Gmax * (N / (KN + N)) * (L / (KL + L))

Se impusieron límites mínimos (> 0.01) para KN y KL.

El nuevo modelo evita colapsar los parámetros a cero y mantiene coherencia fisiológica.

Resultados obtenidos:
Gmax = 45.75
KN = 0.10
KL = 2.99
Gmin = 10.00

Interpretación:

KN y KL ya no son irreales; representan concentraciones de semisaturación razonables.

El modelo ahora predice un nivel basal mínimo de PUFA incluso bajo condiciones limitantes.

La gráfica 3D mantiene su forma general, pero la superficie no cae a cero sino a Gmin, lo cual es biológicamente más consistente.

MODELO DINÁMICO — SISTEMA C–N–O–B

Objetivo:
Modelar de manera dinámica cómo la disponibilidad de carbono (C), nitrógeno (N) y oxígeno (O) influye en la biomasa (B) a lo largo del tiempo.

Modelo usado (ecuaciones diferenciales):

dC/dt = I − aC − bBC
dN/dt = J − cN − dBN
dB/dt = eBC − fB

Tres componentes → modelo dinámico válido según los requisitos del curso.

Qué se hizo:

Simular el sistema usando solve_ivp.

Explorar la dinámica para distintos parámetros.

Generar gráficas comparando cómo cambian C, N y B dependiendo de las tasas de entrada y consumo.

Crear una versión con parámetros ajustables para estudiar el comportamiento del sistema bajo distintos escenarios.

Resultados:

Si la entrada de recursos es baja (I y J pequeños), la biomasa se estabiliza en valores bajos.

Si los recursos son altos, la biomasa aumenta hasta alcanzar un equilibrio determinado por su propia tasa de consumo.

El modelo dinámico permite observar comportamientos temporales que los modelos Monod (estáticos) no pueden capturar.

COMPARACIÓN GENERAL ENTRE MODELOS

Modelo Monod I/II (estático):

Predice valores finales de producción de PUFA bajo distintas condiciones.

No incluye tiempo.

Adecuado para describir respuestas fisiológicas de equilibrio.

Modelo dinámico C–N–O–B:

Simula la evolución temporal del sistema.

Incluye recursos y biomasa interactuando.

Permite analizar estabilidad, puntos de equilibrio y sensibilidad a parámetros.

CONCLUSIÓN

El proyecto integra:

Exploración multivariada (PCA, clustering)

Modelos biológicos de respuesta (Monod y su versión corregida)

Modelos dinámicos basados en ecuaciones diferenciales

Cada parte aporta una visión distinta: estructura global de los datos, relación entre factores ambientales y biomasa, y dinámica temporal de los recursos.
