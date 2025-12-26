# Gen-Searcher-
Repositorio con el código fuente y los conjuntos de datos generados durante el período de prácticas externas del Máster en Bioinformática y Biología Computacional.
# DESCRIPCIÓN DE LA CARPETA



Esta carpeta contiene varias libretas de código (RMarkdown) dedicadas a la descarga, procesamiento y análisis funcional de datos provenientes de diferentes bases de datos públicas.



##### NOTEBOOKS INCLUIDAS

* **BioPortal.Rmd:** describe el protocolo utilizado para descargar datos de proyectos alojados en *BioPortal* y cómo combinarlos con otros tipos de datos.
* **GEO.Rmd:** explica el proceso de descarga de datos desde *GEO*, tanto mediante código como de forma manual (útil cuando la descarga automatizada falla).
* **UCSC\_Xena.Rmd** y **UCSC\_Xena\_doble\_cop.Rmd**: detallan el protocolo seguido para descargar y procesar datos provenientes de *UCSC Xena.*



##### PROCESAMIENTO Y ANÁLISIS POSTERIOR



En las cuatro libretas se realizan los siguientes pasos después de la descarga:



* Cálculo de expresión diferencial entre distintas comparaciones (según lo permitan los datos).
* Cálculo de estadísticas descriptivas para cada condición: media, mediana y desviación estándar.
* Análisis funcional usando Reactome y Gene Ontology para identificar las rutas asociadas a los genes diferencialmente expresados.
* Los resultados finales se almacenan en archivos .rds.



##### UNIFICACIÓN DE DATOS



La notebook **Unificar\_data.Rmd** describe el procedimiento para combinar y condensar archivos RDS de un mismo tipo de cáncer. En el ejemplo incluido, se unifican datos provenientes de UCSC Xena y GEO.



##### CONSIDERACIONES IMPORTANTES



Los notebooks fueron desarrollados usando un *dataset* de melanoma. Si se analizan otros tipos de cáncer, será necesario modificar:

* Nombres de variables
* Nombres de columnas
* La definición de los contrastes



Estos cambios están explicados dentro de cada notebook.



### CARPETA: Gene\_Searcher



Esta carpeta contiene la herramienta Gene\_Searcher, compuesta por:



* **App**: Script que debe ejecutarse para iniciar la interfaz gráfica de Gene\_Searcher.
* **data\_tipo\_de\_cancer.rds**: Archivos RDS con información procesada para los siguientes cánceres: Melanoma (Melanoma\_vs\_Metástasis, Tipo\_Metástasis\_vs\_Tumor, Metástasis\_vs\_Metástasis), Mama e Hígado



###### Gene Searcher:





La app de Gene Searcher tiene tres pestañas principales:



* **Visualización**: esta pestaña permite seleccionar un tipo de cáncer y buscar un gen concreto. Si el gen está presente en el conjunto de datos, se genera un violin plot que muestra su nivel de expresión entre las condiciones seleccionadas. Además, se incluye información funcional del gen, como las rutas de Reactome asociadas y los términos de Gene Ontology (GO).
* **Visualización Global**: esta pestaña permite seleccionar un tipo de cáncer y muestra un volcano plot interactivo con los genes sobreexpresados (up-regulados) y subexpresados (down-regulados). También presenta las rutas de Reactome más representadas en el análisis.
* Ayuda: pestaña explicativa de la App así como los links a las bases de datos empleadas.



































