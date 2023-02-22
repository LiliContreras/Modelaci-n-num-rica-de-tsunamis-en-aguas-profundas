# Modelacion-numerica-de-tsunamis-en-aguas-profundas
Programas empleados en la modelación numérica de tsunamis, modelos batimétricos y scripts de condiciones iniciales para la simulación de modelos de aguas profundas en Geoclaw.

DesplazamientoVertical_Okada_Geoclaw:
  Se muestra de manera esquemática el funcionamiento del archivo "maketopo.py" empleado en Geoclaw para crear un archivo de condiciones iniciales que corresponden a la fuente tsunamigénica. Para su correcto funcionamiento, se debe tener instalado Clawpack. 
  El archivo de salida contenido igualmente en este repositorio, corresponde al desplazamiento vertical usado como condición inicial en el método de ecuaciones integrales (DespVertical.txt).

Archivos de condiciones iniciales para llevar a cabo la simulación en Geoclaw:
  - dtopo_ : archivo de condiciones iniciales correspondiente a la fuente
  - Modelo 1 y batimetría: representan los modelos batimétricos correspondientes a los casos 1 y 2.

Ecuaciones Integrales:

  main: EcIntegrales
  
        - Se elige el modelo a sobre el que se aplicará el método y los puntos de observación en superficie.
        - Se elige cómo se presentan las condiciones (0-Imagenes, 1-condición P=0, 2-condición de Gravedad)
        -Se elige la fuente tsunamigénica (1-Desplazamiento horizontal, 2-Desplazamiento vertical unitario, 3- Desplazamiento vertical con método de Okada)
        - De elegir el desplazamiento horizontal como fuente, se da la opción de obtener las presiones en el contacto lateral o los desplazamientos en superficie.
        
  Funciones:
  
        - refine_boundarymesh: se crea refinamiento en el modelo batimétrico, considerando si se usa el método de imágenes o las condiciones en superficie.
        - build_image: se construye la imágen del modelo en el caso de usar ese método.
        -build_solve_ie: se construye y resuelve el sistema de ecuaciones considerando el método elegido (imágenes, P=0 o condición de gravedad)
        - compute_unormal: Cálculo de presiones y desplazamiento vertical sobre los puntos de observación.
        - Westergaard: Se calcula presión sobre pared lateral del modelo con ayuda de un método analítico.
        
  Archivos de texto:
  
        - DespVertical.txt: Contiene el desplazamiento vertical obtenido con el modelo de okada dada una falla previamente definida. Este resultado se obtuvo con la función de Okada del software Geoclaw.
  
Archivos .mat:

        - Contienen los resultados de la aplicación del método de ecuaciones integrales para el cálculo de presiones sobre una pared vertical, bajo tres consideraciones de aplicación de condiciones iniciales.
