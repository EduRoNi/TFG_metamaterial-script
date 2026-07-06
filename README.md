# TFG_metamaterial-script

Código desarrollado para la caracterización mecánica homogénea de estructuras de metamaterial mediante simulación por Elementos Finitos (ANSYS) y automatización en Python. El repositorio incluye el conjunto de scripts necesarios para calcular la convergencia del RVE (*Representative Volume Element*), las propiedades elásticas y plásticas del material, y la superficie de fluencia (criterio de Hill).

Este código forma parte del Trabajo de Fin de Grado (TFG) *"Caracterización del comportamiento elastoplástico de metamateriales: scripting de modelos de elementos finitos para el cálculo y generación de datos para el uso de IA"*.

## Tabla de contenidos

- [Requisitos](#requisitos)
- [Estructura del repositorio](#estructura-del-repositorio)
- [Descripción de los módulos](#descripción-de-los-módulos)
- [Instrucciones de uso](#instrucciones-de-uso)
- [Resultados](#resultados)

## Requisitos

### Software

| Componente | Versión |
|---|---|
| ANSYS | ≥ 2024 R1 (con licencia de *solver*) |
| Python | 3.12.10 |

### Librerías nativas de Python

- `os` — Gestión de rutas de archivos y creación automática de directorios.
- `sys` — Conexión de variables locales hacia el lanzamiento de otros códigos.
- `subprocess` — Permite lanzar otros códigos desde `Main.py`.
- `shutil` — Operaciones avanzadas con archivos (copiado, movimiento y borrado de directorios).
- `time` — Monitorización y registro de los tiempos de cómputo de las simulaciones.

### Librerías externas


- `numpy` — Operaciones matriciales para el tratamiento de datos estructurales.
- `pandas` — Estructuración de datos y automatización de la lectura/escritura en Excel. 
- `ansys.mapdl.core` — Interfaz PyMAPDL para la conexión, control y *scripting* directo con el *solver* de ANSYS. 
- `matplotlib.pyplot` — Generación y exportación de gráficas (curvas *stress-strain*, superficies de fluencia). 
- `h5py` — Almacenamiento y manipulación de datos masivos en formato HDF5. 
- `skimage` — Procesamiento de imágenes 3D de la superficie de fluencia. 

## Estructura del repositorio

Los *scripts* principales deben coexistir en la misma carpeta de trabajo:

```
Directorio_Raiz/
├── Main.py
├── Funciones.py
├── Calculo_convergencia_uniaxial.py
├── Calculo_cortante_elastico.py
├── Calculo_uniaxial_plastico.py
├── Calculo_cortante_plastico.py
├── Criterio_fallo.py
└── mi_modelo_celda.x_t

temporal/

Carpeta_Resultados/
├── Base_Datos_Resultados.xlsx
└── Gráficas y figuras
```

- **`Directorio_Raiz/`**: contiene el código y el modelo geométrico de la celda de estudio.
- **`temporal/`**: carpeta vacía utilizada para los archivos intermedios generados durante las simulaciones (se eliminan automáticamente).
- **`Carpeta_Resultados/`**: carpeta de salida donde se almacenan los resultados finales (Excel, gráficas y figuras).

## Descripción de los módulos

| Archivo / Módulo | Función principal | Resultados / Salidas |
|---|---|---|
| `Main.py` | Control del flujo global. | — |
| `Funciones.py` | Recopilación de funciones utilizadas. | — |
| `Calculo_convergencia_uniaxial.py`* | Convergencia del tamaño del RVE con ensayo elástico uniaxial. | Dimensiones de RVE, $E$ y $\nu$ por dirección. |
| `Calculo_cortante_elastico.py`* | Ensayo elástico a cortante puro. | Módulo $G$. |
| `Calculo_uniaxial_plastico.py`* | Ensayos uniaxiales en régimen plástico. | Curvas $\sigma$-ε y F-u, límites elásticos y estructurales. |
| `Calculo_cortante_plastico.py`* | Ensayos a cortante en régimen plástico. | Curvas $\tau$-γ y F-u a cortante, límites elásticos y estructurales. |
| `Criterio_fallo.py` | Cálculo de la superficie de fluencia (criterio de Hill). | Superficies de fluencia 3D y proyecciones. |

> \* Cada uno de estos módulos cuenta con una versión adicional para *Modelos Reales*, que incluye modificaciones por imperfección geométrica.

## Instrucciones de uso

Gracias al alto grado de automatización del algoritmo, el usuario solo necesita interactuar con el archivo principal `Main.py`. Antes de ejecutar cualquier simulación:

1. **Preparar el modelo geométrico.**
   Guardar el archivo de la celda de estudio del metamaterial en formato Parasolid (`.x_t`) dentro de la misma carpeta de trabajo. Si el mallado automático falla, es necesario incluir un archivo `ds.cdb` con la malla previamente generada, externo al algoritmo.

2. **Configurar las variables en `Main.py`.**
   Modificar únicamente los siguientes parámetros de inicialización:

   - `file_name`: nombre exacto del archivo del modelo, entre comillas y sin extensión (ej. `'mi_modelo_celda'`).
   - `reports_dir`: ruta del directorio donde el código creará una carpeta con el nombre `file_name`, en la que se guardarán todos los resultados (archivos Excel e imágenes de las distintas gráficas).
   - `temporal`: ruta de una carpeta vacía utilizada para los archivos intermedios de las simulaciones, que se eliminarán al finalizar.

3. **Ejecutar el código.**
   Con todo lo anterior configurado, basta con ejecutar `Main.py` para lanzar el flujo completo de simulación y post-procesado.

## Resultados

Al finalizar la ejecución, `Carpeta_Resultados/` contendrá:

- `Base_Datos_Resultados.xlsx`: valores numéricos de todos los ensayos (propiedades elásticas, curvas tensión-deformación, límites estructurales, etc.).
- Gráficas y figuras exportadas (curvas *stress-strain*, superficies de fluencia 3D y sus proyecciones).


## Autor

Eduardo Rodríguez Nieves
