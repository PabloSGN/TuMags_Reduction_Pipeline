# TuMags Reduction Pipeline — Quick Start Guide

Este documento describe los pasos necesarios para ejecutar el **TuMags Reduction Pipeline**, configurar su archivo YAML y preparar el entorno Python.

---

## 1. Descripción general

El *TuMags Reduction Pipeline* procesa observaciones solares identificadas mediante un `obs_ID` definido en el archivo `config_tumag.yaml`.  
Cada ejecución produce varios niveles de reducción, desde calibraciones básicas hasta restauración por *phase diversity*.

| Nivel | Descripción | Función |
|:------|:-------------|:--------|
| **0.5** | Corrección de *flat*, recorte y guardado FITS | `reduce_image_0_5` |
| **0.7** | Alineación y demodulación | `reduce_image_0_7` |
| **1.0** | Corrección de *cross-talk* y normalización | `reduce_image_1_0` |
| **1.1** | Restauración por *phase diversity* | `reduce_image_1_1` |

Los archivos de salida tienen este formato: 
`<obsID>*<om_name>*<Nlambda>_<date>*LV*<level>_v<proc_version>.fits`

---

## 2. Instalación y configuración inicial

### 2.1. Requisitos
- Python ≥ 3.8  
- Paquetes indicados en `requirements.txt`  
- Datos organizados en las rutas (ver más abajo como configurar TuMag en un archivo Yalm):
  - `cfg.tumag_data_location`
  - `cfg.Organized_files_local_folder_name`
- El `obs_ID` debe existir en `process_data_timelines.obs_dict`

### 2.2. Crear entorno con Conda
```bash
conda create --name TuMag python=3.10
conda activate TuMag
pip install -r requirements.txt
````

---

## 3. Configuración del directorio del código

Antes de ejecutar el pipeline, **es necesario modificar** el archivo `process_data_main.py` para incluir el directorio local donde se encuentra el código del pipeline.

Busca la sección de importación (al principio del archivo) y añade la siguiente línea **ajustando la ruta a tu sistema**:

```python
# Location of the TuMag software:
import sys
sys.path.append("/Users/orozco/IdAdA Dropbox/David orozco suárez/Python/TuMAG_codes/TuMags_Reduction_Pipeline")
```

⚠️ Esta línea permite que Python encuentre los módulos internos del pipeline, incluso si se ejecuta desde un directorio distinto al del código fuente.

---

## 4. Archivos principales del pipeline

| Archivo                                 | Descripción                                                   |
| :-------------------------------------- | :------------------------------------------------------------ |
| `process_data_main.py`                  | Script principal; ejecuta todas las etapas del pipeline.      |
| `config_tumag.yaml`                     | Archivo de configuración con parámetros de ejecución.         |
| `process_data_utils.py`                 | Cargador de configuración y utilidades.                       |
| `process_data_timelines.py`             | Define `obs_dict` y *timelines*.                              |
| `fits_files_handling.py`                | Lectura/escritura de FITS (`generate_fits`, `update_header`). |
| `master_dark.py`, `master_flatfield.py` | Generación de calibraciones.                                  |
| `alignment.py`, `demodulation.py`       | Procesamiento de alineación y demodulación.                   |
| `xtalk_jaeggli.py`                      | Rutinas de corrección de *cross-talk*.                        |
| `pd_functions_v22.py`                   | Funciones de *phase diversity*.                               |

---

## 5. Configuración YAML (`config_tumag.yaml`)

A continuación se muestra el contenido de ejemplo del archivo `config_tumag.yaml`, que define todos los parámetros usados por el pipeline.

```yaml
# ==========================================================
# TuMags Reduction Pipeline Configuration (YAML)
# ==========================================================
# Example:
#   python3 ../../TuMags_Reduction_Pipeline/process_data_main.py \
#       -f ../../TuMags_Reduction_Pipeline/config_tumag.yaml
#
# IMPORTANT:
# - Parameter names must be unique (no duplicates across levels).
# - Boolean and None keywords are lowercase: false, true, none
# - Results can be checked using:
#     python3 visor.py /Volumes/TuMag_d1/reduccion/.../file.fits
# ==========================================================

# ----------------------------------------------------------
# GLOBAL PARAMETERS
# ----------------------------------------------------------
obs_ID: "06_SPOT_AR1"                   # Observation ID
proc_version: "0.x"                     # Processing version (standard cross-talk)
process_line: "2.02"                    # Observation mode to process
process_ocs: ":"                        # ":" -> all, or range, or single ID
process_files: ":"                      # ":" -> all, or range, or single ID
parallel: true
max_workers: 3
tumag_data_location: "./"               # Path to Organized_files_local folder
Organized_files_local_folder_name: "Organized_files_local"
output_folder: "/Volumes/TuMag_d1/reduccion/"   # Output results folder

# ----------------------------------------------------------
# REPROCESSING OPTIONS
# ----------------------------------------------------------
force_redo:
  redo_flat: false
  redo_dark: false
  level0_5: false
  level0_7: true
  level1_0: false
  level1_1: false

# ----------------------------------------------------------
# PLOTTING OPTIONS
# ----------------------------------------------------------
plots:
  plot_darks: false
  plot_flats: false
  plot_level0_5: true
  plot_level0_7: true
  plot_level1_0: true
  plot_level1_1: true
  roi_plots: [200, -200, 200, -200]

# ----------------------------------------------------------
# DARK AND FLAT OPTIONS
# ----------------------------------------------------------
darks_indexes:
  dark_from:
  dark_to: -2

flat_norm_roi: [250, -250, 250, -250]  # Suggested by Fbailen

# ----------------------------------------------------------
# LEVEL 0.5 OPTIONS
# ----------------------------------------------------------
level_05:
  filtering: true
  save_level_05: true

size: 2016
centro: [973, 996]
corte: 800

# ----------------------------------------------------------
# LEVEL 0.7 OPTIONS
# ----------------------------------------------------------
level_07:
  align_mode: "fourier"
  align_roi: [600,1100,100,600]
  align_accuracy: 0.001
  align_rot_data_filter: "pinholes_results.csv"
  align_quadrants: 0
  align_order: 0

# ----------------------------------------------------------
# LEVEL 1.0 OPTIONS
# ----------------------------------------------------------
level_10:
  crosst_mode: "standard"
  crosst_roi: [600,1100,100,600]
  crosst_region: [600,1100,100,600]
  crosst_threshold: 0.01
  crosst_intensity_threshold: 0.8
  crosst_last_wave: -1
  crosst_quadrants: 0
  crosst_interference: 0
  mmatrix: 0
  normalization: 0
  plot_crosst_method: false
  crosst_verbose: false
  add_level_10_label: ''

# ----------------------------------------------------------
# LEVEL 1.1 OPTIONS
# ----------------------------------------------------------
zernike_id: "06_SPOT_Fe2.02_0"

# ----------------------------------------------------------
# OBSERVATION MODES (REFERENCE)
# ----------------------------------------------------------
# observation_modes = {
#     0   : "0s",
#     1   : "0p",
#     2   : "1",
#     3   : "2.02",
#     4   : "2.06",
#     5   : "3.02",
#     6   : "3.06",
#     7   : "4",
#     8   : "5.02",
#     9   : "5.06",
#     64  : "PD_calibration",
#     65  : "Spectral_calibration",
#     66  : "Polarimetric_calibration",
#     128 : "Test"
# }
```

---

## 6. Ejecución del pipeline

Ejemplo de ejecución básica:

```bash
python3 process_data_main.py -f config_tumag.yaml
```

> Si el parámetro `-f` no se especifica, el script usa `config_tumag.yaml` por defecto.
> hay que incluir la ruta del program aprincipal

### Ejecución paralela

Configura en el YAML:

```yaml
parallel: true
max_workers: 4
```

---

## 7. Flujo de trabajo recomendado

1. **Editar `config_tumag.yaml`**

   * Ajustar `obs_ID`, rutas y opciones de niveles.
   * Revisar flags de `force_redo` según el nivel que se quiera recalcular.

2. **Verificar estructura de datos**
   Los datos deben estar organizados bajo `tumag_data_location/Organized_files_local_folder_name`.

3. **Ejecutar el pipeline**

   ```bash
   python3 process_data_main.py -f config_tumag.yaml
   ```

4. **Comprobar resultados**
   Visualizar archivos FITS con:

   ```bash
   python3 visor.py /Volumes/TuMag_d1/reduccion/.../file.fits
   ```

---

## 8. Solución de problemas comunes

| Problema                          | Causa probable / Solución                                                                              |
| :-------------------------------- | :----------------------------------------------------------------------------------------------------- |
| **“No files found” o OCs vacíos** | Revisar `obs_ID` y `process_line` en `process_data_timelines.obs_dict`.                                |
| **Recalcular flats/darks**        | Activar `force_redo.redo_flat` o `force_redo.redo_dark`.                                               |
| **Fallo de alineación**           | Asegurarse de que el CSV de rotación (`pinholes_results.csv`) coincide con las fechas de los archivos. |
| **ImportError al ejecutar**       | Confirmar que la ruta añadida en `sys.path.append()` apunta correctamente al directorio del código.    |

