


# **TuMags Data Reduction Pipeline — Quick Start Guide**

This document describes the basic steps required to run the **TuMags Data Reduction Pipeline**, configure its YAML file, and prepare the Python environment.

---

## **1. Overview**

The **TuMags Data Reduction Pipeline** processes solar observations (timelines) identified by an `obs_ID` as defined in the `config_tumag.yaml` file.  
In the reduction process several reduction levels can be run, from basic calibrations to **phase diversity** restoration.

| Level | Description | Function |
|-------|-------------|----------|
| **0.5** | Flat correction, cropping, and FITS saving | `reduce_image_0_5` |
| **0.7** | Alignment and demodulation | `reduce_image_0_7` |
| **1.0** | Cross-talk correction and normalization | `reduce_image_1_0` |
| **1.1** | Phase diversity restoration | `reduce_image_1_1` |

Output files follow this format:  
`<obsID>*<om_name>*<Nlambda>_<date>*LV*<level>_v<proc_version>.fits`

Level **0.5** requires access to the raw data and a special folder configuration while the rest of the levels can be run provided the previous level is accessible. All levels, except raw, will be accessible at the TuMag database.

---

## **2. Installation and Initial Setup**

### **2.1. Requirements**

- Python ≥ 3.8  
- Packages listed in `requirements.txt`  
- Data organized in the following paths (see below how to configure TuMag in a YAML file):  
  - `tumag_data_location` is the path where the Organized_files_local folder is present (the name of the folder is provided in the `Organized_files_local_folder_name` field). This folder contains the `csv` files with the information of the raw data location. This is not necessary if **0.5** data level is accessible.
- The `obs_ID` must exist in `process_data_timelines.obs_dict` parameter which is defined in `process_data_timelines.py`. This file contains only few timelines. The information required in this file can be found in the [TuMag LOGBOOK web](https://docs.google.com/spreadsheets/d/1RJ5KIgxMN6B-1xDe9gRfoTh1L_uTDMbZw0ajLK6S0so/edit?), although it will be fully upgraded eventually.

### **2.2. Create Conda Environment**

```bash
conda create --name TuMag python=3.10
conda activate TuMag
pip install -r requirements.txt
```

---

## **3. Code Directory Configuration**

Before running the pipeline, **you must modify** the `process_data_main.py` file to include the local directory where the pipeline code is located.  
Look for the import section (at the beginning of the file) and add the following line **adjusting the path to your system**:

```python
# Location of the TuMag software:
import sys
sys.path.append("/Users/orozco/IdAdA Dropbox/David orozco suárez/Python/TuMAG_codes/TuMags_Reduction_Pipeline")
```

⚠️ This line allows Python to find the pipeline’s internal modules, even if executed from a different directory. This behaviour will be change soon.

---

## **4. Main Pipeline Files**

| File | Description |
|------|-------------|
| `process_data_main.py` | Main script; runs all pipeline stages. |
| `config_tumag.yaml` | Configuration file with execution parameters. **copy anywhere to use it**. |
| `process_data_utils.py` | Configuration loader and utilities. |
| `process_data_timelines.py` | Defines `obs_dict` and *timelines*. |
| `fits_files_handling.py` | FITS read/write (`generate_fits`, `update_header`). |
| `master_dark.py`, `master_flatfield.py` | Calibration generation. |
| `alignment.py`, `demodulation.py` | Alignment and demodulation processing. |
| `xtalk_jaeggli.py` | Cross-talk correction routines. |
| `pd_functions_v22.py` | Phase diversity functions. |

---

## **5. YAML Configuration (`config_tumag.yaml`)**

Below is an example of the `config_tumag.yaml` file, which defines all parameters used by the pipeline:

<details>
  <summary>📂 Check YAML config</summary>

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
</details>

The file includes:

- **Global parameters**: `obs_ID`, `proc_version`, `process_line`, paths, parallel settings.
- **Reprocessing options**: Flags to force recalculation of specific levels.
- **Plotting options**: Enable/disable plots for each level.
- **Dark and flat options**: Indexes and normalization regions.
- **Level-specific options**: Parameters for levels 0.5, 0.7, 1.0, and 1.1.
- **Observation modes**: Reference list of available modes.

---

## **6. Running the Pipeline**

Basic execution example:

```bash
python3 process_data_main.py -f config_tumag.yaml
```

> If the `-f` parameter is not specified, the script defaults to `config_tumag.yaml`.  
> You must include the path to the main script.

### **Parallel Execution**

Set in the YAML:

```yaml
parallel: true
max_workers: 4
```

---

## **7. Recommended Workflow**

1. **Edit `config_tumag.yaml`**  
   - Adjust `obs_ID`, paths, and level options.  
   - Review `force_redo` flags for levels to recalculate.

2. **Verify Data Structure**  
   - Data must be organized under `tumag_data_location/Organized_files_local_folder_name`.

3. **Run the Pipeline**

   ```bash
   python3 process_data_main.py -f config_tumag.yaml
   ```

4. **Check Results**  
   - Visualize FITS files with:

     ```bash
     python3 visor.py /Volumes/TuMag_d1/reduccion/.../file.fits
     ```

---

## **8. Common Troubleshooting**

| Problem | Likely Cause / Solution |
|---------|--------------------------|
| **“No files found” or empty OCs** | Check `obs_ID` and `process_line` in `process_data_timelines.obs_dict`. |
| **Need to recalculate flats/darks** | Set `force_redo.redo_flat` or `force_redo.redo_dark` to `true`. |
| **Alignment failure** | Ensure `pinholes_results.csv` matches file dates. |
| **ImportError on execution** | Confirm the `sys.path.append()` path correctly points to the code directory. |
