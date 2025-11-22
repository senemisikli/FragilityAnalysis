# Seismic Fragility Analysis (MATLAB)

This repository contains MATLAB scripts for generating seismic fragility curves based on building-level peak ground velocity (PGV) and drift data.  
The workflow implements a power-law seismic demand model and lognormal fragility formulation following Hueste & Bai (2004).

## Features
- PGV–Drift scatter plots (log-scale)
- Regression fitting (individual and group-based)
- Computation of C, k, β_fs, R² and limit-state capacities
- Median PGV estimation (P=0.50)
- Automatic export of:
  - Regression figures  
  - Fragility curves  
  - Summary tables (CSV)

## Output Structure

After running the script, the following folders/files are created:

### **Figures**
- `driftvspgv_figs/`
- `powerlawfit_fits/`
- `fragility_figs/`

### **Tables**
- `Fragility_Summary_MAIN.csv`
- `Median_PGV_Table.csv`
- `Regression_Stats.csv`
- `Group_Average_Stats.csv`

## Requirements
- MATLAB R2021 or newer  
- *Statistics and Machine Learning Toolbox* (required for `normcdf`)  

If the toolbox is not available, the script can be adapted using an alternative CDF formulation.

## Author
**Senem Işıklı**  
Graduate Researcher – METU Civil Engineering  
Project focus: seismic hazard, vulnerability analysis, site response, and sequential earthquakes.

## License
MIT License – see `LICENSE` file for details.
