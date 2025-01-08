# README: Network of Neuron-Astrocyte Mass Models (NNAMM) Simulations

This README file assists in understanding the codes and data accompanying the [manuscript](https://doi.org/10.1371/journal.pcbi.1012683) titled "Dialogue mechanisms between astrocytic and neuronal networks: a whole-brain modelling approach".

## Table of Contents

1. [Requirements](#requirements)
2. [Main Scripts](#main-scripts)
3. [Data](#data)
4. [Data Use Terms](#data-use-terms)
5. [Main Script Details](#main-script-details)
   - [Simulations](#simulations)
   - [Bifurcation Diagrams](#bifurcation-diagrams)
   - [MatCont Analysis](#matcont-analysis)
   - [PyDSTool Analysis](#pydstool-analysis)
   - [Illustrations](#illustrations)
6. [Helper Scripts](#helper-scripts)
7. [Codes and Data Not Included](#codes-and-data-not-included)
8. [Citation](#citation)
9. [Contact Information](#contact-information)

---

## Requirements

To run the scripts, the following software and packages are required:

- **MATLAB** (tested with version R2022a)
- **[MatCont](https://gitlab.utwente.nl/m7686441/matcont)** (tested with version 7.3)
- **Python** (tested with version 3.6.15)
- **[PyDSTool](https://github.com/robclewley/pydstool)** (tested with version 0.91.0)

For detailed methodologies, please refer to the manuscript.

---

## Main Scripts

- [`./codes/nnamm_simulations.m`](./codes/nnamm_simulations.m): Runs the main simulations.
- [`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m): Draws bifurcation diagrams.
- [`./codes/matcont_analysis.m`](./codes/matcont_analysis.m): Performs MatCont analyses.
- [`./codes/pydstool_analysis.py`](./codes/pydstool_analysis.py): Performs PyDSTool analyses.
- [`./codes/illustrations_basic.mlx`](./codes/illustrations_basic.mlx): An interactive MATLAB Notebook for handling data and creating figures.

Refer to [Main Script Details](#main-script-details) for more information.

---

## Data

- **Post-processed neuronal connectome for simulations:**  
  [`./data/connectomes/HCP-average_laus18-3_sc.csv`](./data/connectomes/HCP-average_laus18-3_sc.csv)

- **Post-processed astrocytic connectome for simulations:**  
  [`./data/connectomes/ICBM152A09C_laus18-3_mid-fs-311558V_weights.csv`](./data/connectomes/ICBM152A09C_laus18-3_mid-fs-311558V_weights.csv)

- **Labels for post-processed connectomes:**  
  [`./data/connectomes/hierarchy.laus18-3.lobes.xls`](./data/connectomes/hierarchy.laus18-3.lobes.xls)

- **Raw connectomes and their labels:**  
  [`./data/connectomes/raw/`](./data/connectomes/raw/)

- **Post-processed whole-brain mean data after running calibration simulations:**  
  [`./data/analyses/activity/papaya-r0/activity.whole-brain.mean.xls`](./data/analyses/activity/papaya-r0/activity.whole-brain.mean.xls)

- **Post-processed whole-brain mean data after running main simulations:**  
  [`./data/analyses/activity/papaya-r1/activity.whole-brain.mean.xls`](./data/analyses/activity/papaya-r1/activity.whole-brain.mean.xls)

- **Post-processed limit cycle curves data from PyDSTool analyses:**  
  [`./data/bifurcation-diagrams/alpha/C_Pyr_to_Pyr-7.5.lcc.mat`](./data/bifurcation-diagrams/alpha/C_Pyr_to_Pyr-7.5.lcc.mat)

---

## Data Use Terms

This repository includes a structural connectome derived from the Human Connectome Project (HCP) Young Adult dataset.
Users of this data must comply with the [HCP Data Use Terms](./data/connectomes/raw/DataUseTerms-HCP-Open-Access-26Apr2013.pdf).
For more information, see the [HCP Young Adult Data Use Terms](https://www.humanconnectome.org/study/hcp-young-adult/data-use-terms).

---

## Main Script Details

### Simulations

[`./codes/nnamm_simulations.m`](./codes/nnamm_simulations.m)

- Performs network simulations using dependencies located in [`./codes/util/`](./codes/util/), [`./codes/functions/`](./codes/functions/), and [`./codes/nnamm/`](./codes/nnamm/). The variable `proj_dir` (line 12) specifies these dependencies.
- Utilises post-processed connectomes from [`./data/connectomes/`](./data/connectomes/) (see line 38).
- Outputs are saved relative to `proj_dir` (see line 29).
- Does not check available RAM or disk memory before running. Guidance for configuration is provided (see line 50).

### Bifurcation Diagrams

[`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m)

- Quickly draws bifurcation diagrams.
- Depends on files in [`./codes/util/`](./codes/util/) and [`./codes/functions/`](./codes/functions/). Ensure the variable `proj_dir` (line 6) is correctly specified.
- Outputs are saved relative to `proj_dir` (see line 15).

### MatCont Analysis

[`./codes/matcont_analysis.m`](./codes/matcont_analysis.m)

- Continues saddle-node and Hopf branches using MatCont.
- Dependencies are located in [`./codes/util/`](./codes/util/) and [`./codes/functions/`](./codes/functions/). The variable `proj_dir` (line 9) specifies these dependencies.
- The MatCont library is not included in [`./codes/`](./codes/) and must be provided (see link in [Codes and Data Not Included](#codes-and-data-not-included)). The variable `matcont_bin` (line 23) specifies the path to the root directory of the MatCont library (which contains the function `init.m`).
- Requires outputs from [`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m) as inputs (see line 34).
- Outputs are saved in the same directory as the input files with consistent naming.

### PyDSTool Analysis

[`./codes/pydstool_analysis.py`](./codes/pydstool_analysis.py)

- Continues limit cycles using PyDSTool (see link in [Codes and Data Not Included](#codes-and-data-not-included)).
- Uses outputs from [`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m) as inputs (see line 159).
- Outputs are saved in the same directory as the input files with consistent naming.

### Illustrations

[`./codes/illustrations_basic.mlx`](./codes/illustrations_basic.mlx)  
[`./codes/illustrations_basic.html`](./codes/illustrations_basic.html)

- Guides users in handling data and creating figures (e.g., structural connectomes, simulation outputs, and bifurcation data).
- Relies on files in [`./data/`](./data/), [`./codes/util/`](./codes/util/), and [`./codes/functions/`](./codes/functions/). Configure `proj_dir` (line 1) for dependencies.

---

## Helper Scripts

- [`./codes/nnamm/`](./codes/nnamm/)
- [`./codes/functions/`](./codes/functions/)
- [`./codes/util/`](./codes/util/)

---

## Codes and Data Not Included

- Brain Connectivity Toolbox, 2019-03-03 release, [link](https://www.nitrc.org/projects/bct)
- Connectome Mapper 3, version 3.0.3, [link](https://github.com/connectomicslab/connectomemapper3)
- FreeSurfer, version 6.0.0, [link](https://github.com/freesurfer/freesurfer)
- MatCont, version 7.3, [link](https://gitlab.utwente.nl/m7686441/matcont)
- MuxViz, version 3.1, [link](https://github.com/manlius/muxViz)
- PyDSTool, version 0.91.0, [link](https://github.com/robclewley/pydstool)
- Scientific Colour Maps, version 7.0.1, [link](https://zenodo.org/records/5501399)
- Scilpy Python library, version 1.3.0, [link](https://github.com/scilus/scilpy)
- SET, version 1.0, [link](https://set-documentation.readthedocs.io)
- Tractoflow, version 2.2.1, [link](https://github.com/scilus/tractoflow)

---

## Citation

If you use this code or data in your research, please cite:

Ali et al., 2025, "Dialogue mechanisms between astrocytic and neuronal networks: a whole-brain modelling approach", *PLOS Computational Biology*, DOI: [10.1371/journal.pcbi.1012683](https://doi.org/10.1371/journal.pcbi.1012683)

---

## Contact Information

For questions or issues, contact Obaï Bin Ka'b Ali at [ali.obaibk (at) gmail (dot) com].
