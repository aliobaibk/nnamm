# README: Network of Neuron-Astrocyte Mass Models (NNAMM) Simulations

This README file assists in understanding the codes and data accompanying the manuscript titled "Dialogue mechanisms between astrocytic and neuronal networks: a whole-brain modelling approach".

## Citation

If you use this code or data in your research, please cite the following paper:

Ali et al., 2024, "Dialogue mechanisms between astrocytic and neuronal networks: a whole-brain modelling approach", PLOS Computational Biology.

## Requirements

To run the scripts, the following software and packages are required:

- **MATLAB** (tested with version R2022a)
- **[MatCont](https://gitlab.utwente.nl/m7686441/matcont)** (tested with version 7.3)
- **Python** (tested with version 3.6.15)
- **[PyDSTool](https://github.com/robclewley/pydstool)** (tested with version 0.91.0)

For detailed information and methodologies, please refer to the manuscript.

## Main Scripts

- [`./codes/nnamm_simulations.m`](./codes/nnamm_simulations.m)
- [`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m)
- [`./codes/matcont_analysis.m`](./codes/matcont_analysis.m)
- [`./codes/pydstool_analysis.py`](./codes/pydstool_analysis.py)
- [`./codes/illustrations_basic.mlx`](./codes/illustrations_basic.mlx)

## Helper Scripts

- [`./codes/nnamm/`](./codes/nnamm/)
- [`./codes/functions/`](./codes/functions/)
- [`./codes/util/`](./codes/util/)

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

## Main Scripts Details

### Simulations

[`./codes/nnamm_simulations.m`](./codes/nnamm_simulations.m)

- This script details how the network simulations were done.
- This script has dependencies in [`./codes/util/`](./codes/util/), [`./codes/functions/`](./codes/functions/), and [`./codes/nnamm/`](./codes/nnamm/). The variable `proj_dir` on line 12 specifies these dependencies.
- This script also uses the post-processed connectomes in [`./data/connectomes/`](./data/connectomes/) (see line 31).
- By default, the outputs of this script are saved relative to the directory `proj_dir` (see line 22).
- The script does not check available RAM or disk memory before running. Some guidance is given in the script. See also line 43 to configure the simulations.

### Bifurcation Diagrams

[`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m)

- This script helps to quickly draw bifurcation diagrams.
- This script has dependencies in [`./codes/util/`](./codes/util/) and [`./codes/functions/`](./codes/functions/).
- Before running this script, line 6 (the variable `proj_dir`) must be edited to correctly specify the dependencies.
- By default, outputs are saved relative to the directory `proj_dir` (see line 7).

### MatCont Analysis

[`./codes/matcont_analysis.m`](./codes/matcont_analysis.m)

- This script uses MatCont to continue saddle-node and Hopf branches.
- This script has dependencies in [`./codes/util/`](./codes/util/) and [`./codes/functions/`](./codes/functions/). The variable `proj_dir` on line 9 specifies these dependencies.
- The MatCont library is not included in [`./codes/`](./codes/) and must be provided (see link below). The variable `matcont_bin` on line 16 specifies the path to the root directory of the MatCont library (which contains the function `init.m`).
- This script uses some outputs of [`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m) as inputs (see line 27).
- By default, the outputs of this script are saved in the same directory as the input files, using names consistent with input file names.

### PyDSTool Analysis

[`./codes/pydstool_analysis.py`](./codes/pydstool_analysis.py)

- This script uses PyDSTool (see link below) to continue limit cycles.
- This script uses some outputs of [`./codes/bd_draw_alpha.m`](./codes/bd_draw_alpha.m) as inputs (see line 159).
- By default, the outputs of this script are saved in the same directory as the input files, using names consistent with input file names.

### Illustrations

[`./codes/illustrations_basic.mlx`](./codes/illustrations_basic.mlx)  
[`./codes/illustrations_basic.html`](./codes/illustrations_basic.html)

- This script shows how to handle different data (e.g., structural connectomes, simulation data, bifurcation data) to perform different plots as in the manuscript.
- This script has dependencies in [`./data/`](./data/), [`./codes/util/`](./codes/util/) and [`./codes/functions/`](./codes/functions/). The variable `proj_dir` on line 1 specifies these dependencies.

## Codes and Data Not Included

- Brain Connectivity Toolbox, 2019-03-03 release, [Brain Connectivity Toolbox](https://www.nitrc.org/projects/bct)
- Connectome Mapper 3, version 3.0.3, [Connectome Mapper 3](https://github.com/connectomicslab/connectomemapper3)
- FreeSurfer, version 6.0.0, [FreeSurfer](https://github.com/freesurfer/freesurfer)
- MatCont, version 7.3, [MatCont](https://gitlab.utwente.nl/m7686441/matcont)
- MuxViz, version 3.1, [MuxViz](https://github.com/manlius/muxViz)
- PyDSTool, version 0.91.0, [PyDSTool](https://github.com/robclewley/pydstool)
- Scientific Colour Maps, version 7.0.1, [Scientific Colour Maps](https://zenodo.org/records/5501399)
- Scilpy Python library, version 1.3.0, [Scilpy](https://github.com/scilus/scilpy)
- SET, version 1.0, [SET Documentation](https://set-documentation.readthedocs.io)
- Tractoflow, version 2.2.1, [Tractoflow](https://github.com/scilus/tractoflow)

## Contact Information

For any questions or issues, please contact Obaï Bin Ka'b Ali at [ali.obaibk (at) gmail (dot) com].
