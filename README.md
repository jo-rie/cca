# Code for the paper "Decomposition and graphical correspondence analysis of checkerboard copulas"

This repository contains the code for the paper "Decomposition and graphical correspondence analysis of checkerboard copulas".
The data is generated in Matlab using the scripts `simulation_runs.m` and `matin_turbine_data_analysis.m`.
They use the [ReadYAML package](https://github.com/llerussell/ReadYAML).

The figures are generated in Python using the scripts `run_plots.py` and `run_plots_turbine.py`.
Further Python source code is in the `python_src`, for Matlab in `matlab_src`.

The Decomposition is organized into classes, the corresponding files start with `Copula`. 
The other names should be self-explanatory.

In `config_local.yml`, the path to the [turbine data](https://doi.org/10.1111/rssc.12421) should be specified for the corresponding analysis.