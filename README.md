# Code repo for *Characteristics of dynamic thickness change across diverse outlet glacier geometries and basal conditions*

Author: Donglai Yang

Contact: donglaiy@buffalo.edu or dyang379@gatech.edu

## Overview
This repository contains scripts to create [ISSM](https://issm.jpl.nasa.gov) model that generate the testbeds (Fig. 1 in the paper) and run the experiments (Fig. 2 in the paper). This repo does not contain the simulation data. The simulation data will be made available prior to publication.

## Getting started
After cloning to your local repository, add `functions` to your MATLAB path. Parts of the post-analysis scripts in `analysis_scripts` require [Climate Data Toolbox](https://chadagreene.com/CDT/CDT_Getting_Started.html).

### Add simulation data
Once simulation data is made available, you can put all testbed output (folders separated by each testbed) in `flatMISMIP_testbeds/long_models_yang`.

### Folders descriptions

`analysis_scripts`: it stores all the post-processing scripts.

`exp_files`: ISSM model domain coordinates (for flat bed models, to be differentiated from the rough fractal bed). No need to change.

`functions`: user-defined functions essential to both model running and post-processing steps.

`long_par_files`: parameterization files for ISSM models (flat bed models, to be differentiated from the rough fractal bed). If you intend to change certain experiment parameters, change `long_par_files/MISMIP_template.par` and then run `generate_pars.m` (one directory up). It will update all parameterization files.

`plots`: it stores images from post-processing analysis.

`wavybed_par_file`: parameterization files for ISSM models (rough fractal bed).

`wavybed_exp_files`: ISSM model model domain coordinates (rough fractal bed).

`long_models_yang`: the folder to store simulation output. Currently empty.

### Files descriptions

`runme_free_front_long.m`: main file, containing all experiments for all testbeds.

`runme_wavybed.m`: file for running rough bed experiment (Fig. 2D in the paper).

`rumme_param.csv`: parameter files essential for model running (i.e., `runme_free_front_long.m`).

`downsize_all_models.m`: downsize model output by removing non-essential fields, particularly `md.results.TransientSolution` from ISSM model.

`generate_model_combinations.m`: this file generates the 18 testbeds specs stored in `md_var_combinations.csv` which is used in simulations.

> Warning: in this script, as well as in many output name nomenclature, one may see fjord width written in "W5000_...","W8000_...", "W11000_...", these correspond to width of 4000m, 6000m, and 8000m in the paper.

`generate_random_bed.m`: generate the rough fractal bed with [artifical surface generator](https://www.mathworks.com/matlabcentral/fileexchange/60817-surface-generator-artificial-randomly-rough-surfaces). 

`generate_pars.m`: this generates the .par files in `long_par_files`.

`generate_wavy_pars.m`: this generates the .par files in `wavy_par_files`.


## Run analysis scripts
---
To run analysis scripts
