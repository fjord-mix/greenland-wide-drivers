# model-runs

Repository with scripts for setting up, running, and plotting the model outputs in Mas e Braga et al. (XXXX). 

The model itself can be downloaded [here](https://github.com/fjord-mix/box-model), and the release version published can be found [(not yet uploaded)](https://github.com/fjord-mix/box-model)

## Required datasets

- Supporting datasets published together with the manuscript
- OMG datasets (using the PO.DAAC API)
    - [CTD](https://doi.org/10.5067/OMGEV-CTDS1)
    - [AXCTD](https://doi.org/10.5067/OMGEV-AXCT1)
- Subglacial discharge from [Karlsson et al. (2023)](https://doi.org/10.22008/FK2/BOVBVR)

## Basic walkthrough
Driver files for running the model ensembles:
- `driver_rpm_grl_fjords.m` is the file to run for the large ensemble presented in the manuscript
    - all accessory functions and plotting should be easy to find following this script
- `driver_assess_misfits.m` is the file to run for evaluating whether the parameter space covered was enough to find the misfit minima
    - these can be plotted from `plotting_scripts/plot_misfits_per_parameter.m`
- all plotting functions are in the homonimous folder
- the folder `z_legacy` contains previous iterations and code that might be useful depending on what the user wants to do

