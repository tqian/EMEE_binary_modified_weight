# EMEE_binary_modified_weight
For paper with Yihan Bao 2021


This repository contains the code to replicate all numerical results in our paper: [title][link]

## Scripts
* [dgm_simulation.R](dgm_simulation.R): a data-generative model for MRT data, used for simulation presented in the paper. (Described in the beginning of Section 5 of paper)

* [WCLS_modified.R](WCLS_modified.R): the modified-EMEE estimator proposed in paper. (Section 4 of paper)

* [WCLS_original.R](WCLS_original.R): Original EMEE estimator by [Qian, et al](https://arxiv.org/abs/1906.00528). 

* [GEE_estimators.R](GEE_estimators.R): Implementation of GEE estimators, used for comparison with modified-EMEE.
    * log_linear_GEE: fitting GEE for MRT.
    * log_linear_GEE_geepack: fitting GEE using [geepack](https://www.jstatsoft.org/article/view/v015i02) package.
    * brm_DR, brm_MLE: brm estimator developed by [Richardson et al](https://arxiv.org/abs/1510.02430). This function is NOT used in the simulations presented in the paper.

* [simulations_R](simulations_R): Code for reproducing the simulation results on consistency and efficiency tables of the paper. 

* [plot_delta.R](plot_delta.R): Code for reproducing efficiency plot vs different Deltas in paper.

* [plot_proba.R](plot_proba.R): Code for reproducing efficiency plot vs different randomization probability of treatment in paper.

## Replicating Results

### Simulations
* Run [simulations.R](simulations.R) to reproduce the simulation results on consistency and efficiency tables of the paper. 

* Run [plot_delta.R](plot_delta.R) to reproduce the efficiency plot vs different Deltas of the paper. 

* Run [plot_proba.R](plot_proba.R) to reproduce the efficiency plot vs different randomization probability of treatment of the paper. 

### Applications
