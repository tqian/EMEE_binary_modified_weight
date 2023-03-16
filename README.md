# EMEE_binary_modified_weight
For paper with Yihan Bao 2021


This repository contains the code to replicate all numerical results in our paper: [title][link]

## Scripts
### folder: estimator implementation
* [WCLS_modified.R](estimator_implementation/WCLS_modified.R): the modified-EMEE estimator proposed in paper. (Section 4 of paper)

* [WCLS_original.R](estimator_implementation/WCLS_original.R): Original EMEE estimator by [Qian, et al](https://arxiv.org/abs/1906.00528). 

* [GEE_estimators.R](estimator_implementation/GEE_estimators.R): Implementation of GEE estimators, used for comparison with modified-EMEE.
    * log_linear_GEE: fitting GEE for MRT.
    * log_linear_GEE_geepack: fitting GEE using [geepack](https://www.jstatsoft.org/article/view/v015i02) package.
    * brm_DR, brm_MLE: brm estimator developed by [Richardson et al](https://arxiv.org/abs/1510.02430). This function is NOT used in the simulations presented in the paper.

### folder: simulation
* [simulation_results(marginal).R](simulations/simulation_results(marginal).R): Code for comparing the  
performance of pd-EMEE, EMEE and GEE for fully marginal effect, setting delta = 3 and = 10 respectively.

* [simulation_results(with_moderator).R](simulations/simulation_results(with_moderator).R): Code for comparing the  
performance of pd-EMEE, EMEE and GEE for causal effect with moderator, setting delta = 3 and = 10 respectively.

* [all_simulation_plots.R](simulations/all_simulation_plots.R): Code for reproducing the efficiency plot vs different Deltas and vs different randomization probability of treatment (and vs different reference regime) in the paper (and appendix H). 

* [plot_delta.R](simulations/plot_delta.R): Code for reproducing the raw simulation result of relative efficiency across different deltas, preparation for Figure 2 (left).

* [plot_proba.R](simulations/plot_proba.R): Code for reproducing the raw simulation result of relative efficiency across different deltas, preparation for Figure 2 (right). 

* [plot_K.R](simulations/plot_K.R): Code for reproducing the raw simulation result of relative efficiency across different reference regime, preparation for Appendix H.

* [utils_forK.R](simulations/utils_forK.R): Code for helper functions of plot_K.R.


* [dgm_simulation.R](simulations/dgm_simulation.R): a data-generative model for MRT data, used for simulation presented in the paper. (Described in the beginning of Section 4 of paper)

### folder: data analysis
* [analysis.R](data_analysis/analysis.R): Code for Drinkless data analysis of treatment in paper.

## Replicating Results

### Simulations
* Run [all_simulation_plots.R](simulations/all_simulation_plots.R) to reproduce the efficiency plot vs different Deltas and vs different randomization probability of treatment in the paper (Figure 2), and Appendix H. 

* Run [simulation_results(marginal).R](simulations/simulation_results(marginal).R) to reproduce the table of
estimating performance across pd-EMEE, EMEE and GEE for fully marginal effect, setting delta = 3 and 10 for generating Table 2 respectively.

* Run [simulation_results(with_moderator).R](simulations/simulation_results(with_moderator).R) to reproduce the table of estimating performance across pd-EMEE, EMEE and GEE for causal effect with moderator, setting delta = 3 and 10 for generating Table 3 respectively.

### Applications
* Run [analysis.R](data_analysis/analysis.R) to reproduce the data analysis section of paper.
