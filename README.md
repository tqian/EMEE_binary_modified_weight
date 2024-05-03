# EMEE_binary_modified_weight
For paper with Yihan Bao 2021


This repository contains the code to replicate all numerical results in our paper: [title][link]

## Scripts
### folder: estimator implementation
* [WCLS_modified.R](estimator_implementation/WCLS_modified.R): the pd-EMEE estimator proposed in the paper. (Section 3 of paper)

* [WCLS_original.R](estimator_implementation/WCLS_original.R): Original EMEE estimator by [Qian, et al](https://arxiv.org/abs/1906.00528). 

* [GEE_estimators.R](estimator_implementation/GEE_estimators.R): Implementation of GEE estimators, used for comparison with modified-EMEE.
    * log_linear_GEE: fitting GEE for MRT.
    * log_linear_GEE_geepack: fitting GEE using [geepack](https://www.jstatsoft.org/article/view/v015i02) package.
    * brm_DR, brm_MLE: brm estimator developed by [Richardson et al](https://arxiv.org/abs/1510.02430). This function is NOT used in the simulations presented in the paper.

 * [EMEE_ImprovedEffByProjection.R](estimator_implementation/EMEE_ImprovedEffByProjection.R): the pd-EMEE2 estimator proposed in the paper. (Section 4 of paper)
 
### folder: simulation
* [simulation_results(marginal).R](simulations/simulation_results(marginal).R): Code for comparing the  
performance of EMEE, pd-EMEE and pd-EMEE2 for fully marginal effect, setting delta = 3 and = 10 respectively.

* [simulation_results(with_moderator).R](simulations/simulation_results(with_moderator).R): Code for comparing the  
performance of EMEE, pd-EMEE and pd-EMEE2 for causal effect with moderator, setting delta = 3 and = 10 respectively.

* [simulation_results(marginal).R](simulations/simulation_results(marginal)_GEE.R): Code for comparing the  
performance of EMEE, pd-EMEE and GEE for fully marginal effect, setting delta = 3 and = 10 respectively.

* [simulation_results(with_moderator).R](simulations/simulation_results(with_moderator)_GEE.R): Code for comparing the  
performance of EMEE, pd-EMEE and GEE for causal effect with moderator, setting delta = 3 and = 10 respectively.

* [simulation_plots.R](simulations/simulation_plots.R): Code for reproducing the efficiency plot vs different Deltas and vs different randomization probability of treatment (and vs different reference regime) in the paper (and appendix J). 

* [plot_delta.R](simulations/plot_delta.R): Code for reproducing the raw simulation result of relative efficiency across different deltas, preparation for Figure 2 (left).

* [plot_proba.R](simulations/plot_proba.R): Code for reproducing the raw simulation result of relative efficiency across different deltas, preparation for Figure 2 (right). 

* [plot_K.R](simulations/plot_K.R): Code for reproducing the raw simulation result of relative efficiency across different reference regimes, preparation for Appendix J.

* [utils_forK.R](simulations/utils_forK.R): Code for helper functions of plot_K.R.

* [dgm_simulation.R](simulations/dgm_simulation.R): a data-generative model for MRT data, used for simulation presented in the paper. (Described in Appendix H)

### folder: HeartSteps_analysis
* [EMEE_userinput_ipw.R](HeartSteps_analysis/EMEE_userinput_ipw.R): Code for the generalized version of pd-EMEE in the paper.

* [HeartSteps - 1. preprocessing, analysis.R](HeartSteps_analysis/HeartSteps - 1. preprocessing, analysis.R): Code for preprocessing and analyzing HeartSteps dataset in paper.

* [HeartSteps - 2. plotting.R](HeartSteps_analysis/HeartSteps - 2. plotting.R): Code for plotting HeartSteps analysis in paper.

### folder: DrinkLess_analysis
* [DrinkLess - 1. plotting.R](DrinkLess_analysis/DrinkLess - 1.  preprocessing, analysis.R):  Code for preprocessing and analyzing DrinkLess dataset in the paper.

* [DrinkLess - 2. plotting.R](DrinkLess_analysis/DrinkLess - 2. plotting.R): Code for plotting Drinkless data analysis in the paper.


## Replicating Results

### Simulations
* Run [simulation_plots.R](simulations/simulation_plots.R) to reproduce the efficiency plot vs different Deltas and vs different randomization probabilities of treatment in the paper (Figure 2). 

* Run [simulation_results(marginal).R](simulations/simulation_results(marginal).R) to reproduce the table of
estimating performance across EMEE, pd-EMEE and pd-EMEE2 for fully marginal effect, setting delta = 3 and 10 for generating Table 1 respectively.

* Run [simulation_results(with_moderator).R](simulations/simulation_results(with_moderator).R) to reproduce the table of estimating performance across pd-EMEE, EMEE and pd-EMEE2 for causal effect with moderator, setting delta = 3 and 10 for generating Table 2 respectively.

* Run [simulation_results(marginal)_GEE.R](simulations/simulation_results(marginal).R) to reproduce the table of
estimating performance across EMEE, pd-EMEE and GEE for fully marginal effect, setting delta = 3 and 10 for generating Table 1 respectively.

* Run [simulation_results(with_moderator)_GEE.R](simulations/simulation_results(with_moderator).R) to reproduce the table of estimating performance across pd-EMEE, EMEE and GEE for causal effect with moderator, setting delta = 3 and 10 for generating Table 2 respectively.

### Applications
* Run codes in HeartSteps_analysis and DrinkLess_analysis in sequence to generate Figure 3a and 3b respectively.
