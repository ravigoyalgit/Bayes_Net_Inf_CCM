# COVID/Results-v1
Run with base set of simulation scripts (Bayes_Net_Inf_CCM)
Only change was to BayesInf_v4_func_read_params.R
 - beta_l = 5.5/49.99153*14 *0.15
 - gamma_a = 1/4
 - gamma_l = 1/14
 - Prob_Distr_Params[[2]] = ppois(..., lambda=49.99153, ...)

Degree
MCMC_wgt = 0, 0.25, 0.5, 0.75, 1
num_samples = 25, 50, 100, 200
init_P_truth_bool = FALSE
