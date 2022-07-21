#' @title Bayes_Net_Inf_CCM
#' @description generates data and estimates contact network properties using survey, clinical, and genetic data
#' @description sets parameters for simulating network data and performing Bayesian inference, which is run with Bayes_Net_Inf_CCM()
#' @param population number of individuals in the network
#' @param beta_a probability of infection for acute phase
#' @param beta_l probability of infection for long-term phase
#' @param gamma_a length of infection for acute phase
#' @param gamma_l length of infection for long-term phase
#' @param degree_dist degree distribution of contacts, can be "poisson" (requires lambda parameter) or "nbinom" (requires size and mu parameters)
#' @param lambda mean degree for the Poisson degree distribution
#' @param size scale for the Negative Binomial degree distribution
#' @param mu mean degree for the Negative Binomial degree distribution
#' @param num_samples number of behavior survey samples
#' @param genetic_bits length of genetic sequence
#' @param n_mcmc_trials number of MCMC trials to perform
#' @param MCMC_wgt weight from 0 to 1 indicating how much weight to put on the prior versus the MCMC weight
#' @param ratio_orig_mean_deg ratio of how much to inflate or deflate number of partners reported; is no increase/decrease, 0.5 is half the mean degree reported, 2 is twice the mean degree reported
#' @param epidemic select a set of epidemic parameters to use instead of entering custom parameters, can be set to "HIV" or "COVID" (overrides and parameters set in the function's arguments); set to "NULL" to use custom parameters via the function's arguments
#' @return results_list = list(G_stats.df, Prob_Distr_Params, Prob_Distr_Params_hyperprior, n_mcmc_trials, G_stats_truth, init_seed), results of data simulation and Bayesian inference, including estimated contact network statistics
#' @export

Bayes_Net_Inf_CCM <- function(population = 1000,
                              beta_a = 0,
                              beta_l = (1-(1-0.128)^(1/14)),
                              gamma_a = 1/4,
                              gamma_l = 1/14,
                              lambda = 10,
                              size = NULL,
                              mu = NULL,
                              degree_dist = "poisson",
                              num_samples = 25,
                              genetic_bits = 1024,
                              n_mcmc_trials = 2000,
                              MCMC_wgt = 0.75,
                              ratio_orig_mean_deg = 1,
                              epidemic = NULL) {

  ################## Setup ##################

  CCMnetpy::CCMnet_python_setup()

  ################## Parameters ##################

  #HIV epidemic parameters
  if(epidemic=="COVID") {
    population = 1000
    beta_a = 0
    beta_l = (1-(1-0.128)^(1/14))
    gamma_a = 1/4
    gamma_l = 1/14
    lambda = 20 #
    degree_dist = "poisson"
    num_samples = 100 #
    genetic_bits = 1024
    n_mcmc_trials = 2000
    MCMC_wgt = 0.75
    ratio_orig_mean_deg = 1
  }

  #COVID epidemic parameters
  if(epidemic=="HIV") {
    population = 1000
    beta_a = 1-((1-(0.01505))^3)
    beta_l = 1-((1-(0.008911817))^3)
    gamma_a = 3/2
    gamma_l = 3/(36+4.2011+3.8709-2)
    size = 1.017340
    mu = 6.192894
    degree_dist = "nbinom"
    num_samples = 100 #
    genetic_bits = 1024
    n_mcmc_trials = 2000
    MCMC_wgt = 0.75
    ratio_orig_mean_deg = 1
  }

  init_seed = runif(1,1,1000)

  #population parameters
  population = population
  covPattern = c(rep(0,population/2), rep(1,population/2)) #only use for Mixing

  #epidemic parameters
  num_init_infected = 1
  beta_a = beta_a
  beta_l = beta_l
  gamma_a = gamma_a
  gamma_l = gamma_l

  #network model parameters
  #Network_stats = list(c("Mixing"))
  Network_stats = list(c("Degree"))
  Prob_Distr = list(c("Multinomial_Poisson"))
  Prob_Distr_Params = vector("list", 2)
  Prob_Distr_Params[[1]] = c(1000) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
  #Prob_Distr_Params[[2]] = c(0.1, 0.3, 0.6)
  if (degree_dist == "poisson") {Prob_Distr_Params[[2]] = dpois(c(0:(population-1)), lambda=lambda, log=FALSE)}
  if (degree_dist == "nbinom") {Prob_Distr_Params[[2]] = dnbinom(c(0:(population-1)), size=size, mu=mu, log=FALSE)}
  #c(c(0.0001, 0.0001, 0.2, 0.3, 0.2, 0.1, 0.05, 0.025, 0.025), rep(0.0001, 491))

  #hyper prior parameters
  Prob_Distr_Params_hyperprior = vector("list", 3)
  Prob_Distr_Params_hyperprior[[1]][1] = 'Dirichlet_Gamma'
  Prob_Distr_Params_hyperprior[[2]][1] = 1 #gamma k (shape)
  Prob_Distr_Params_hyperprior[[2]][2] = 1000 #gamma theta (scale)
  #Prob_Distr_Params_hyperprior[[3]] = rep(1/length(Prob_Distr_Params[[2]]), length(Prob_Distr_Params[[2]])) #dirichlet alpha
  Prob_Distr_Params_hyperprior[[3]] = rep(1/population, population) #dirichlet alpha

  strong_prior = TRUE

  num_samples = num_samples

  #genetic parameters
  genetic_bits = genetic_bits

  #Bayes Inf MCMC
  n_mcmc_trials = n_mcmc_trials

  #Prior vs. MCMC weight
  MCMC_wgt = MCMC_wgt

  #set transmission network to truth
  init_P_truth_bool = FALSE

  #set contact network
  init_G_truth_bool = 1 #0 random, 1 P_start, 2 G_truth

  #reporting bias
  reporting_bias = TRUE
  ratio_orig_mean_deg = ratio_orig_mean_deg #1 means the same, .5 is half the mean degree reported, 2 twice the mean degree reported

  ################## Bayesian Inference ##################

  Init_nets = Inital_bayes_inf(population,Network_stats,
                               Prob_Distr, Prob_Distr_Params,covPattern,
                               Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist,
                               P_truth, init_P_truth_bool, G_truth, init_G_truth_bool)

  G = Init_nets[[1]] #G_truth
  P = Init_nets[[2]] #P_truth

  P_a = NULL
  G_a = NULL

  G_stats.df = c()
  ecount_P = c()

  for (mcmc_counter in c(1:n_mcmc_trials)) {

    G_info = Update_G(G,P,Ia,Il,R,beta_a,beta_l,Prob_Distr_Params=Prob_Distr_Params, Network_stats = Network_stats, Prob_Distr = Prob_Distr, covPattern = covPattern)
    G = G_info[[1]]
    G_stats = as.numeric(G_info[[2]])

    G_stats.df = rbind(G_stats.df, G_stats)

    P = Update_P(G,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)

    Prob_Distr_Params = Update_Prob_Distr_Params(G,Prob_Distr_Params_hyperprior=Prob_Distr_Params_hyperprior, Network_stats = Network_stats, Prob_Distr = Prob_Distr, Prob_Distr_Params = Prob_Distr_Params, G_stats = G_stats, MCMC_wgt)

    print(mcmc_counter)

  }

  ################## Results ##################

  results_list = list(G_stats.df,
                      Prob_Distr_Params,
                      Prob_Distr_Params_hyperprior,
                      n_mcmc_trials,
                      G_stats_truth,
                      init_seed)

}
