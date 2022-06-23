#' @title read_params
#' @description sets parameters for simulating network data and performing Bayesian inference, which is run with Bayes_Net_Inf_CCM()
#' @param population number of individuals in the network
#' @param beta_a probability of infection for acute phase
#' @param beta_l probability of infection for long-term phase
#' @param gamma_a length of infection for acute phase
#' @param gamma_l length of infection for long-term phase
#' @param num_samples number of behavior survey samples
#' @param genetic_bits length of genetic sequence
#' @param n_mcmc_trials number of MCMC trials to perform
#' @param MCMC_wgt weight from 0 to 1 indicating how much weight to put on the prior versus the MCMC weight
#' @param init_P_truth_bool Boolean, whether to start at the true transmission network; TRUE or FALSE
#' @param init_G_truth_bool where to start at the true contact network; 0 is random, 1 is P_start, 2 is G_truth
#' @param reporting_bias Boolean, whether to simulate reporting bias; TRUE or FALSE
#' @param ratio_orig_mean_deg ratio of how much to inflate or deflate number of partners reported; is no increase/decrease, 0.5 is half the mean degree reported, 2 is twice the mean degree reported
#' @return none
#' @export

read_params <- function(population = 1000,
                        beta_a = 0,
                        beta_l = (1-(1-0.128)^(1/14)),
                        gamma_a = 1/4,
                        gamma_l = 1/14,
                        num_samples = 25,
                        genetic_bits = 1024,
                        n_mcmc_trials = 1000,
                        MCMC_wgt = 0,
                        init_P_truth_bool = FALSE,
                        init_G_truth_bool = 0,
                        reporting_bias = FALSE,
                        ratio_orig_mean_deg = 1) {

init_seed = runif(1,1,1000)

#Population parameters
population = population
covPattern = c(rep(0,population/2), rep(1,population/2)) #only use for Mixing

#Epidemic parameters
num_init_infected = 1
beta_a = beta_a
beta_l = beta_l
gamma_a = gamma_a
gamma_l = gamma_l

#Network model parameters
#Network_stats = list(c("Mixing"))
Network_stats = list(c("Degree"))
Prob_Distr = list(c("Multinomial_Poisson"))
Prob_Distr_Params = vector("list", 2)
Prob_Distr_Params[[1]] = c(1000) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
#Prob_Distr_Params[[2]] = c(0.1, 0.3, 0.6)
Prob_Distr_Params[[2]] = dpois(c(0:(population-1)), lambda=10, log = FALSE)
  #c(c(0.0001, 0.0001, 0.2, 0.3, 0.2, 0.1, 0.05, 0.025, 0.025), rep(0.0001, 491))


#Hyper prior parameters
Prob_Distr_Params_hyperprior = vector("list", 3)
Prob_Distr_Params_hyperprior[[1]][1] = 'Dirichlet_Gamma'
Prob_Distr_Params_hyperprior[[2]][1] = 1 #gamma k (shape)
Prob_Distr_Params_hyperprior[[2]][2] = 1000 #gamma theta (scale)
#Prob_Distr_Params_hyperprior[[3]] = rep(1/length(Prob_Distr_Params[[2]]), length(Prob_Distr_Params[[2]])) #dirichlet alpha
Prob_Distr_Params_hyperprior[[3]] = rep(1/population, population) #dirichlet alpha


strong_prior = TRUE
num_samples = num_samples

#Genetic parameters
genetic_bits = genetic_bits

#Bayes Inf MCMC
n_mcmc_trials = n_mcmc_trials

#Prior vs. MCMC weight
MCMC_wgt = MCMC_wgt

#set transmission network to truth
init_P_truth_bool = init_P_truth_bool

#set contact network
init_G_truth_bool = init_G_truth_bool #0 random, 1 P_start, 2 G_truth

#Reporting bias
reporting_bias = reporting_bias
ratio_orig_mean_deg = ratio_orig_mean_deg #1 means the same, .5 is half the mean degree reported, 2 twice the mean degree reported

}
