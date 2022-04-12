#load parameter input values
#load("combos.rda")

deg_dist="nbinom"
size=1.017340
mu=6.192894

alpha =  1

#sim_no = combos[id,1]

init_seed = runif(1,1,1000)

#Population parameters
population = 1000
covPattern = c(rep(0,population/2), rep(1,population/2)) #only use for Mixing

#Epidemic parameters - for 3 month timestep
num_init_infected = 1
beta_a = 1-((1-(0.01505*alpha))^3)
beta_l = 1-((1-(0.008911817*alpha))^3)
gamma_a = 3/2
gamma_l = 3/(36+4.2011+3.8709-2)

#Network model parameters
#Network_stats = list(c("Mixing"))
Network_stats = list(c("Degree"))
Prob_Distr = list(c("Multinomial_Poisson"))
Prob_Distr_Params = vector("list", 2)
Prob_Distr_Params[[1]] = c(1000) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
#Prob_Distr_Params[[2]] = c(0.1, 0.3, 0.6) 
Prob_Distr_Params[[2]] =  dnbinom(c(0:(population-1)),size=1.017340,mu=6.192894,log=FALSE)
  #c(c(0.0001, 0.0001, 0.2, 0.3, 0.2, 0.1, 0.05, 0.025, 0.025), rep(0.0001, 491))

#Hyper prior parameters
Prob_Distr_Params_hyperprior = vector("list", 3)
Prob_Distr_Params_hyperprior[[1]][1] = 'Dirichlet_Gamma'
Prob_Distr_Params_hyperprior[[2]][1] = 1 #gamma k (shape)
Prob_Distr_Params_hyperprior[[2]][2] = 1000 #gamma theta (scale)
#Prob_Distr_Params_hyperprior[[3]] = rep(1/length(Prob_Distr_Params[[2]]), length(Prob_Distr_Params[[2]])) #dirichlet alpha
Prob_Distr_Params_hyperprior[[3]] = rep(1/population, population) #dirichlet alpha

strong_prior = TRUE
#num_samples = combos[id,1]
num_samples = 25

#Genetic parameters
genetic_bits = 1024

#Bayes Inf MCMC
#n_mcmc_trials = 5000
n_mcmc_trials = 1000

#Prior vs. MCMC weight
#MCMC_wgt = combos[id,3]
MCMC_wgt = 0.75

#set transmission network to truth
init_P_truth_bool = FALSE

#set contact network
init_G_truth_bool = 0 #0 random, 1 P_start, 2 G_truth

#reporting bias
reporting_bias = TRUE
#ratio_orig_mean_deg = combos[id,4] #1 means the same, .5 is half the mean degree reported, 2 twice the mean degree reported
ratio_orig_mean_deg = 1
