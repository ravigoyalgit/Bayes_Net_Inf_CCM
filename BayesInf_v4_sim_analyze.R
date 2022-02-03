
library('tidyverse')

################################
################################
#Degree
################################
################################

#prior_mse_all = c()
#post_mse_all = c()

directory_list = c("0.25", "0.5", "0.75", "1") #"0", "025", 
directory_list = c("0.5", "0.75", "1") #"0", "025", 


mse_all = data.frame(distr = NULL,
                     mse = NULL,
                     mcmc_wgt = NULL,
                     sample_size = NULL)

for (directory in directory_list) {
  
  for (sample_size in c(25, 50, 100, 200)) {
    
    #load_file = paste("C:/Users/ravij/OneDrive/Desktop/Network Research/network inference/NetBayes_git/NetBayes_git/Simulations3-init_P_truth_bool-FALSE/", directory, "/combined-results-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/COVIDSimulations/MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    load_file = paste("C:/Users/ravij/Downloads/v4/beta_l-0.0097_lambda-10_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v5-partial/beta_l-0.0097_lambda-20_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    
    load(load_file)
    
    prior_mse = c()
    post_mse = c()
    
    population = 1000
    
    for (i in c(1:length(results))) {
      
      G_stats.df = results[[i]][[1]]
      Prob_Distr_Params = results[[i]][[2]]
      Prob_Distr_Params_hyperprior = results[[i]][[3]]
      n_mcmc_trials = results[[i]][[4]]
      G_stats_truth = results[[i]][[5]]
      
      ProbDistr_stats.df = c()
      
      for (p_counter in c(1:200)) {
        x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
        ProbDistr_stats.df = rbind(ProbDistr_stats.df, population*x2)
      }
      
      
      net_results_full.df = data.frame(mean_val = c(apply(G_stats.df[c(801:1000),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                                       var_val = c(apply(G_stats.df[c(801:1000),], 2, var), apply(ProbDistr_stats.df, 2, var)),
                                       stat_type = c(rep("Posterior", population), rep("Prior", population)),
                                       net_prop = rep(paste("Degree", c(0:(population-1)), sep="_"), 2),
                                       truth = rep(as.numeric(G_stats_truth), 2)
      ) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
        mutate(bias_val = abs(mean_val - truth)) %>%
        mutate(mse_val = bias2_val + var_val)
      
      net_results.df = net_results_full.df #%>% filter(net_prop %in% paste("Degree", c(2:7), sep="_"))
      
      prior_mse = c(prior_mse, net_results.df %>% filter(stat_type == "Prior") %>% pull(mse_val) %>% sum())
      post_mse = c(post_mse, net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum())
      
      #if ((net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum()) > 300000) {
      #  print(i)
      #  print(net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum())
      #}
    }
    
    #post_mse %>% hist()
    
    prior_0 = prior_mse %>% median()
    post_0 = post_mse %>% median()
    
    mse_ind = data.frame(distr = c("prior", "post"),
                         mse= c(prior_0, post_0),
                         mcmc_wgt = rep(directory, 2),
                         sample_size = rep(sample_size, 2))
    
    mse_all = bind_rows(mse_all, mse_ind)
    
    #prior_mse_all = c(prior_mse_all, prior_0)
    #post_mse_all = c(post_mse_all, post_0)
    
  }
  print(directory)
}



ggplot(data=mse_all, aes(x=as.factor(sample_size), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~mcmc_wgt, nrow = 2)


ggplot(data=mse_all %>% filter(mcmc_wgt == 0.75), aes(x=as.factor(sample_size), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Sexual Behavior Survey Samples", y = "Mean Squared Error", color = "Data Sources") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_fill_discrete(name = "Data Sources", labels = c("Multiple Sources\nIntegrated", "Sexual Behavior Survey\nOnly"))

mse_all %>% filter(mcmc_wgt == 0.75) %>%
  pivot_wider(names_from = distr, values_from = mse) %>%
  mutate(MSE_improve = (prior-post)/prior * 100)

############

directory = 0.5
sample_size = 25
load_file = paste("C:/Users/ravij/Downloads/v4/beta_l-0.0097_lambda-20_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
load(load_file)

G_stats.df = results[[i]][[1]]
Prob_Distr_Params = results[[i]][[2]]
Prob_Distr_Params_hyperprior = results[[i]][[3]]
n_mcmc_trials = results[[i]][[4]]
G_stats_truth = results[[i]][[5]]

############

par(mfrow = c(5,5))

for (i in c(15:39)) {

plot(G_stats.df[,i], ylab = paste("Degree", i, sep = " "), xlab = "Sample",
     ylim = c(min(c(G_stats.df[,i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[1]))),
              max(c(G_stats.df[,i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[1])))))
abline(h=G_stats_truth[i], col = 'red')
abline(h=mean(G_stats.df[c(801:1000),i]), col = 'blue')
abline(h=mean(ProbDistr_stats.df[,i]), col = 'green')
}

par(mfrow = c(1,1))




###############################
###############################
#MIXING
###############################
###############################


directory_list = c("0.25", "0.5", "0.75", "1") #"0", "025", 

mse_all = data.frame(distr = NULL,
                     mse = NULL,
                     mcmc_wgt = NULL,
                     sample_size = NULL)

for (directory in directory_list) {
  
  for (sample_size in c(25, 50, 100, 200)) {
    
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6/beta_l-0.0058_lambda-50_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6/beta_l-0.0097_lambda-10_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6-vanilla/beta_l-0.0097_lambda-50_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    
    load(load_file)
    
    prior_mse = c()
    post_mse = c()
    
    for (i in c(1:length(results))) {
      
      G_stats.df = results[[i]][[1]]
      Prob_Distr_Params = results[[i]][[2]]
      Prob_Distr_Params_hyperprior = results[[i]][[3]]
      n_mcmc_trials = results[[i]][[4]]
      G_stats_truth = results[[i]][[5]]
      
      ProbDistr_stats.df = c()
      
      for (p_counter in c(1:1000)) {
        x1 = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
        x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
        ProbDistr_stats.df = rbind(ProbDistr_stats.df, x1*x2)
      }
      
      net_results.df = data.frame(mean_val = c(apply(G_stats.df[-c(1:4000),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                                  var_val = c(apply(G_stats.df[-c(1:4000),], 2, var), apply(ProbDistr_stats.df, 2, var)),
                                  stat_type = c(rep("Posterior", 3), rep("Prior", 3)),
                                  net_prop = rep(c("00", "01", "11"), 2),
                                  truth = rep(as.numeric(G_stats_truth), 2)
      ) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
        mutate(bias_val = abs(mean_val - truth)) %>%
        mutate(mse_val = bias2_val + var_val)
      
      prior_mse = c(prior_mse, net_results.df %>% filter(stat_type == "Prior") %>% pull(mse_val) %>% sum())
      post_mse = c(post_mse, net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum())
      
    }
    
    prior_0 = prior_mse %>% median()
    post_0 = post_mse %>% median()
    
    mse_ind = data.frame(distr = c("prior", "post"),
                         mse= c(prior_0, post_0),
                         mcmc_wgt = rep(directory, 2),
                         sample_size = rep(sample_size, 2))
    
    mse_all = bind_rows(mse_all, mse_ind)
    
  }
  print(directory)
}

ggplot(data=mse_all, aes(x=as.factor(sample_size), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~mcmc_wgt)

ggplot(data=mse_all %>% filter(mcmc_wgt == 0.75), aes(x=as.factor(sample_size), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Sexual Behavior Survey Samples", y = "Mean Squared Error", color = "Data Sources") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_fill_discrete(name = "Data Sources", labels = c("Multiple Sources\nIntegrated", "Sexual Behavior Survey\nOnly"))

mse_all %>% filter(mcmc_wgt == 0.75) %>%
  pivot_wider(names_from = distr, values_from = mse) %>%
  mutate(MSE_improve = (prior-post)/prior * 100)

57#############
##Trace plots

directory = 0.75
sample_size = 50
load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6-vanilla/beta_l-0.0097_lambda-50_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
load(load_file)

load(load_file)

prior_mse = c()
post_mse = c()

for (i in c(1:length(results))) {
  
  G_stats.df = results[[i]][[1]]
  Prob_Distr_Params = results[[i]][[2]]
  Prob_Distr_Params_hyperprior = results[[i]][[3]]
  n_mcmc_trials = results[[i]][[4]]
  G_stats_truth = results[[i]][[5]]
  
  ProbDistr_stats.df = c()
  
  for (p_counter in c(1:1000)) {
    x1 = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
    x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
    ProbDistr_stats.df = rbind(ProbDistr_stats.df, x1*x2)
  }
  
  net_results.df = data.frame(mean_val = c(apply(G_stats.df[-c(1:4000),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                              var_val = c(apply(G_stats.df[-c(1:4000),], 2, var), apply(ProbDistr_stats.df, 2, var)),
                              stat_type = c(rep("Posterior", 3), rep("Prior", 3)),
                              net_prop = rep(c("00", "01", "11"), 2),
                              truth = rep(as.numeric(G_stats_truth), 2)
  ) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
    mutate(bias_val = abs(mean_val - truth)) %>%
    mutate(mse_val = bias2_val + var_val)
  
  prior_mse = c(prior_mse, net_results.df %>% filter(stat_type == "Prior") %>% pull(mse_val) %>% sum())
  post_mse = c(post_mse, net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum())
  
}


i = which(order(post_mse) == 25)

G_stats.df = results[[i]][[1]]
Prob_Distr_Params = results[[i]][[2]]
Prob_Distr_Params_hyperprior = results[[i]][[3]]
n_mcmc_trials = results[[i]][[4]]
G_stats_truth = results[[i]][[5]]

ProbDistr_stats.df = c()

for (p_counter in c(1:1000)) {
  x1 = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
  x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
  ProbDistr_stats.df = rbind(ProbDistr_stats.df, x1*x2)
}

############

par(mfrow = c(1,3))

for (i in c(1:3)) {
  
  if (i == 1) {
    y_axis_lab = "Mixing 1-1"
  }
  if (i == 2) {
    y_axis_lab = "Mixing 1-2"
  }
  if (i == 3) {
    y_axis_lab = "Mixing 2-2"
  }
  
  plot(G_stats.df[,i], ylab = y_axis_lab, xlab = "MCMC Iteration",
       ylim = c(min(c(G_stats.df[,i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[1]))),
                max(c(G_stats.df[,i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[1])))))
  abline(v=4000, col = 'red')
  
  #abline(h=G_stats_truth[i], col = 'red')
  #abline(h=mean(G_stats.df[c(4001:5000),i]), col = 'blue')
  #abline(h=mean(ProbDistr_stats.df[,i]), col = 'green')
}

par(mfrow = c(1,1))
