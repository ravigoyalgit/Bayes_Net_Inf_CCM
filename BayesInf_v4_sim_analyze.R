
library('tidyverse')
memory.limit(size = 12000)
################################
################################
#Degree
################################
################################

#directory_list = c("0.25", "0.5", "0.75", "1") #"0", "025", 
#directory_list = c("0.5", "0.75", "1") #"0", "025", 
directory_list = c("0.75") #"0", "025", 
#directory_list = c("0.5", "0.6", "0.7", "0.8", "0.9") #"0", "025", 

lambda = 50
range_start = 1001
range_end = 2000
mse_sim.df = data.frame(distr = NULL,
                     mse = NULL,
                     mcmc_wgt_df = NULL,
                     sample_size_df = NULL,
                     index = NULL,
                     num_sim = NULL)

directory = "0.75"
sample_size = 50

for (directory in directory_list) {
  
  for (sample_size in c(25, 50, 100, 200)) {
    
    #load_file = paste("C:/Users/ravij/OneDrive/Desktop/Network Research/network inference/NetBayes_git/NetBayes_git/Simulations3-init_P_truth_bool-FALSE/", directory, "/combined-results-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/COVIDSimulations/MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Downloads/v4/beta_l-0.0097_lambda-10_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v5-partial/beta_l-0.0097_lambda-20_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results_bias/v1/beta_l-0.0097_lambda-10_MCMC_wgt-", directory, "_num_samples-", sample_size, "_ratio_orig_mean_deg-0.5.rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results_bias/v1/beta_l-0.0097_lambda-10_MCMC_wgt-0.75_num_samples-", sample_size, "_ratio_orig_mean_deg-", directory, ".rda", sep = "") 
    #load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7/beta_l-0.0097_lambda-", lambda, "_MCMC_wgt-0.75_num_samples-", sample_size, ".rda", sep = "") 
    load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/beta_l-0.0097_lambda-", lambda, "_MCMC_wgt-0.75_num_samples-", sample_size, ".rda", sep = "") 
    
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
      
      if (((G_stats.df[c(range_start:range_end),lambda] == 0) %>% sum()) == 0) {
        
        for (p_counter in c(1:(range_end - range_start + 1))) {
          x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
          ProbDistr_stats.df = rbind(ProbDistr_stats.df, population*x2)
        }
        
        
        net_results_full.df = data.frame(mean_val = c(apply(G_stats.df[c(range_start:range_end),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                                         var_val = c(apply(G_stats.df[c(range_start:range_end),], 2, var), apply(ProbDistr_stats.df, 2, var)),
                                         stat_type = c(rep("Posterior", population), rep("Prior", population)),
                                         net_prop = rep(paste("Degree", c(0:(population-1)), sep="_"), 2),
                                         truth = rep(as.numeric(G_stats_truth), 2)
        ) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
          mutate(bias_val = abs(mean_val - truth)) %>%
          mutate(mse_val = bias2_val + var_val)
        
        net_results.df = net_results_full.df 
        
        prior_mse = net_results.df %>% filter(stat_type == "Prior") %>% pull(mse_val) %>% sum()
        post_mse = net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum()
        
        mse_sim_ind.df = data.frame(distr = c("prior", "post"),
                                    mse = c(prior_mse, post_mse),
                                    mcmc_wgt_df = rep(directory, 2),
                                    sample_size_df = rep(sample_size, 2),
                                    index = rep(i,2),
                                    num_sim = rep(length(results), 2)
                                    )
        mse_sim.df = bind_rows(mse_sim.df, mse_sim_ind.df)
          
      }
      print(paste("Sample_size: ", sample_size, "Index: ", i))
    }

    #print(paste("Sample_size: ", sample_size))

  }
  print(directory)
}

ggplot(data=mse_all, aes(x=as.factor(sample_size_df), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~mcmc_wgt_df, nrow = 2)


p1 = ggplot(data=mse_all %>% filter(mcmc_wgt_df == 0.75), aes(x=as.factor(sample_size_df), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Sexual Behavior Survey Samples", y = "Mean Squared Error", color = "Data Sources") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_fill_discrete(name = "Data Sources", labels = c("Multiple Sources\nIntegrated", "Sexual Behavior Survey\nOnly"))

svg(filename = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/fig_lambda_", lambda, ".svg", sep = ""),
    width = 20, height = 20)
p1
dev.off()

mse_all %>% filter(mcmc_wgt_df == 0.75) %>%
  pivot_wider(names_from = distr, values_from = mse) %>%
  mutate(MSE_improve = (prior-post)/prior * 100)

save.image(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/workspace_lambda", lambda, ".RData", sep =""))
saveRDS(mse_sim.df , file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))

############

lambda = 10
mse_sim_10.df = readRDS(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))
mse_sim_10.df$lambda = lambda

lambda = 20
mse_sim_20.df = readRDS(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))
mse_sim_20.df$lambda = lambda

lambda = 30
mse_sim_30.df = readRDS(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))
mse_sim_30.df$lambda = lambda

lambda = 40
mse_sim_40.df = readRDS(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))
mse_sim_40.df$lambda = lambda

lambda = 50
mse_sim_50.df = readRDS(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))
mse_sim_50.df$lambda = lambda

mse_sim.df = 
  bind_rows(
    bind_rows(
      bind_rows(
        bind_rows(mse_sim_10.df, mse_sim_20.df), 
      mse_sim_30.df),
    mse_sim_40.df),
  mse_sim_50.df) %>% as_tibble()

mse_all = mse_sim.df %>%
  group_by(distr, sample_size_df, mcmc_wgt_df, lambda) %>%
  summarise(mse = median(mse)) %>%
  ungroup()

p = ggplot(data=mse_all, aes(x=as.factor(sample_size_df), y=mse, fill=distr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Sexual Behavior Survey Samples", y = "Mean Squared Error", color = "Data Sources") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_fill_discrete(name = "Data Sources", labels = c("Multiple Sources\nIntegrated", "Sexual Behavior Survey\nOnly")) +
  facet_wrap(~lambda, nrow = 2)

svg(filename = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/fig_mse_all.svg", sep = ""),
    width = 20, height = 15)
p
dev.off()

############


lambda = 20
mse_all.df = readRDS(paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/mse_lambda", lambda, ".rds", sep =""))

sample_size = 50
median.df = mse_all.df %>% 
  filter(sample_size_df == sample_size,
         distr == "post") %>%
  arrange(mse)

i_index = median.df[round(nrow(median.df)/2),"index"]

#load_file = paste("C:/Users/ravij/Downloads/v4/beta_l-0.0097_lambda-20_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/beta_l-0.0097_lambda-", lambda, "_MCMC_wgt-0.75_num_samples-", sample_size, ".rda", sep = "") 
load(load_file)

G_stats.df = results[[i_index]][[1]]
Prob_Distr_Params = results[[i_index]][[2]]
Prob_Distr_Params_hyperprior = results[[i_index]][[3]]
n_mcmc_trials = results[[i_index]][[4]]
G_stats_truth = results[[i_index]][[5]]

############

G_stats_lim.df = as_tibble(G_stats.df[c(1:2000), c((lambda-5):(6+lambda))])
colnames(G_stats_lim.df) =  paste("Degree", c((lambda-5):(6+lambda)), sep = " ")
G_stats_lim.df$iter = c(1:nrow(G_stats_lim.df))
G_stats_long.df = pivot_longer(G_stats_lim.df,
             cols = starts_with("Degree"),
             names_to = "Prop",
             values_to = "Count")

p = ggplot(data=G_stats_long.df, aes(x=iter, y=Count, fill=Prop)) +
  geom_point() +
  labs(x="MCMC Iteration", y = "Number of Nodes") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  geom_vline(xintercept = 1000, col = "red") +
  facet_wrap(vars(Prop), scales = "free_y", nrow = 4) 


svg(filename = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v7x250/fig_trace_lambda_", lambda, ".svg", sep = ""),
    width = 20, height = 20)
p
dev.off()

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
    load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6/beta_l-0.0097_lambda-50_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
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
      
      net_results.df = data.frame(mean_val = c(apply(G_stats.df[c(4001:5000),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                                  var_val = c(apply(G_stats.df[c(4001:5000),], 2, var), apply(ProbDistr_stats.df, 2, var)),
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

#############
##Trace plots

directory = 0.75
sample_size = 50
#load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6-vanilla/beta_l-0.0097_lambda-50_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
load_file = paste("C:/Users/ravij/Dropbox/Academic/Research/Projects/Bayes_Net_Inf_CCM/Bayes_Net_Inf_CCM_results/v6/beta_l-0.0097_lambda-50_MCMC_wgt-", directory, "_num_samples-", sample_size, ".rda", sep = "") 
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

par(mfrow = c(3,5))

for (i in c(1:15)) {
  
  if (i == 1) {
    y_axis_lab = "Mixing 1-1"
  }
  if (i == 2) {
    y_axis_lab = "Mixing 1-2"
  }
  if (i == 3) {
    y_axis_lab = "Mixing 2-2"
  }
  
  plot(G_stats.df[c(1:2000),i], ylab = y_axis_lab, xlab = "MCMC Iteration",
       ylim = c(min(c(G_stats.df[c(1:2000),i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[1]))),
                max(c(G_stats.df[c(1:2000),i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[1])))))
  abline(v=1000, col = 'red')
  
  #abline(h=G_stats_truth[i], col = 'red')
  #abline(h=mean(G_stats.df[c(4001:5000),i]), col = 'blue')
  #abline(h=mean(ProbDistr_stats.df[,i]), col = 'green')
}

par(mfrow = c(1,1))
