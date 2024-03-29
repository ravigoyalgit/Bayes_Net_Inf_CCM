generate_epidemic_data <- function(beta_a = 1/5,
                                   beta_l = 1/2.5,
                                   gamma_a = 1/10,
                                   gamma_l = 1/100,
                                   population = 500,
                                   Prob_Distr_Params = list(list(rep(1,  10))),
                                   Prob_Distr = 'DirMult',
                                   num_init_infected = 1,
                                   genetic_bits = 1024,
                                   Network_stats = NULL,
                                   strong_prior = NULL,
                                   num_samples = 100,
                                   covPattern = NULL) {
  
  invalid_epidemic = TRUE
  counter = 1
  max_try = 5
  
  while(invalid_epidemic & counter < max_try) {

    CCMnet_Result = CCMnetpy::CCMnet_constr(Network_stats=Network_stats,
                                  Prob_Distr=Prob_Distr,
                                  Prob_Distr_Params=Prob_Distr_Params, 
                                  samplesize = as.integer(1),
                                  burnin=as.integer(1000000), 
                                  interval=as.integer(10),
                                  statsonly=TRUE,
                                  G=NULL,
                                  P=NULL,
                                  population=as.integer(population), 
                                  covPattern = as.integer(covPattern),
                                  bayesian_inference = FALSE,
                                  Ia = NULL, 
                                  Il = NULL, 
                                  R = NULL, 
                                  epi_params = NULL,
                                  print_calculations = FALSE) 
    
    G = CCMnet_Result[[1]]
    
    Initial_Data = Initialize_G_P_Ia_Il_R(population, beta_a, beta_l, gamma_a, gamma_l, G = G, num_init_infected = num_init_infected) 
    
    G = Initial_Data[[1]]
    P = Initial_Data[[2]]
    Ia = Initial_Data[[3]]
    Il = Initial_Data[[4]]
    R = Initial_Data[[5]]
    
    if (igraph::ecount(P) > 0) {
      invalid_epidemic = FALSE
      print("No Epidemic: restarting...")
    }
    
    counter = counter + 1
  }
  
  
  if (counter < max_try) {
  
    #####Calculate Mean and Variance of Graph properties#######
    
    if (Network_stats == 'Mixing') { 
      if (strong_prior == TRUE) {
        gamma_kappa = Prob_Distr_Params_hyperprior[[2]][1]
        gamma_theta = Prob_Distr_Params_hyperprior[[2]][2]
        
        alpha = Prob_Distr_Params_hyperprior[[3]]
        
        g_edgecount = igraph::gsize(G)
        
        sampled_v = sample(igraph::V(G), num_samples, replace=FALSE)
        
        sampled_v_deg =igraph::degree(G, sampled_v, mode = "all")
        g_edgecount = sum(sampled_v_deg) * .5
        
        sampled_covMatrix = matrix(rep(0, length(unique(covPattern))^2), 
                                   nrow = length(unique(covPattern)),
                                   ncol = length(unique(covPattern)))
        for (j in sampled_v) {
          neighbors_j = igraph::neighbors(G, j) #takes vertex id, returns names = vertex_id - 1
          neighbors_j = as.numeric(V(G)[neighbors_j])
          cov_j = covPattern[as.numeric(V(G)[j])] +1
          if (length(neighbors_j) > 0) {
            for (i in neighbors_j) {
              cov_i = covPattern[i] + 1
              sampled_covMatrix[cov_i, cov_j] = sampled_covMatrix[cov_i, cov_j] + 1
              if (cov_i != cov_j) {
                sampled_covMatrix[cov_j, cov_i] = sampled_covMatrix[cov_j, cov_i] + 1
              }
            }
          }
        }
        
        g_mixingcount = sampled_covMatrix[upper.tri(sampled_covMatrix, diag = TRUE)]
        
        Prob_Distr_Params_hyperprior[[2]][1] = g_edgecount
        Prob_Distr_Params_hyperprior[[2]][2] = population/num_samples
        
        Prob_Distr_Params_hyperprior[[3]] = alpha + g_mixingcount#as.numeric(CCMnet_Result[[2]])
        
      }
      Prob_Distr_Params[[1]][1] = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
      Prob_Distr_Params[[2]] = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
    }
    
    if (Network_stats == 'Degree') { 
      if (strong_prior == TRUE) {

        alpha = Prob_Distr_Params_hyperprior[[3]]
        sampled_v = sample(igraph::V(G), num_samples, replace=FALSE)
        sampled_v_deg =igraph::degree(G, sampled_v, mode = "all")
        
        if(reporting_bias == TRUE) {
          sampled_v_deg = gen_reporting_bias(sampled_v_deg)
        }
        
        g_deghist = (sampled_v_deg+1) %>% tabulate(nbins = population)

        Prob_Distr_Params_hyperprior[[3]] = alpha + g_deghist
        
      }
      Prob_Distr_Params[[1]][1] = 10
      Prob_Distr_Params[[2]] = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
      min_val = min(Prob_Distr_Params[[2]][which(Prob_Distr_Params[[2]]>0)])
      Prob_Distr_Params[[2]][which(Prob_Distr_Params[[2]]==0)] = min_val
      Prob_Distr_Params[[2]] = Prob_Distr_Params[[2]]/sum(Prob_Distr_Params[[2]])
    }
    
    PG_Data = Genetic_Seq_Data(P=P, Ia=Ia,final_time = max(Ia[which(Ia < Inf)]), genetic_bits = genetic_bits)
    T_dist = hamming.distance(PG_Data[,-(genetic_bits+1)])
    
    return(list(G, P, Ia, Il, R, T_dist, Prob_Distr_Params_hyperprior, Prob_Distr_Params, CCMnet_Result[[2]]))
  } else {
    print("ERROR: FAILED")
    return(list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL))
  }
}

SIIR.simulator <- function(g, population, beta_a, beta_l, gamma_a, gamma_l, num_init_infected = 1) {
  init_infected = sample(c(1:population),num_init_infected,replace=FALSE)
  
  Ia = array(Inf, dim=population)
  Il = array(Inf, dim=population)
  R = array(Inf, dim=population)
  Parent = array(NA, dim=population)
  
  for (init in init_infected) {
    Ia[init] = 0
    Il[init] = rexp(1, gamma_a)
    R[init] = Il[init] + rexp(1, gamma_l)
    Parent[init] = NA
  }
  
  Infected_Nodes = init_infected
  while(length(Infected_Nodes) > 0) {
    
    Node_index = min(which(Ia[Infected_Nodes] == min(Ia[Infected_Nodes]))) #identify and remove more current infected
    Node_id = Infected_Nodes[Node_index]
    Infected_Nodes = Infected_Nodes[-Node_index]
    
    Node_neighbors = which(g[Node_id,]==1)
    if (length(Node_neighbors) > 0) {
      for (poss_infect in Node_neighbors) {
        if (beta_a > 0) {
          acute_time = rexp(1, beta_a) #time of infectious contact during acute
          acute_phase_bool = ((acute_time + Ia[Node_id]) < Il[Node_id])  && ((acute_time + Ia[Node_id]) < Ia[poss_infect])
        } else {
          acute_phase_bool = FALSE
        }
        if (acute_phase_bool)  { #still in acute phase
          Ia[poss_infect] = acute_time + Ia[Node_id]
          Il[poss_infect] = rexp(1, gamma_a) + Ia[poss_infect]
          R[poss_infect] = Il[poss_infect] + rexp(1, gamma_l)
          Parent[poss_infect] = Node_id
          Infected_Nodes = unique(c(Infected_Nodes, poss_infect))
        } else {
          longterm_time = rexp(1, beta_l) #time of infectious contact during long-term
          if ( ((longterm_time + Il[Node_id]) < R[Node_id]) && ((longterm_time + Il[Node_id]) < Ia[poss_infect]) ) { #still in long-term phase
            Ia[poss_infect] = longterm_time + Il[Node_id]
            Il[poss_infect] = rexp(1, gamma_a) + Ia[poss_infect]
            R[poss_infect] = Il[poss_infect] + rexp(1, gamma_l)
            Parent[poss_infect] = Node_id
            Infected_Nodes = unique(c(Infected_Nodes, poss_infect))
          }
        }
      }
    }
  }
  infected_nodes = which(Ia < Inf)
  exampleepidemic = matrix(nrow = length(infected_nodes), ncol = 5)
  counter = 1
  for (infected_node in infected_nodes) {
    exampleepidemic[counter, ] = c(infected_node, Parent[infected_node], Ia[infected_node], Il[infected_node], R[infected_node])
    counter = counter + 1
  }
  colnames(exampleepidemic) = c("Infected", "Parent", "Acute", "Long-Term","Trt")
  return(exampleepidemic)
}


Initialize_G_P_Ia_Il_R <- function(population, beta_a, beta_l, gamma_a, gamma_l, G = NULL, num_init_infected = 1) {

  if (beta_a == 0) { 
    #Run SEIR model
    
    library('epinet')
    exampleepidemic <- SEIR.simulator(M = network::as.edgelist(intergraph::asNetwork(G)),
                                      N = population, beta = beta_l, 
                                      ki = 1, thetai = 1/gamma_l,
                                      ke = 1, thetae = 1/gamma_a,
                                      latencydist='gamma')
    detach("package:epinet")
    
  } else {
    exampleepidemic <- SIIR.simulator(G, population, beta_a, beta_l, gamma_a, gamma_l, num_init_infected)
  }
  
  nodes_attr_df = data.frame(name = c(0:(population-1)))
  p_df = data.frame(from = exampleepidemic[-which(is.na(exampleepidemic[,2])),2]-1, 
             to = exampleepidemic[-which(is.na(exampleepidemic[,2])),1]-1)
  P = graph_from_data_frame(p_df, directed=TRUE, vertices = nodes_attr_df)

  Ia = array(Inf, dim = population)
  for (i in c(1:length(exampleepidemic[,3]))) {
    Ia[exampleepidemic[i,1]] = exampleepidemic[i,3] 
  }
  
  R = array(Inf, dim = population)
  for (i in c(1:length(exampleepidemic[,3]))) {
    R[exampleepidemic[i,1]] = exampleepidemic[i,5]
  }	
  
  Il = array(Inf, dim = population)
  for (i in c(1:length(exampleepidemic[,3]))) {
    Il[exampleepidemic[i,1]] = exampleepidemic[i,4] #+ abs(exampleepidemic[1,3])
  }

  
  return(list(G,P,Ia,Il,R))
  
}


gen_reporting_bias <- function(sampled_v_deg) {
  
  orig_sampled_v_deg_total = sum(sampled_v_deg)
  bias_sampled_v_deg_total = round(orig_sampled_v_deg_total * ratio_orig_mean_deg)
  
  bias_sampled_v_deg = sampled_v_deg
  
  if (ratio_orig_mean_deg < 1) {
    num_remove = orig_sampled_v_deg_total - bias_sampled_v_deg_total
    for (i in c(1:num_remove)) {
      id = sample(c(1:length(sampled_v_deg)), 1, prob = bias_sampled_v_deg)
      bias_sampled_v_deg[id] = bias_sampled_v_deg[id] - 1
    }
  }
  
  if (ratio_orig_mean_deg > 1) {
    num_add = bias_sampled_v_deg_total - orig_sampled_v_deg_total
    for (i in c(1:num_add)) {
      id = sample(c(1:length(sampled_v_deg)), 1, prob = NULL)
      bias_sampled_v_deg[id] = bias_sampled_v_deg[id] + 1
    }
  }
  
  return(bias_sampled_v_deg)
}

