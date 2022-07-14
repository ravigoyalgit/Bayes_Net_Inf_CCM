#' @title Inital_bayes_inf
#' @description setup starting contact and transmission networks for Bayesian inference
#' @param population number of individuals in the network
#' @param Network_stats ?
#' @param Prob_Distr ?
#' @param Prob_Distr_Params ?
#' @param covPattern ?
#' @param Ia vector of times of infection for each individual
#' @param Il vector of times of transition to long-term infection for each individual
#' @param R vector of times of transition to recovery for each individual
#' @param beta_a probability of infection for acute phase
#' @param beta_l probability of infection for long-term phase
#' @param gamma_a length of infection for acute phase
#' @param gamma_l length of infection for long-term phase
#' @param T_dist genetic distances
#' @param P_truth the true transmission network
#' @param init_P_truth_bool Boolean, whether to start at the true transmission network; TRUE or FALSE
#' @param G_truth the true contact network
#' @param init_G_truth_bool where to start at the true contact network; 0 is random, 1 is P_start, 2 is G_truth
#' @return list(G_start, P_start), starting contact and transmission networks
#' @export

Inital_bayes_inf <- function(population,Network_stats, Prob_Distr, Prob_Distr_Params,
                             covPattern,
                             Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist,
                             P_truth, init_P_truth_bool, G_truth, init_G_truth_bool) {

  if (init_P_truth_bool) {
    P_start = P_truth
  } else {
    G_full = igraph::graph.full(n = population, directed = FALSE)
    G_full = G_full %>% igraph::set_vertex_attr("name", value = c(0:(population-1)))
    P_start = Update_P(G_full,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)
  }

  if (init_G_truth_bool == 0) {

    CCMnet_Result = CCMnetpy::CCMnet_constr(Network_stats=Network_stats,
                                            Prob_Distr=Prob_Distr,
                                            Prob_Distr_Params=Prob_Distr_Params,
                                            samplesize = as.integer(1),
                                            burnin=as.integer(100000),
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

    G_start2 = CCMnet_Result[[1]]

    U_P_Start = igraph::as.undirected(P_start, mode = "collapse")
    igraph::V(U_P_Start)$name <- as.character(c(0:(population-1)))
    G_start = igraph::union(G_start2, U_P_Start)
  } else if (init_G_truth_bool == 1) {
    G_start = igraph::as.undirected(P_start, mode = "collapse")
  } else {
    G_start = G_truth
  }

  return(list(G_start, P_start))
}
