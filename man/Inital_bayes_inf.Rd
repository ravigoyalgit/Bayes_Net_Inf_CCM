% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesInf_v4_func_initial_bayes_inf.R
\name{Inital_bayes_inf}
\alias{Inital_bayes_inf}
\title{Inital_bayes_inf}
\usage{
Inital_bayes_inf(
  population,
  Network_stats,
  Prob_Distr,
  Prob_Distr_Params,
  covPattern,
  Ia,
  Il,
  R,
  beta_a,
  beta_l,
  gamma_a,
  gamma_l,
  T_dist,
  P_truth,
  init_P_truth_bool,
  G_truth,
  init_G_truth_bool
)
}
\arguments{
\item{population}{number of individuals in the network}

\item{Network_stats}{?}

\item{Prob_Distr}{?}

\item{Prob_Distr_Params}{?}

\item{covPattern}{?}

\item{Ia}{vector of times of infection for each individual}

\item{Il}{vector of times of transition to long-term infection for each individual}

\item{R}{vector of times of transition to recovery for each individual}

\item{beta_a}{probability of infection for acute phase}

\item{beta_l}{probability of infection for long-term phase}

\item{gamma_a}{length of infection for acute phase}

\item{gamma_l}{length of infection for long-term phase}

\item{T_dist}{genetic distances}

\item{P_truth}{the true transmission network}

\item{init_P_truth_bool}{Boolean, whether to start at the true transmission network; TRUE or FALSE}

\item{G_truth}{the true contact network}

\item{init_G_truth_bool}{where to start at the true contact network; 0 is random, 1 is P_start, 2 is G_truth}
}
\value{
list(G_start, P_start), starting contact and transmission networks
}
\description{
setup starting contact and transmission networks for Bayesian inference
}
