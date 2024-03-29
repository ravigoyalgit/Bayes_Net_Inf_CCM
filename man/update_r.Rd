% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesInf_v4_func_bayes_inf.R
\name{update_r}
\alias{update_r}
\title{update_r}
\usage{
update_r(r, g, P, Ia, Il, R_times, edge1, beta_a, beta_l)
}
\arguments{
\item{r}{?}

\item{g}{?}

\item{P}{transmission network}

\item{Ia}{vector of times of infection for each individual}

\item{Il}{vector of times of transition to long-term infection for each individual}

\item{R_times}{?}

\item{edge1}{?}

\item{beta_a}{probability of infection for acute phase}

\item{beta_l}{probability of infection for long-term phase}
}
\value{
p_noinfect / (1 - p_noinfect) or (1-p_noinfect) / p_noinfect, probability of acceptance or rejection
}
\description{
gives the probability of the Metropolis-Hastings algorithm accepting or rejecting the proposal network
}
