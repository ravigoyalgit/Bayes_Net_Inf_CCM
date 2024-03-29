% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesInf_v4_func_gen_genetic_data.R
\name{genetic_seq_root_node}
\alias{genetic_seq_root_node}
\title{genetic_seq_root_node}
\usage{
genetic_seq_root_node(genetic_bits)
}
\arguments{
\item{genetic_bits}{length of genetic sequence to generate}
}
\value{
sample(c(0,1), genetic_bits, replace = TRUE), a genetic sequence
}
\description{
generates a genetic sequence for an infected individual that has no parent by randomly sampling each genetic bit
}
