#' @title Genetic_Seq_Data
#' @description generates phylogenetic data by assigning genetic sequences to all individuals who are infected at the end of simulation
#' @param P transmission network
#' @param Ia vector of times of infection for each individual
#' @param final_time latest observed infection time
#' @param genetic_bits length of genetic sequences to generate
#' @return Genetic_Seq_mat, matrix of genetic sequences; rows are individuals, columns are genetic bits
#' @export

Genetic_Seq_Data <- function(P, Ia, final_time = max(Ia[which(Ia < Inf)]), genetic_bits = 128) {

  infection_order = order(Ia)
  Genetic_Seq_mat = matrix(0, nrow = length(Ia), ncol = genetic_bits)
  update_times_a = array(0, dim = length(Ia))
  Genetic_Seq_mat = cbind(Genetic_Seq_mat, update_times_a)
  for(i in infection_order[c(1:length(which(Ia < Inf)))]) {
    parent_i = which(P[,i]==1)
    if (length(parent_i) == 0) { #root node
      Genetic_Seq_mat[i,] = c(genetic_seq_root_node(genetic_bits), Ia[i])
    } else { #has a parent
      Genetic_Seq_mat[parent_i,] = c(genetic_seq_parent_node(Genetic_Seq_mat[parent_i,-(genetic_bits+1)], Genetic_Seq_mat[parent_i,(genetic_bits+1)], Ia[i], genetic_bits), Ia[i])
      Genetic_Seq_mat[i,] = c(genetic_seq_inf_node(Genetic_Seq_mat[parent_i,-(genetic_bits+1)], Genetic_Seq_mat[parent_i,(genetic_bits+1)], Ia[i], genetic_bits), Ia[i])
    }
  }

  # Update everyone to be sampled at end of trial

  for(i in infection_order[c(1:length(which(Ia < Inf)))]) {
    Genetic_Seq_mat[i,] = c(genetic_seq_parent_node(Genetic_Seq_mat[i,-(genetic_bits+1)], Genetic_Seq_mat[i,(genetic_bits+1)], final_time, genetic_bits), final_time)
  }
  return(Genetic_Seq_mat)

}

#' @title genetic_seq_root_node
#' @description generates a genetic sequence for an infected individual that has no parent by randomly sampling each genetic bit
#' @param genetic_bits length of genetic sequence to generate
#' @return sample(c(0,1), genetic_bits, replace = TRUE), a genetic sequence
#' @export

genetic_seq_root_node <- function(genetic_bits) {
  return(sample(c(0,1), genetic_bits, replace = TRUE))
}

#' @title genetic_seq_parent_node
#' @description generates a genetic sequence for an infected individual by mutating genetic bits proportional to how long infected
#' @param Genetic_Seq_a genetic sequence that is inherited, to be mutated
#' @param last_updatetime last update time
#' @param ctime current time
#' @param genetic_bits length of genetic sequences to generate
#' @return Genetic_Seq_a, a mutated genetic sequence
#' @export

genetic_seq_parent_node <- function(Genetic_Seq_a, last_updatetime, ctime, genetic_bits) {
  nbit = round(ctime - last_updatetime)
  change_bits = array(0, dim = genetic_bits)
  change_bits[sample(c(1:genetic_bits), nbit, replace = FALSE)] = 1
  Genetic_Seq_a = (Genetic_Seq_a + change_bits) %% 2
  return(Genetic_Seq_a)
}

#' @title genetic_seq_inf_node
#' @description generates a genetic sequence for an infected individual by changing a fixed number of bits
#' @param Genetic_Seq_a genetic sequence that is inherited, to be mutated
#' @param last_updatetime last update time
#' @param ctime current time
#' @param genetic_bits length of genetic sequences to generate
#' @return Genetic_Seq_a, a mutated genetic sequence
#' @export

genetic_seq_inf_node <- function(Genetic_Seq_a, last_updatetime, ctime, genetic_bits) {
  nbit = 10
  change_bits = array(0, dim = genetic_bits)
  change_bits[sample(c(1:genetic_bits), nbit, replace = FALSE)] = 1
  Genetic_Seq_a = (Genetic_Seq_a + change_bits) %% 2
  return(Genetic_Seq_a)
}
