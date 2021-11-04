#function to read data and generate primary school network from http://www.sociopatterns.org/datasets/primary-school-cumulative-networks/

#required packages
library(tidyverse)
library(rgexf)
library(network)
library(intergraph)
library(ergm)

gen_school_net <- function() {
  
  dir_loc_data = "/Users/sallyslipher/Desktop/MSU Research Associate/Projects/Bayesian Modeling/Paper 1 Network Application/"
  
  school.df <- read_tsv(paste0(dir_loc_data, "metadata.txt"), col_names = FALSE)
  names(school.df) = c("name", "classroom", "gender")
  
  school.df = school.df %>% 
    mutate(role = if_else(classroom == "Teachers", "Teacher", "Student"))
  
  school.df$id = c(1:nrow(school.df))
  
  teacher.df = school.df %>% 
    filter(role == "Teacher")
  
  school.gexf <- read.gexf(paste0(dir_loc_data, "sp_data_school_day_1_g.gexf"))
  
  school.ig <- gexf.to.igraph(school.gexf)
  school.net <- asNetwork(school.ig)
  
  network::set.vertex.attribute(school.net, "classroom", value = school.df$classroom)
  network::set.vertex.attribute(school.net, "gender", value = school.df$gender)
  network::set.vertex.attribute(school.net, "role", value = school.df$role)
  
  #network::get.vertex.attribute(school.net, "vertex.names") Same order as data.frame
  
  for (i in c(1:nrow(teacher.df))) {
    i = 8
    classcount.df = network::get.vertex.attribute(school.net, "classroom")[
      get.neighborhood(school.net, v = teacher.df$id[i], type = c("combined"), na.omit = TRUE)] %>% 
      table() %>%
      as_tibble() 
    
    names(classcount.df) = c("classroom", "n") 
    
    school.df$classroom[teacher.df$id[i]] = classcount.df %>% 
      slice_max(n = 1, order_by = n, with_ties = FALSE) %>% pull(classroom) %>% as.character()
  }
  
  network::set.vertex.attribute(school.net, "classroom", value = school.df$classroom)
  
  return(school.net)
  
}

#generate network
school_net <- gen_school_net()

#network summary
school_net
plot(school_net)
summary(school_net) #vertex attributes

#density
network.density(school_net)

#network size
network.size(school_net)

#degree distribution
sna::degree(school_net)/2 #network package counts both directions = doubled
hist(sna::degree(school_net)/2)

#mean degree
mean(sna::degree(school_net)/2)

#assortativity
#knitr::kable(summary(school_net ~ nodemix('classroom'))) #table
summary(school_net ~ nodemix('classroom')) #mixing of classrooms
