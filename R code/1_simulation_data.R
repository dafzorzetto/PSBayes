###################################################################################
# ---             SIMULATION STUDY:             ---
# ---                5 settings                 ---
# ---         PRINCIPAL STRATIFICATION          ---
###################################################################################

#library
#library(ggplot2)
#library(ggExtra)

# upload functions
source("src/simulation_functions.R")

###################################################################################

# --- generation parameters ---

n=500          # units for each samples
samples=100    # number of samples

# --- 5 settings ---

# 2 confounders
# 2 strata: 1 DIS, 1 ASS+                                                                                                                                                          in common - allocation is shifted
scenario_1=lapply(1:samples, function(s) 
  prova=setup_sim_2cov(seed=s,
                       eta=c(1,2,3),
                       sigma_p=rep(0.05,3),
                       allocation_0=c(1,2,2),
                       allocation_1=c(1,3,3),
                       beta_0=c(1,2),
                       beta_1=c(1,2,-1,1.5),
                       sigma_y=c(-0.5,0.1)))

# 2 confounders
# 3 strata: 1 DIS, 1 ASS+, 1 ASS-  
scenario_2=lapply(1:samples, function(s) 
  prova=setup_sim_2cov(seed=s,
                       eta=c(1,2,3),
                       sigma_p=rep(0.05,3),
                       allocation_0=c(2,2,2),
                       allocation_1=c(1,2,3),
                       beta_0=c(1,2),
                       beta_1=c(1,2,-1,1.5),
                       sigma_y=c(-0.5,0.1)))

# case 1 with closer strata and different variance
scenario_3=lapply(1:samples, function(s) 
  prova=setup_sim_2cov(seed=s,
                       eta=c(1.5,2,2.75),
                       sigma_p=c(0.2,0.1,0.15),
                       allocation_0=c(1,2,3),
                       allocation_1=c(2,3,3),
                       beta_0=c(1,2),
                       beta_1=c(1,2,-1,1.5),
                       sigma_y=c(-0.5,0.1)))

# case 2 with closer strata and different variance
scenario_4=lapply(1:samples, function(s) 
  prova=setup_sim_2cov(seed=s,
                       eta=c(1.5,2,2.5),
                       sigma_p=c(0.2,0.1,0.15),
                       allocation_0=c(2,2,2),
                       allocation_1=c(1,2,3),
                       beta_0=c(1,2),
                       beta_1=c(1,2,-1,1.5),
                       sigma_y=c(-0.5,0.1)))

# case 2 with 5 covariates
scenario_5=lapply(1:samples, function(s) 
  prova=setup_sim_5cov(seed=s,
                       eta=c(1,2,3),
                       sigma_p=rep(0.05,3),
                       allocation_0=c(2,2,2),
                       allocation_1=c(1,2,3),
                       beta_0=c(1,2),
                       beta_1=c(1,2,-1,1.5),
                       sigma_y=c(-0.5,0.1)))


save.image("data_simulations.RData")

