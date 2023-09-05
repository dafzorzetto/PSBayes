###################################################################################
# ---             SIMULATION STUDY:             ---
# ---                5 settings                 ---
# ---         PRINCIPAL STRATIFICATION          ---
###################################################################################

# upload functions
source("src/simulation_functions.R")
source("src/plot_simulation.R")

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
                       beta_1=c(1,2,-1,0.5),
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
                       beta_1=c(1,1.2,-1,1),
                       sigma_y=c(-0.5,0.1)))

# case 1 with closer strata and different variance
scenario_3=lapply(1:samples, function(s) 
  prova=setup_sim_2cov(seed=s,
                       eta=c(1.5,2,2.5),
                       sigma_p=c(0.12,0.1,0.08),
                       allocation_0=c(1,2,2),
                       allocation_1=c(1,3,3),
                       beta_0=c(1,2),
                       beta_1=c(1,1.2,-0.8,0.5),
                       sigma_y=c(-0.5,0.1)))

# case 2 with closer strata and different variance
scenario_4=lapply(1:samples, function(s) 
  prova=setup_sim_2cov(seed=s,
                       eta=c(1.5,2,2.5),
                       sigma_p=c(0.12,0.1,0.08),
                       allocation_0=c(2,2,2),
                       allocation_1=c(1,2,3),
                       beta_0=c(1,2),
                       beta_1=c(1,1.2,-0.8,0.5),
                       sigma_y=c(-0.5,0.1)))

# case 2 with 5 covariates
scenario_5=lapply(1:samples, function(s) 
  prova=setup_sim_5cov(seed=s,
                       eta=c(1,2,3,4),
                       sigma_p=rep(0.05,4),
                       allocation_0=c(2,2,3),
                       allocation_1=c(1,2,4),
                       beta_0=c(1,2),
                       beta_1=c(1,1.2,-1,0.5),
                       sigma_y=c(-0.5,0.1)))


save.image("data_simulations.RData")

###################################################################################
# ---     plots: histograms of post-treatment or outcome distributions     ---
###################################################################################

# POST_TREATMENT variable: divider by treatment level
hist_P_sim(data=scenario_1[[1]],s_n=1)
hist_P_sim(data=scenario_2[[1]],s_n=2)
hist_P_sim(data=scenario_3[[1]],s_n=3)
hist_P_sim(data=scenario_4[[1]],s_n=4)
hist_P_sim(data=scenario_5[[1]],s_n=5)

# POST_TREATMENT variable: divider by treatment level
hist_diff_P_sim(data=scenario_1[[1]],s_n=1)
hist_diff_P_sim(data=scenario_2[[1]],s_n=2)
hist_diff_P_sim(data=scenario_3[[1]],s_n=3)
hist_diff_P_sim(data=scenario_4[[1]],s_n=4)
hist_diff_P_sim(data=scenario_5[[1]],s_n=5)

# OUTCOME variable: divider by treatment level
hist_Y_sim(data=scenario_1[[1]],s_n=1)
hist_Y_sim(data=scenario_2[[1]],s_n=2)
hist_Y_sim(data=scenario_3[[1]],s_n=3)
hist_Y_sim(data=scenario_4[[1]],s_n=4)
hist_Y_sim(data=scenario_5[[1]],s_n=5)

# POST_TREATMENT variable: divider by treatment level
hist_diff_Y_sim(data=scenario_1[[1]],s_n=1)
hist_diff_Y_sim(data=scenario_2[[1]],s_n=2)
hist_diff_Y_sim(data=scenario_3[[1]],s_n=3)
hist_diff_Y_sim(data=scenario_4[[1]],s_n=4)
hist_diff_Y_sim(data=scenario_5[[1]],s_n=5)


