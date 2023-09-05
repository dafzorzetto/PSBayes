###################################################################################
# ---             ESTIMATION MODELS:             ---
# ---                5 settings                 ---
# ---  1- CONFONDER-DEPENDENT SHARED ATOMS MODEL   
# ---  2- competitor
###################################################################################

# upload functions
source("src/CASDMM.R")

# upload data
load("data_simulations.RData")

###################################################################################
# --- general parameters for Gibbs samplers ----

# iterations
R=3000
R_burnin=1500

# max number clusters
n_cluster=10

#########################################

# estimation CASDMM for the 5 scenarios
#CASDMM_1=mclapply(1:samples, Gibbs_CASDMM, sim=scenario_1, mc.cores = 6)
CASDMM_2=mclapply(1:samples, Gibbs_CASDMM, sim=scenario_2, mc.cores = 6)
CASDMM_3=mclapply(1:samples, Gibbs_CASDMM, sim=scenario_3, mc.cores = 6)
CASDMM_4=mclapply(1:samples, Gibbs_CASDMM, sim=scenario_4, mc.cores = 6)
#CASDMM_5=mclapply(1:samples, Gibbs_CASDMM, sim=scenario_5, mc.cores = 6)

save.image("estimated_models.RData")