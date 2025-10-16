###################################################################################
# ---             ESTIMATION MODELS:             ---
# ---                5 settings                 ---
# ---  1- CONFONDER-DEPENDENT SHARED ATOMS MODEL   
# ---  2- competitor
###################################################################################

# upload functions
source("src/CASBAH.R")
source("src/competitor_SLM.R")

# upload data
load("data_simulations.RData")

# code for LJD (modified to obtain Y)
source("part_I_all_functions.R")

###################################################################################
# --- general parameters for Gibbs samplers ----

# iterations
R=3000
R_burnin=1500

# max number clusters
n_cluster=10

#########################################

# estimation proposed model for the 5 scenarios
CASDMM_1 <- mclapply(1:samples, Gibbs_CASDMM, sim=scenario_1, mc.cores = 5)
CASDMM_2 <- mclapply(1:samples, Gibbs_CASDMM, sim=scenario_2, mc.cores = 5)
CASDMM_3 <- mclapply(1:samples, Gibbs_CASDMM, sim=scenario_3, mc.cores = 5)
CASDMM_4 <- mclapply(1:samples, Gibbs_CASDMM, sim=scenario_4, mc.cores = 5)
CASDMM_5 <- mclapply(1:samples, Gibbs_CASDMM, sim=scenario_5, mc.cores = 5)

save.image("estimated_models.RData")

# estimation SLM model 
SLM_1 <- mclapply(1:samples, SLM_GIbbs, sim=scenario_1, mc.cores = 6)
SLM_2 <- mclapply(1:samples, SLM_GIbbs, sim=scenario_2, mc.cores = 6)
SLM_3 <- mclapply(1:samples, SLM_GIbbs, sim=scenario_3, mc.cores = 6)
SLM_4 <- mclapply(1:samples, SLM_GIbbs, sim=scenario_4, mc.cores = 6)
SLM_5 <- mclapply(1:samples, SLM_GIbbs, sim=scenario_5, mc.cores = 6)

save.image("estimated_models_SLM.RData")

# estimation LJD model 
LJD_1=lapply(1:samples,function(i) point_estimator(Z=scenario_1[[i]]$data$T, 
                  X=data.frame(scenario_1[[i]]$data$X), 
                  S=scenario_1[[i]]$data$T*(scenario_1[[i]]$data$P_1)+
                    (1-scenario_1[[i]]$data$T)*(scenario_1[[i]]$data$P_0), 
                  Y=scenario_1[[i]]$data$T*(scenario_1[[i]]$data$Y_1)+
                    (1-scenario_1[[i]]$data$T)*(scenario_1[[i]]$data$Y_0), 
                  n_divisions=100, copula_type='gaussian', rho=0.25))
LJD_2=lapply(1:samples,function(i) point_estimator(Z=scenario_2[[i]]$data$T, 
                                                   X=data.frame(scenario_2[[i]]$data$X), 
                                                   S=scenario_2[[i]]$data$T*(scenario_2[[i]]$data$P_1)+
                                                     (1-scenario_2[[i]]$data$T)*(scenario_2[[i]]$data$P_0), 
                                                   Y=scenario_2[[i]]$data$T*(scenario_2[[i]]$data$Y_1)+
                                                     (1-scenario_2[[i]]$data$T)*(scenario_2[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0.25))
LJD_3=lapply(1:samples,function(i) point_estimator(Z=scenario_3[[i]]$data$T, 
                                                   X=data.frame(scenario_3[[i]]$data$X), 
                                                   S=scenario_3[[i]]$data$T*(scenario_3[[i]]$data$P_1)+
                                                     (1-scenario_3[[i]]$data$T)*(scenario_3[[i]]$data$P_0), 
                                                   Y=scenario_3[[i]]$data$T*(scenario_3[[i]]$data$Y_1)+
                                                     (1-scenario_3[[i]]$data$T)*(scenario_3[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0.25))
LJD_4=lapply(1:samples,function(i) point_estimator(Z=scenario_4[[i]]$data$T, 
                                                   X=data.frame(scenario_4[[i]]$data$X), 
                                                   S=scenario_4[[i]]$data$T*(scenario_4[[i]]$data$P_1)+
                                                     (1-scenario_4[[i]]$data$T)*(scenario_4[[i]]$data$P_0), 
                                                   Y=scenario_4[[i]]$data$T*(scenario_4[[i]]$data$Y_1)+
                                                     (1-scenario_4[[i]]$data$T)*(scenario_4[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0.25))
LJD_5=lapply(1:samples,function(i) point_estimator(Z=scenario_5[[i]]$data$T, 
                                                   X=data.frame(scenario_5[[i]]$data$X), 
                                                   S=scenario_5[[i]]$data$T*(scenario_5[[i]]$data$P_1)+
                                                     (1-scenario_5[[i]]$data$T)*(scenario_5[[i]]$data$P_0), 
                                                   Y=scenario_5[[i]]$data$T*(scenario_5[[i]]$data$Y_1)+
                                                     (1-scenario_5[[i]]$data$T)*(scenario_5[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0.25))


# estimation LJD model  rho=0
LJD_1_0=lapply(1:samples,function(i) point_estimator(Z=scenario_1[[i]]$data$T, 
                                                   X=data.frame(scenario_1[[i]]$data$X), 
                                                   S=scenario_1[[i]]$data$T*(scenario_1[[i]]$data$P_1)+
                                                     (1-scenario_1[[i]]$data$T)*(scenario_1[[i]]$data$P_0), 
                                                   Y=scenario_1[[i]]$data$T*(scenario_1[[i]]$data$Y_1)+
                                                     (1-scenario_1[[i]]$data$T)*(scenario_1[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0))
LJD_2_0=lapply(1:samples,function(i) point_estimator(Z=scenario_2[[i]]$data$T, 
                                                   X=data.frame(scenario_2[[i]]$data$X), 
                                                   S=scenario_2[[i]]$data$T*(scenario_2[[i]]$data$P_1)+
                                                     (1-scenario_2[[i]]$data$T)*(scenario_2[[i]]$data$P_0), 
                                                   Y=scenario_2[[i]]$data$T*(scenario_2[[i]]$data$Y_1)+
                                                     (1-scenario_2[[i]]$data$T)*(scenario_2[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0))
LJD_3_0=lapply(1:samples,function(i) point_estimator(Z=scenario_3[[i]]$data$T, 
                                                   X=data.frame(scenario_3[[i]]$data$X), 
                                                   S=scenario_3[[i]]$data$T*(scenario_3[[i]]$data$P_1)+
                                                     (1-scenario_3[[i]]$data$T)*(scenario_3[[i]]$data$P_0), 
                                                   Y=scenario_3[[i]]$data$T*(scenario_3[[i]]$data$Y_1)+
                                                     (1-scenario_3[[i]]$data$T)*(scenario_3[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0))
LJD_4_0=lapply(1:samples,function(i) point_estimator(Z=scenario_4[[i]]$data$T, 
                                                   X=data.frame(scenario_4[[i]]$data$X), 
                                                   S=scenario_4[[i]]$data$T*(scenario_4[[i]]$data$P_1)+
                                                     (1-scenario_4[[i]]$data$T)*(scenario_4[[i]]$data$P_0), 
                                                   Y=scenario_4[[i]]$data$T*(scenario_4[[i]]$data$Y_1)+
                                                     (1-scenario_4[[i]]$data$T)*(scenario_4[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0))
LJD_5_0=lapply(1:samples,function(i) point_estimator(Z=scenario_5[[i]]$data$T, 
                                                   X=data.frame(scenario_5[[i]]$data$X), 
                                                   S=scenario_5[[i]]$data$T*(scenario_5[[i]]$data$P_1)+
                                                     (1-scenario_5[[i]]$data$T)*(scenario_5[[i]]$data$P_0), 
                                                   Y=scenario_5[[i]]$data$T*(scenario_5[[i]]$data$Y_1)+
                                                     (1-scenario_5[[i]]$data$T)*(scenario_5[[i]]$data$Y_0), 
                                                   n_divisions=100, copula_type='gaussian', rho=0))

save.image("estimated_models_LJD.RData")