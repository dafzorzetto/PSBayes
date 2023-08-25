###################################################################################
# ---     PRINCIPAL STRATIFICATION      ---
# ---        P-model + Y-model          ---
# ---         analysis results          ---
###################################################################################

#libraries
library(mcclust)

#load data
load("estimated_models.RData")

###################################################################################

# --- PARTITIONS: strata ---
clusters_strata_1=lapply(1:samples,function(c) as.vector(table(CASDMM_1[[c]]$S_strata_cluster)))
clusters_strata_2=lapply(1:samples,function(c) as.vector(table(CASDMM_2[[c]]$S_strata_cluster)))
clusters_strata_3=lapply(1:samples,function(c) as.vector(table(CASDMM_3[[c]]$S_strata_cluster)))
clusters_strata_4=lapply(1:samples,function(c) as.vector(table(CASDMM_4[[c]]$S_strata_cluster)))
clusters_strata_5=lapply(1:samples,function(c) as.vector(table(CASDMM_5[[c]]$S_strata_cluster)))

# --- adjasted RAND index (ARI) ---

ARI_1_CASDMM<-sapply(1:samples, function(c)
  arandi(CASDMM_1[[c]]$S_strata_cluster, scenario_1[[c]]$clusters$S_strata))
ARI_2_CASDMM<-sapply(1:samples, function(c)
  arandi(CASDMM_2[[c]]$S_strata_cluster, scenario_2[[c]]$clusters$S_strata))
ARI_3_CASDMM<-sapply(1:samples, function(c)
  arandi(CASDMM_3[[c]]$S_strata_cluster, scenario_3[[c]]$clusters$S_strata))
ARI_4_CASDMM<-sapply(1:samples, function(c)
  arandi(CASDMM_4[[c]]$S_strata_cluster, scenario_4[[c]]$clusters$S_strata))
ARI_5_CASDMM<-sapply(1:samples, function(c)
  arandi(CASDMM_5[[c]]$S_strata_cluster, scenario_5[[c]]$clusters$S_strata))

ARI<-matrix(c(mean(ARI_1_CASDMM),sd(ARI_1_CASDMM),
              mean(ARI_2_CASDMM),sd(ARI_2_CASDMM),
              mean(ARI_3_CASDMM),sd(ARI_3_CASDMM),
              mean(ARI_4_CASDMM),sd(ARI_4_CASDMM),
              mean(ARI_5_CASDMM),sd(ARI_5_CASDMM)),ncol=5)
row.names(ARI)<-c("mean","sd")
colnames(ARI)<-paste0("scenario ",1:5)
round(ARI,4)

###################################################################################

# --- bias and MSE for P: post-treatment ---

bias_P_CASDMM1=sapply(1:samples, function(c) 
  CASDMM_1[[c]]$post_P_1_imp-CASDMM_1[[c]]$post_P_0_imp-scenario_1[[c]]$data$P_1+scenario_1[[c]]$data$P_0)
bias_P_CASDMM2=sapply(1:samples, function(c) 
  CASDMM_2[[c]]$post_P_1_imp-CASDMM_2[[c]]$post_P_0_imp-scenario_2[[c]]$data$P_1+scenario_2[[c]]$data$P_0)
bias_P_CASDMM3=sapply(1:samples, function(c) 
  CASDMM_3[[c]]$post_P_1_imp-CASDMM_3[[c]]$post_P_0_imp-scenario_3[[c]]$data$P_1+scenario_3[[c]]$data$P_0)
bias_P_CASDMM4=sapply(1:samples, function(c) 
  CASDMM_4[[c]]$post_P_1_imp-CASDMM_4[[c]]$post_P_0_imp-scenario_4[[c]]$data$P_1+scenario_4[[c]]$data$P_0)
bias_P_CASDMM5=sapply(1:samples, function(c) 
  CASDMM_5[[c]]$post_P_1_imp-CASDMM_5[[c]]$post_P_0_imp-scenario_5[[c]]$data$P_1+scenario_5[[c]]$data$P_0)

mse_P_CASDMM1=sapply(1:samples, function(c) 
  (CASDMM_1[[c]]$post_P_1_imp-CASDMM_1[[c]]$post_P_0_imp-scenario_1[[c]]$data$P_1+scenario_1[[c]]$data$P_0)^2)
mse_P_CASDMM2=sapply(1:samples, function(c) 
  (CASDMM_2[[c]]$post_P_1_imp-CASDMM_2[[c]]$post_P_0_imp-scenario_2[[c]]$data$P_1+scenario_2[[c]]$data$P_0)^2)
mse_P_CASDMM3=sapply(1:samples, function(c) 
  (CASDMM_3[[c]]$post_P_1_imp-CASDMM_3[[c]]$post_P_0_imp-scenario_3[[c]]$data$P_1+scenario_3[[c]]$data$P_0)^2)
mse_P_CASDMM4=sapply(1:samples, function(c) 
  (CASDMM_4[[c]]$post_P_1_imp-CASDMM_4[[c]]$post_P_0_imp-scenario_4[[c]]$data$P_1+scenario_4[[c]]$data$P_0)^2)
mse_P_CASDMM5=sapply(1:samples, function(c) 
  (CASDMM_5[[c]]$post_P_1_imp-CASDMM_5[[c]]$post_P_0_imp-scenario_5[[c]]$data$P_1+scenario_5[[c]]$data$P_0)^2)

###################################################################################

# --- bias and MSE for Y: outcome ---

bias_Y_CASDMM1=sapply(1:samples, function(c) 
  CASDMM_1[[c]]$post_Y_1_imp-CASDMM_1[[c]]$post_Y_0_imp-scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)
bias_Y_CASDMM2=sapply(1:samples, function(c) 
  CASDMM_2[[c]]$post_Y_1_imp-CASDMM_2[[c]]$post_Y_0_imp-scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)
bias_Y_CASDMM3=sapply(1:samples, function(c) 
  CASDMM_3[[c]]$post_Y_1_imp-CASDMM_3[[c]]$post_Y_0_imp-scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)
bias_Y_CASDMM4=sapply(1:samples, function(c) 
  CASDMM_4[[c]]$post_Y_1_imp-CASDMM_4[[c]]$post_Y_0_imp-scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)
bias_Y_CASDMM5=sapply(1:samples, function(c) 
  CASDMM_5[[c]]$post_Y_1_imp-CASDMM_5[[c]]$post_Y_0_imp-scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)

mse_Y_CASDMM1=sapply(1:samples, function(c) 
  (CASDMM_1[[c]]$post_Y_1_imp-CASDMM_1[[c]]$post_Y_0_imp-scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)^2)
mse_Y_CASDMM2=sapply(1:samples, function(c) 
  (CASDMM_2[[c]]$post_Y_1_imp-CASDMM_2[[c]]$post_Y_0_imp-scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)^2)
mse_Y_CASDMM3=sapply(1:samples, function(c) 
  (CASDMM_3[[c]]$post_Y_1_imp-CASDMM_3[[c]]$post_Y_0_imp-scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)^2)
mse_Y_CASDMM4=sapply(1:samples, function(c) 
  (CASDMM_4[[c]]$post_Y_1_imp-CASDMM_4[[c]]$post_Y_0_imp-scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)^2)
mse_Y_CASDMM5=sapply(1:samples, function(c) 
  (CASDMM_5[[c]]$post_Y_1_imp-CASDMM_5[[c]]$post_Y_0_imp-scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)^2)


#########################################################################
save.image("estimands.RData")

