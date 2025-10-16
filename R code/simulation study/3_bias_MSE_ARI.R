###################################################################################
# ---     PRINCIPAL STRATIFICATION      ---
# ---        P-model + Y-model          ---
# ---         analysis results          ---
###################################################################################

#libraries
library(mcclust)

#load data
load("estimated_models.RData")
load("sim_LJD.RData")

###################################################################################

# --- PARTITIONS: strata ---
clusters_strata_1=lapply(1:samples,function(c) as.vector(table(CASDMM_1[[c]]$S_strata_cluster)))
clusters_strata_2=lapply(1:samples,function(c) as.vector(table(CASDMM_2[[c]]$S_strata_cluster)))
clusters_strata_3=lapply(1:samples,function(c) as.vector(table(CASDMM_3[[c]]$S_strata_cluster)))
clusters_strata_4=lapply(1:samples,function(c) as.vector(table(CASDMM_4[[c]]$S_strata_cluster)))
clusters_strata_5=lapply(1:samples,function(c) as.vector(table(CASDMM_5[[c]]$S_strata_cluster)))

table(sapply(1:samples, function(c) length(clusters_strata_1[[c]])))
table(sapply(1:samples, function(c) length(clusters_strata_2[[c]])))
table(sapply(1:samples, function(c) length(clusters_strata_3[[c]])))
table(sapply(1:samples, function(c) length(clusters_strata_4[[c]])))
table(sapply(1:samples, function(c) length(clusters_strata_5[[c]])))

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

# CASBAH
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

# Schartz Li Mealli model
bias_P_SLM1=sapply(1:samples, function(c) 
  SLM_1[[c]]$post_P_1_imp-SLM_1[[c]]$post_P_0_imp-scenario_1[[c]]$data$P_1+scenario_1[[c]]$data$P_0)
bias_P_SLM2=sapply(1:samples, function(c) 
  SLM_2[[c]]$post_P_1_imp-SLM_2[[c]]$post_P_0_imp-scenario_2[[c]]$data$P_1+scenario_2[[c]]$data$P_0)
bias_P_SLM3=sapply(1:samples, function(c) 
  SLM_3[[c]]$post_P_1_imp-SLM_3[[c]]$post_P_0_imp-scenario_3[[c]]$data$P_1+scenario_3[[c]]$data$P_0)
bias_P_SLM4=sapply(1:samples, function(c) 
  SLM_4[[c]]$post_P_1_imp-SLM_4[[c]]$post_P_0_imp-scenario_4[[c]]$data$P_1+scenario_4[[c]]$data$P_0)
bias_P_SLM5=sapply(1:samples, function(c) 
  SLM_5[[c]]$post_P_1_imp-SLM_5[[c]]$post_P_0_imp-scenario_5[[c]]$data$P_1+scenario_5[[c]]$data$P_0)

mse_P_SLM1=sapply(1:samples, function(c) 
  (SLM_1[[c]]$post_P_1_imp-SLM_1[[c]]$post_P_0_imp-scenario_1[[c]]$data$P_1+scenario_1[[c]]$data$P_0)^2)
mse_P_SLM2=sapply(1:samples, function(c) 
  (SLM_2[[c]]$post_P_1_imp-SLM_2[[c]]$post_P_0_imp-scenario_2[[c]]$data$P_1+scenario_2[[c]]$data$P_0)^2)
mse_P_SLM3=sapply(1:samples, function(c) 
  (SLM_3[[c]]$post_P_1_imp-SLM_3[[c]]$post_P_0_imp-scenario_3[[c]]$data$P_1+scenario_3[[c]]$data$P_0)^2)
mse_P_SLM4=sapply(1:samples, function(c) 
  (SLM_4[[c]]$post_P_1_imp-SLM_4[[c]]$post_P_0_imp-scenario_4[[c]]$data$P_1+scenario_4[[c]]$data$P_0)^2)
mse_P_SLM5=sapply(1:samples, function(c) 
  (SLM_5[[c]]$post_P_1_imp-SLM_5[[c]]$post_P_0_imp-scenario_5[[c]]$data$P_1+scenario_5[[c]]$data$P_0)^2)


###################################################################################

# --- bias and MSE for Y: outcome ---

# CASBAH
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

# Schartz Li Mealli's model
bias_Y_SLM1=sapply(1:samples, function(c) 
  SLM_1[[c]]$post_Y_1_imp-SLM_1[[c]]$post_Y_0_imp-scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)
bias_Y_SLM2=sapply(1:samples, function(c) 
  SLM_2[[c]]$post_Y_1_imp-SLM_2[[c]]$post_Y_0_imp-scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)
bias_Y_SLM3=sapply(1:samples, function(c) 
  SLM_3[[c]]$post_Y_1_imp-SLM_3[[c]]$post_Y_0_imp-scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)
bias_Y_SLM4=sapply(1:samples, function(c) 
  SLM_4[[c]]$post_Y_1_imp-SLM_4[[c]]$post_Y_0_imp-scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)
bias_Y_SLM5=sapply(1:samples, function(c) 
  SLM_5[[c]]$post_Y_1_imp-SLM_5[[c]]$post_Y_0_imp-scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)

mse_Y_SLM1=sapply(1:samples, function(c) 
  (SLM_1[[c]]$post_Y_1_imp-SLM_1[[c]]$post_Y_0_imp-scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)^2)
mse_Y_SLM2=sapply(1:samples, function(c) 
  (SLM_2[[c]]$post_Y_1_imp-SLM_2[[c]]$post_Y_0_imp-scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)^2)
mse_Y_SLM3=sapply(1:samples, function(c) 
  (SLM_3[[c]]$post_Y_1_imp-SLM_3[[c]]$post_Y_0_imp-scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)^2)
mse_Y_SLM4=sapply(1:samples, function(c) 
  (SLM_4[[c]]$post_Y_1_imp-SLM_4[[c]]$post_Y_0_imp-scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)^2)
mse_Y_SLM5=sapply(1:samples, function(c) 
  (SLM_5[[c]]$post_Y_1_imp-SLM_5[[c]]$post_Y_0_imp-scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)^2)

# LJD's model rho=0.25
bias_Y_LJD1=sapply(1:samples, function(c) 
  LJD_1[[c]]$Mu1-LJD_1[[c]]$Mu0-scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)
bias_Y_LJD2=sapply(1:samples, function(c) 
  LJD_2[[c]]$Mu1-LJD_2[[c]]$Mu0-scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)
bias_Y_LJD3=sapply(1:samples, function(c) 
  LJD_3[[c]]$Mu1-LJD_3[[c]]$Mu0-scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)
bias_Y_LJD4=sapply(1:samples, function(c) 
  LJD_4[[c]]$Mu1-LJD_4[[c]]$Mu0-scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)
bias_Y_LJD5=sapply(1:samples, function(c) 
  LJD_5[[c]]$Mu1-LJD_5[[c]]$Mu0-scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)

mse_Y_LJD1=sapply(1:samples, function(c) 
  (LJD_1[[c]]$Mu1-LJD_1[[c]]$Mu0-scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)^2)
mse_Y_LJD2=sapply(1:samples, function(c) 
  (LJD_2[[c]]$Mu1-LJD_2[[c]]$Mu0-scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)^2)
mse_Y_LJD3=sapply(1:samples, function(c) 
  (LJD_3[[c]]$Mu1-LJD_3[[c]]$Mu0-scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)^2)
mse_Y_LJD4=sapply(1:samples, function(c) 
  (LJD_4[[c]]$Mu1-LJD_4[[c]]$Mu0-scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)^2)
mse_Y_LJD5=sapply(1:samples, function(c) 
  (LJD_5[[c]]$Mu1-LJD_5[[c]]$Mu0-scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)^2)

#########################################################################
save.image("estimands.RData")


round(c(median(apply(bias_Y_LJD1,2,mean)),median(apply(bias_Y_LJD2,2,mean)),
        median(apply(bias_Y_LJD3,2,mean)),median(apply(bias_Y_LJD4,2,mean)),
        median(apply(bias_Y_LJD5,2,mean))),4)
round(c(IQR(apply(bias_Y_LJD1,2,mean)),IQR(apply(bias_Y_LJD2,2,mean)),
        IQR(apply(bias_Y_LJD3,2,mean)),IQR(apply(bias_Y_LJD4,2,mean)),
        IQR(apply(bias_Y_LJD5,2,mean))),4)


#########################################################################

#########################################################################
PCE_true_1 <- sapply(1:samples, function(i) 
  sort(sapply(unique(scenario_1[[i]]$clusters$S_strata),function(s)
    mean((scenario_1[[i]]$data$Y_1-scenario_1[[i]]$data$Y_0)[scenario_1[[i]]$clusters$S_strata==s]))))
PCE_true_2 <- sapply(1:samples, function(i) 
  sort(sapply(unique(scenario_2[[i]]$clusters$S_strata),function(s)
    mean((scenario_2[[i]]$data$Y_1-scenario_2[[i]]$data$Y_0)[scenario_2[[i]]$clusters$S_strata==s]))))
PCE_true_3 <- sapply(1:samples, function(i) 
  sort(sapply(unique(scenario_3[[i]]$clusters$S_strata),function(s)
    mean((scenario_3[[i]]$data$Y_1-scenario_3[[i]]$data$Y_0)[scenario_3[[i]]$clusters$S_strata==s]))))
PCE_true_4 <- sapply(1:samples, function(i) 
  sort(sapply(unique(scenario_4[[i]]$clusters$S_strata),function(s)
    mean((scenario_4[[i]]$data$Y_1-scenario_4[[i]]$data$Y_0)[scenario_4[[i]]$clusters$S_strata==s]))))
PCE_true_5 <- sapply(1:samples, function(i) 
  sort(sapply(unique(scenario_5[[i]]$clusters$S_strata),function(s)
    mean((scenario_5[[i]]$data$Y_1-scenario_5[[i]]$data$Y_0)[scenario_5[[i]]$clusters$S_strata==s]))))


PCE_CASDMM_1 <- sapply(1:samples, function(i) 
  sort(sapply(1:max(CASDMM_1[[i]]$S_strata_cluste),function(s)
    mean((CASDMM_1[[i]]$post_Y_1_imp-CASDMM_1[[i]]$post_Y_0_imp)[CASDMM_1[[i]]$S_strata_cluste==s]))))
PCE_CASDMM_2 <- sapply(1:samples, function(i) 
  sort(sapply(1:max(CASDMM_2[[i]]$S_strata_cluste),function(s)
    mean((CASDMM_2[[i]]$post_Y_1_imp-CASDMM_2[[i]]$post_Y_0_imp)[CASDMM_2[[i]]$S_strata_cluste==s]))))
PCE_CASDMM_3 <- sapply(1:samples, function(i) 
  sort(sapply(1:max(CASDMM_3[[i]]$S_strata_cluste),function(s)
    mean((CASDMM_3[[i]]$post_Y_1_imp-CASDMM_3[[i]]$post_Y_0_imp)[CASDMM_3[[i]]$S_strata_cluste==s]))))
PCE_CASDMM_4 <- sapply(1:samples, function(i) 
  sort(sapply(1:max(CASDMM_4[[i]]$S_strata_cluste),function(s)
    mean((CASDMM_4[[i]]$post_Y_1_imp-CASDMM_4[[i]]$post_Y_0_imp)[CASDMM_4[[i]]$S_strata_cluste==s]))))
PCE_CASDMM_5 <- sapply(1:samples, function(i) 
  sort(sapply(1:max(CASDMM_5[[i]]$S_strata_cluste),function(s)
    mean((CASDMM_5[[i]]$post_Y_1_imp-CASDMM_5[[i]]$post_Y_0_imp)[CASDMM_5[[i]]$S_strata_cluste==s]))))

PCE_LJD_1 <- sapply(1:samples, function(i) sort(LJD_1[[i]]$CE[1:3]))
PCE_LJD_2 <- sapply(1:samples, function(i) sort(LJD_2[[i]]$CE[1:3]))
PCE_LJD_3 <- sapply(1:samples, function(i) sort(LJD_3[[i]]$CE[1:3]))
PCE_LJD_4 <- sapply(1:samples, function(i) sort(LJD_4[[i]]$CE[1:3]))
PCE_LJD_5 <- sapply(1:samples, function(i) sort(LJD_5[[i]]$CE[1:3]))

round(c(apply(PCE_true_1,1,median),apply(PCE_true_2,1,median),
        apply(PCE_true_3,1,median),apply(PCE_true_4,1,median),
        apply(PCE_true_5,1,median)),1)

mean_our <-c(mean(sapply(PCE_CASDMM_1, `[[`, 1)),
             mean(sapply(PCE_CASDMM_1, function(x) if (length(x) > 1) x[[2]] else NA), na.rm = TRUE),
  mean(sapply(PCE_CASDMM_2, `[[`, 1)),mean(sapply(PCE_CASDMM_2, `[[`, 2)),
  mean(sapply(PCE_CASDMM_2, function(x) if (length(x) > 2) x[[3]] else NA), na.rm = TRUE),
  mean(sapply(PCE_CASDMM_3, `[[`, 1)),
  mean(sapply(PCE_CASDMM_3, function(x) if (length(x) > 1) x[[2]] else NA), na.rm = TRUE),
  mean(sapply(PCE_CASDMM_4, `[[`, 1)),mean(sapply(PCE_CASDMM_4, `[[`, 2)),
  mean(sapply(PCE_CASDMM_4, function(x) if (length(x) > 2) x[[3]] else NA), na.rm = TRUE),
  mean(sapply(PCE_CASDMM_5, `[[`, 1)),mean(sapply(PCE_CASDMM_5, `[[`, 2)),
  mean(sapply(PCE_CASDMM_5, function(x) if (length(x) > 2) x[[3]] else NA), na.rm = TRUE))
sd_our <-c(sd(sapply(PCE_CASDMM_1, `[[`, 1)),
           sd(sapply(PCE_CASDMM_1, function(x) if (length(x) > 1) x[[2]] else NA), na.rm = TRUE),
           sd(sapply(PCE_CASDMM_2, `[[`, 1)),sd(sapply(PCE_CASDMM_2, `[[`, 2)),
           sd(sapply(PCE_CASDMM_2, function(x) if (length(x) > 2) x[[3]] else NA), na.rm = TRUE),
           sd(sapply(PCE_CASDMM_3, `[[`, 1)),
           sd(sapply(PCE_CASDMM_3, function(x) if (length(x) > 1) x[[2]] else NA), na.rm = TRUE),
           sd(sapply(PCE_CASDMM_4, `[[`, 1)),sd(sapply(PCE_CASDMM_4, `[[`, 2)),
           sd(sapply(PCE_CASDMM_4, function(x) if (length(x) > 2) x[[3]] else NA), na.rm = TRUE),
           sd(sapply(PCE_CASDMM_5, `[[`, 1)),sd(sapply(PCE_CASDMM_5, `[[`, 2)),
           sd(sapply(PCE_CASDMM_5, function(x) if (length(x) > 2) x[[3]] else NA), na.rm = TRUE))

LJD_mean<-c(apply(PCE_LJD_1,1,mean),apply(PCE_LJD_2,1,mean),
            apply(PCE_LJD_3,1,mean),apply(PCE_LJD_4,1,mean),
            apply(PCE_LJD_5,1,mean))
LJD_sd<-c(apply(PCE_LJD_1,1,sd),apply(PCE_LJD_2,1,sd),
            apply(PCE_LJD_3,1,sd),apply(PCE_LJD_4,1,sd),
            apply(PCE_LJD_5,1,sd))

xtable(cbind(mean_our,sd_our))
xtable(cbind(LJD_mean,LJD_sd))
