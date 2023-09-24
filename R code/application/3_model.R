##############################################################
# ---       estimation model     ----
##############################################################


#load matched data
load("dataset_matched.RData")

#library
library(mvtnorm)
library(CholWishart)
library(truncnorm)
library(invgamma)
library(BNPmix)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(fmsb)

###################################################################

# iterations
R=4000
R_burnin=3000

# max number clusters
n_cluster=10

source("Gibbs.R")

###################################################################

# prepar dataset --- common part
matrix_X=dataset_matched[,c("PctBlack","PctHisp","PctHighSchool",
                            "PctUrban","PctFemale","PctPoor")]
matrix_X=as.matrix(matrix_X)

T_var=dataset_matched$a
P_obs=dataset_matched$pm_diff
Y_obs=dataset_matched$all_causes_ADJ

system.time(results_CASBAH <- Gibbs_CASDMM(1,matrix_X,T_var,P_obs,Y_obs))
system.time(results_CASBAH_cov <- Gibbs_CASDMM_cov(1,matrix_X,T_var,P_obs,Y_obs))
system.time(results_CASBAH_cov_constrains <- Gibbs_CASDMM_cov_constrains(1,matrix_X,T_var,P_obs,Y_obs))


results=results_CASBAH
results=results_CASBAH_cov
results=results_CASBAH_cov_constrains

save.image(results,"results_application.RData")


###################################################################
###################################################################

# Provare:

# - tutte le covariate
# - diversi death causes

###################################################################

# Y-model with ALL covariates
matrix_REGR=dataset_matched[,c("PctBlack","PctHisp","PctHighSchool",
                               "PctUrban","PctFemale","PctPoor","MedianHHInc",
                               "PctOccupied","PctMovedIn5","MedianHValue",
                               "smokerate2000","avgdewpt","avgtemp",
                               "avgrelhum","Population")]
matrix_REGR=as.matrix(matrix_REGR)


###################################################################

# ----  BCF  -------
p.score <- glm(T_var ~ matrix_X,
               family = binomial,
               data = as.data.frame(cbind(T_var, matrix_X)))
pihat <- predict(p.score, as.data.frame(matrix_X))

# Perform the Bayesian Causal Forest for the Proportion of Compliers (pic)
bcf_tau <- bcf(P_obs, T_var, matrix_X, matrix_X, pihat, nburn=1000, nsim=1000)
bcf_tau$tau

plot(density(bcf_tau$tau))
