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
R=6000  #7000
R_burnin=3500   

# max number clusters
n_cluster=12  #8

source("Gibbs.R")

comparison_P<-sapply(4:18, function(i) 
  (c(mean(dataset_matched$pm_diff[dataset_matched[,i]==0]),
    mean(dataset_matched$pm_diff[dataset_matched[,i]==1]))-
     mean(dataset_matched$pm_diff))/sd(dataset_matched$pm_diff))

comparison_Y<-sapply(4:18, function(i) 
  (c(mean(dataset_matched$all_causes_ADJ[dataset_matched[,i]==0]),
  mean(dataset_matched$all_causes_ADJ[dataset_matched[,i]==1]))-
    mean(dataset_matched$all_causes_ADJ))/sd(dataset_matched$all_causes_ADJ))

colnames(comparison_P)<-colnames(dataset_matched[,4:18])
colnames(comparison_Y)<-colnames(dataset_matched[,4:18])

round(cbind(t(comparison_P),t(comparison_Y)),2)

###################################################################

# prepar dataset --- common part
matrix_X=dataset_matched[,c("PctBlack","PctHisp","PctHighSchool",
                            "PctUrban","PctFemale","PctPoor")]

matrix_X=sapply(1:(dim(matrix_X)[2]),function(c) 
  (matrix_X[,c]-mean(matrix_X[,c]))/sd(matrix_X[,c]))
matrix_X=as.matrix(matrix_X)

matrix_COV=dataset_matched[,4:18]
#matrix_COV=sapply(1:(dim(matrix_COV)[2]),function(c) 
#  (matrix_COV[,c]-mean(matrix_COV[,c]))/sd(matrix_COV[,c]))
matrix_COV=as.matrix(matrix_COV)

#dataset_matched[,c("PctBlack","PctHisp","PctHighSchool","PctUrban","PctFemale","PctPoor","MedianHHInc",
#                   "PctOccupied","PctMovedIn5","MedianHValue","smokerate2000","avgdewpt","avgtemp", "avgrelhum","Population")]


#matrix_X=apply(matrix_X,2,scale)

T_var=dataset_matched$a
P_obs=dataset_matched$pm_diff
Y_obs=dataset_matched$all_causes_ADJ

system.time(results_CASBAH <- Gibbs_CASDMM(1,matrix_X,T_var,P_obs,Y_obs))
system.time(results_CASBAH_cov <- Gibbs_CASDMM_cov(1,matrix_X,T_var,P_obs,Y_obs))
system.time(results_CASBAH_cov_constrains <- Gibbs_CASDMM_cov_constrains(1,matrix_X,T_var,P_obs,Y_obs))
system.time(results_CASBAH_P0_constrains <- Gibbs_CASDMM_cov_P0(1,matrix_X,T_var,P_obs,Y_obs))
system.time(results_CASBAH_counf <- Gibbs_CASDMM_cov_conf(1,matrix_X,T_var,P_obs,Y_obs,matrix_COV))
system.time(results_CASBAH_counf_P0 <- Gibbs_CASDMM_cov_conf_P0(1,matrix_X,T_var,P_obs,Y_obs,matrix_COV))


save(results,file="results_application.RData")


###################################################################
###################################################################

# NO matched data

# load dataset:
load("Dataset_merged.RData")

data_temp <- data_merged %>% 
  select(a,                                          #treatment
         pm_diff,                                    # post-treatment
         all_causes_ADJ,                             #outcome
         PctUrban, PctBlack, PctHisp, PctHighSchool,
         MedianHHInc, PctPoor, PctFemale,
         PctOccupied, PctMovedIn5, MedianHValue,
         smokerate2000, avgdewpt, avgtemp, avgrelhum, Population,
         FIPS)
data_temp_rnames <- row.names(data_temp)

#dicotomize the confounders
data_temp$PctUrban=1*(data_temp$PctUrban>median(data_temp$PctUrban))
data_temp$PctBlack=1*(data_temp$PctBlack>median(data_temp$PctBlack))
data_temp$PctHisp=1*(data_temp$PctHisp>median(data_temp$PctHisp))
data_temp$PctHighSchool=1*(data_temp$PctHighSchool>median(data_temp$PctHighSchool))
data_temp$MedianHHInc=1*(data_temp$MedianHHInc>median(data_temp$MedianHHInc))
data_temp$PctPoor=1*(data_temp$PctPoor>median(data_temp$PctPoor))
data_temp$PctFemale=1*(data_temp$PctFemale>median(data_temp$PctFemale))
data_temp$PctOccupied=1*(data_temp$PctOccupied>median(data_temp$PctOccupied))
data_temp$PctMovedIn5=1*(data_temp$PctMovedIn5>median(data_temp$PctMovedIn5))
data_temp$MedianHValue=1*(data_temp$MedianHValue>median(data_temp$MedianHValue))
data_temp$smokerate2000=1*(data_temp$smokerate2000>median(data_temp$smokerate2000))
data_temp$avgdewpt=1*(data_temp$avgdewpt>median(data_temp$avgdewpt))
data_temp$avgtemp=1*(data_temp$avgtemp>median(data_temp$avgtemp))
data_temp$avgrelhum=1*(data_temp$avgrelhum>median(data_temp$avgrelhum))
data_temp$Population=1*(data_temp$Population>median(data_temp$Population))

dataset_all<-data_temp

matrix_X=dataset_all[,c("PctBlack","PctHisp","PctHighSchool","PctUrban","PctFemale","PctPoor")]

matrix_X=sapply(1:(dim(matrix_X)[2]),function(c) 
  (matrix_X[,c]-mean(matrix_X[,c]))/sd(matrix_X[,c]))
matrix_X=as.matrix(matrix_X)

matrix_COV=dataset_all[,4:18]
matrix_COV=sapply(1:(dim(matrix_COV)[2]),function(c)  
  (matrix_COV[,c]-mean(matrix_COV[,c]))/sd(matrix_COV[,c]))
matrix_COV=as.matrix(matrix_COV)

T_var=dataset_all$a
P_obs=dataset_all$pm_diff
Y_obs=dataset_all$all_causes_ADJ

system.time(results_CASBAH_P_all_Y_all <- Gibbs_CASDMM_cov(1,matrix_COV,T_var,P_obs,Y_obs))
system.time(results_CASBAH_P_all_Y_all_constrainsP <- Gibbs_CASDMM_cov_P0(1,matrix_COV,T_var,P_obs,Y_obs))

save(results,file="results_application.RData")


###################################################################
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
