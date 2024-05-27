##############################################################
# ---       estimation model     ----
##############################################################

# loading Gibbs sampler function
source("Gibbs.R")

# load dataset (no matching):
load("Dataset_merged.RData")

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
R=6000 
R_burnin=3500   

# max number clusters
n_cluster=12  

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

