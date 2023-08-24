###################################################################################
# ---        SHARED ATOMS MODEL         ---
# ---        P-model + Y-model          ---
# ---         analysis results          ---
###################################################################################

#load data
load("C:/Users/dafne/Desktop/DOTTORATO/PROGETTO PhD _2/P2/results_shared_atoms.RData")

#libraries
library(mcclust)

###################################################################################


# --- PARTITIONS ---

# S_all_cluster
clusters_all_1=lapply(1:samples,function(c) as.vector(table(model_1[[c]]$S_all_cluster)))
clusters_all_2=lapply(1:samples,function(c) as.vector(table(model_2[[c]]$S_all_cluster)))
clusters_all_3=lapply(1:samples,function(c) as.vector(table(model_3[[c]]$S_all_cluster)))
clusters_all_4=lapply(1:samples,function(c) as.vector(table(model_4[[c]]$S_all_cluster)))
clusters_all_5=lapply(1:samples,function(c) as.vector(table(model_5[[c]]$S_all_cluster)))

# S_all_cluster
clusters_strata_1=lapply(1:samples,function(c) as.vector(table(model_1[[c]]$S_strata_cluster)))
clusters_strata_2=lapply(1:samples,function(c) as.vector(table(model_2[[c]]$S_strata_cluster)))
clusters_strata_3=lapply(1:samples,function(c) as.vector(table(model_3[[c]]$S_strata_cluster)))
clusters_strata_4=lapply(1:samples,function(c) as.vector(table(model_4[[c]]$S_strata_cluster)))
clusters_strata_5=lapply(1:samples,function(c) as.vector(table(model_5[[c]]$S_strata_cluster)))

# number of clusters
n_cl_all_1=unlist(lapply(clusters_all_1,function(c) length(c)))
n_cl_all_2=unlist(lapply(clusters_all_2,function(c) length(c)))
n_cl_all_3=unlist(lapply(clusters_all_3,function(c) length(c)))
n_cl_all_4=unlist(lapply(clusters_all_4,function(c) length(c)))
n_cl_all_5=unlist(lapply(clusters_all_5,function(c) length(c)))

n_cl_strata_1=unlist(lapply(clusters_strata_1,function(c) length(c)))
n_cl_strata_2=unlist(lapply(clusters_strata_2,function(c) length(c)))
n_cl_strata_3=unlist(lapply(clusters_strata_3,function(c) length(c)))
n_cl_strata_4=unlist(lapply(clusters_strata_4,function(c) length(c)))
n_cl_strata_5=unlist(lapply(clusters_strata_5,function(c) length(c)))

table(n_cl_all_1)
table(n_cl_all_2)
table(n_cl_all_3)
table(n_cl_all_4)
table(n_cl_all_5)

table(n_cl_strata_1)
table(n_cl_strata_2)
table(n_cl_strata_3)
table(n_cl_strata_4)
table(n_cl_strata_5)

# --- RAD index ---

RAND_1_all<-sapply(1:samples, function(c)
  arandi(model_1[[c]]$S_all_cluster,
         scenario_1[[c]]$clusters_all[,1]*(scenario_1[[c]]$clusters_all[,2]+4)))
RAND_2_all<-sapply(1:samples, function(c)
  arandi(model_2[[c]]$S_all_cluster,
         scenario_2[[c]]$clusters_all[,1]*(scenario_2[[c]]$clusters_all[,2]+4)))
RAND_3_all<-sapply(1:samples, function(c)
  arandi(model_3[[c]]$S_all_cluster,
         scenario_3[[c]]$clusters_all[,1]*(scenario_3[[c]]$clusters_all[,2]+4)))
RAND_4_all<-sapply(1:samples, function(c)
  arandi(model_4[[c]]$S_all_cluster,
         scenario_4[[c]]$clusters_all[,1]*(scenario_4[[c]]$clusters_all[,2]+4)))
RAND_5_all<-sapply(1:samples, function(c)
  arandi(model_5[[c]]$S_all_cluster,
         scenario_5[[c]]$clusters_all[,1]*(scenario_5[[c]]$clusters_all[,2]+4)))

RAND_1_strata<-sapply(1:samples, function(c)
  arandi(model_1[[c]]$S_strata_cluster,
         1+1*(scenario_1[[c]]$clusters_all[,1]==scenario_1[[c]]$clusters_all[,2])+
           2*(scenario_1[[c]]$clusters_all[,1]>scenario_1[[c]]$clusters_all[,2])))
RAND_2_strata<-sapply(1:samples, function(c)
  arandi(model_2[[c]]$S_strata_cluster,
         1+1*(scenario_2[[c]]$clusters_all[,1]==scenario_2[[c]]$clusters_all[,2])+
           2*(scenario_2[[c]]$clusters_all[,1]<scenario_2[[c]]$clusters_all[,2])))
RAND_3_strata<-sapply(1:samples, function(c)
  arandi(model_3[[c]]$S_strata_cluster,
         1+1*(scenario_3[[c]]$clusters_all[,1]==scenario_3[[c]]$clusters_all[,2])+
           2*(scenario_3[[c]]$clusters_all[,1]>scenario_3[[c]]$clusters_all[,2])))
RAND_4_strata<-sapply(1:samples, function(c)
  arandi(model_4[[c]]$S_strata_cluster,
         1+1*(scenario_4[[c]]$clusters_all[,1]==scenario_4[[c]]$clusters_all[,2])+
           2*(scenario_4[[c]]$clusters_all[,1]>scenario_4[[c]]$clusters_all[,2])))
RAND_5_strata<-sapply(1:samples, function(c)
  arandi(model_5[[c]]$S_strata_cluster,
         1+1*(scenario_5[[c]]$clusters_all[,1]==scenario_5[[c]]$clusters_all[,2])+
           2*(scenario_5[[c]]$clusters_all[,1]<scenario_5[[c]]$clusters_all[,2])))

rand_all<-matrix(c(mean(RAND_1_all),sd(RAND_1_all),
                 mean(RAND_2_all),sd(RAND_2_all),
                 mean(RAND_3_all),sd(RAND_3_all),
                 mean(RAND_4_all),sd(RAND_4_all),
                 mean(RAND_5_all),sd(RAND_5_all)),ncol=5)
row.names(rand_all)<-c("mean","sd")
rand_strata<-matrix(c(mean(RAND_1_strata),sd(RAND_1_strata),
                   mean(RAND_2_strata),sd(RAND_2_strata),
                   mean(RAND_3_strata),sd(RAND_3_strata),
                   mean(RAND_4_strata),sd(RAND_4_strata),
                   mean(RAND_5_strata),sd(RAND_5_strata)),ncol=5)
row.names(rand_strata)<-c("mean","sd")

rand_all
rand_strata

###################################################################################

# --- P estimation ---

diff_P_1=sapply(1:samples, function(c) 
  model_1[[c]]$post_P_1_imp-model_1[[c]]$post_P_0_imp-
    scenario_1[[c]]$data$P_1+scenario_1[[c]]$data$P_0)
diff_P_2=sapply(1:samples, function(c) 
  model_2[[c]]$post_P_1_imp-model_2[[c]]$post_P_0_imp-
    scenario_2[[c]]$data$P_1+scenario_2[[c]]$data$P_0)
diff_P_3=sapply(1:samples, function(c) 
  model_3[[c]]$post_P_1_imp-model_3[[c]]$post_P_0_imp-
    scenario_3[[c]]$data$P_1+scenario_3[[c]]$data$P_0)
diff_P_4=sapply(1:samples, function(c) 
  model_4[[c]]$post_P_1_imp-model_4[[c]]$post_P_0_imp-
    scenario_4[[c]]$data$P_1+scenario_4[[c]]$data$P_0)
diff_P_5=sapply(1:samples, function(c) 
  model_5[[c]]$post_P_1_imp-model_5[[c]]$post_P_0_imp-
    scenario_5[[c]]$data$P_1+scenario_5[[c]]$data$P_0)

# bias
par(mfrow=c(1,5))
boxplot(apply(diff_P_1, 2, mean), #ylim=c(-0.8,0.1), 
        main="setting 1", ylab="Bias", xlab="P(1)-P(0)")
abline(h=0, col="red")
boxplot(apply(diff_P_2, 2, mean), #ylim=c(-0.8,0.1), 
        main="setting 2", ylab="Bias", xlab="P(1)-P(0)")
abline(h=0, col="red")
boxplot(apply(diff_P_3, 2, mean), #ylim=c(-0.8,0.1), 
        main="setting 3", ylab="Bias", xlab="P(1)-P(0)")
abline(h=0, col="red")
boxplot(apply(diff_P_4, 2, mean), #ylim=c(-0.8,0.1), 
        main="setting 4", ylab="Bias", xlab="P(1)-P(0)")
abline(h=0, col="red")
boxplot(apply(diff_P_5, 2, mean), #ylim=c(-0.8,0.1), 
        main="setting 5", ylab="Bias", xlab="P(1)-P(0)")
abline(h=0, col="red")

# MSE
par(mfrow=c(1,5))
boxplot(sqrt(apply(diff_P_1^2, 2, mean)), ylim=c(0.2,1.2), 
        main="setting 1", ylab="MSE", xlab="P(1)-P(0)")
boxplot(sqrt(apply(diff_P_2^2, 2, mean)), ylim=c(0.2,1.2),
        main="setting 2", ylab="MSE", xlab="P(1)-P(0)")
boxplot(sqrt(apply(diff_P_3^2, 2, mean)), ylim=c(0.2,1.2),
        main="setting 3", ylab="MSE", xlab="P(1)-P(0)")
boxplot(sqrt(apply(diff_P_4^2, 2, mean)), ylim=c(0.2,1.2),
        main="setting 4", ylab="MSE", xlab="P(1)-P(0)")
boxplot(sqrt(apply(diff_P_5^2, 2, mean)), ylim=c(0.2,1.2),
        main="setting 5", ylab="MSE", xlab="P(1)-P(0)")

median_bias_P<-c(median(apply(diff_P_1, 2, mean)),
                 median(apply(diff_P_2, 2, mean)),
                 median(apply(diff_P_3, 2, mean)),
                 median(apply(diff_P_4, 2, mean)),
                 median(apply(diff_P_5, 2, mean)))
range_bias_P<-c(IQR(apply(diff_P_1, 2, mean)),
                IQR(apply(diff_P_2, 2, mean)),
                IQR(apply(diff_P_3, 2, mean)),
                IQR(apply(diff_P_4, 2, mean)),
                IQR(apply(diff_P_5, 2, mean)))

###################################################################################

# --- P cond CLUSTER ---

identify_cluster_P<-function(part,P_1,P_0,cl_est){
  n_cl=table(part)
  n_cl=n_cl[order(n_cl,decreasing=TRUE)[1:cl_est]]
  est=sapply(as.numeric(names(n_cl)),function(l)
    mean(P_1[part==l]-P_0[part==l]))
  return(est[order(est)])
}

P_cl_all_1=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_1[[c]]$S_all_cluster,
                     P_1=model_1[[c]]$post_P_1_imp,
                     P_0=model_1[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_all_1)))
P_cl_all_2=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_2[[c]]$S_all_cluster,
                     P_1=model_2[[c]]$post_P_1_imp,
                     P_0=model_2[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_all_2)))
P_cl_all_3=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_3[[c]]$S_all_cluster,
                     P_1=model_3[[c]]$post_P_1_imp,
                     P_0=model_3[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_all_3)))
P_cl_all_4=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_4[[c]]$S_all_cluster,
                     P_1=model_4[[c]]$post_P_1_imp,
                     P_0=model_4[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_all_4)))
P_cl_all_5=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_5[[c]]$S_all_cluster,
                     P_1=model_5[[c]]$post_P_1_imp,
                     P_0=model_5[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_all_5)))

P_cl_strata_1=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_1[[c]]$S_strata_cluster,
                     P_1=model_1[[c]]$post_P_1_imp,
                     P_0=model_1[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_strata_1)))
P_cl_strata_2=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_2[[c]]$S_strata_cluster,
                     P_1=model_2[[c]]$post_P_1_imp,
                     P_0=model_2[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_strata_2)))
P_cl_strata_3=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_3[[c]]$S_strata_cluster,
                     P_1=model_3[[c]]$post_P_1_imp,
                     P_0=model_3[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_strata_3)))
P_cl_strata_4=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_4[[c]]$S_strata_cluster,
                     P_1=model_4[[c]]$post_P_1_imp,
                     P_0=model_4[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_strata_4)))
P_cl_strata_5=sapply(1:samples, function(c) 
  identify_cluster_P(part=model_5[[c]]$S_strata_cluster,
                     P_1=model_5[[c]]$post_P_1_imp,
                     P_0=model_5[[c]]$post_P_0_imp,
                     cl_est=median(n_cl_strata_5)))

# conditional causal effect 
par(mfrow=c(1,5))
boxplot(c(P_cl_all_1)~rep(1:(median(n_cl_all_1)),samples), 
        main="setting 1", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=scenario_1[[1]]$par_P[scenario_1[[1]]$par_P[,4],1]-
         scenario_1[[1]]$par_P[scenario_1[[1]]$par_P[,3],1],
       col="red")
boxplot(c(P_cl_all_2)~rep(1:(median(n_cl_all_2)),samples), 
        main="setting 2", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=scenario_2[[1]]$par_P[scenario_2[[1]]$par_P[,4],1]-
         scenario_2[[1]]$par_P[scenario_2[[1]]$par_P[,3],1],
       col="red")
boxplot(c(P_cl_all_3)~rep(1:(median(n_cl_all_3)),samples),  
        main="setting 3", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=scenario_3[[1]]$par_P[scenario_3[[1]]$par_P[,4],1]-
         scenario_3[[1]]$par_P[scenario_3[[1]]$par_P[,3],1],
       col="red")
boxplot(c(P_cl_all_4)~rep(1:(median(n_cl_all_4)),samples), 
        main="setting 4", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=scenario_4[[1]]$par_P[scenario_4[[1]]$par_P[,4],1]-
         scenario_4[[1]]$par_P[scenario_4[[1]]$par_P[,3],1],
       col="red")
boxplot(c(P_cl_all_5)~rep(1:(median(n_cl_all_5)),samples), 
        main="setting 5", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=scenario_5[[1]]$par_P[scenario_5[[1]]$par_P[,4],1]-
         scenario_5[[1]]$par_P[scenario_5[[1]]$par_P[,3],1],
       col="red")

par(mfrow=c(1,5))
boxplot(c(P_cl_strata_1)~rep(1:(median(n_cl_strata_1)),samples), 
        main="setting 1", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=c(0,1),col="red")
boxplot(c(P_cl_strata_2)~rep(1:(median(n_cl_strata_2)),samples), 
        main="setting 2", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=c(0,1),col="red")
boxplot(c(P_cl_strata_3)~rep(1:(median(n_cl_strata_3)),samples), 
        main="setting 3", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=c(-1,0,1),col="red")
boxplot(c(P_cl_strata_4)~rep(1:(median(n_cl_strata_4)),samples),  
        main="setting 4", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=c(0,0.5),col="red")
boxplot(c(P_cl_strata_5)~rep(1:(median(n_cl_strata_5)),samples), 
        main="setting 5", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=c(-0.5,0,0.5),col="red")

###################################################################################

# --- Y estimation ---

diff_Y_1=sapply(1:samples, function(c) 
  model_1[[c]]$post_Y_1_imp-model_1[[c]]$post_Y_0_imp-
    scenario_1[[c]]$data$Y_1+scenario_1[[c]]$data$Y_0)
diff_Y_2=sapply(1:samples, function(c) 
  model_2[[c]]$post_Y_1_imp-model_2[[c]]$post_Y_0_imp-
    scenario_2[[c]]$data$Y_1+scenario_2[[c]]$data$Y_0)
diff_Y_3=sapply(1:samples, function(c) 
  model_3[[c]]$post_Y_1_imp-model_3[[c]]$post_Y_0_imp-
    scenario_3[[c]]$data$Y_1+scenario_3[[c]]$data$Y_0)
diff_Y_4=sapply(1:samples, function(c) 
  model_4[[c]]$post_Y_1_imp-model_4[[c]]$post_Y_0_imp-
    scenario_4[[c]]$data$Y_1+scenario_4[[c]]$data$Y_0)
diff_Y_5=sapply(1:samples, function(c) 
  model_5[[c]]$post_Y_1_imp-model_5[[c]]$post_Y_0_imp-
    scenario_5[[c]]$data$Y_1+scenario_5[[c]]$data$Y_0)

# bias
par(mfrow=c(1,5))
boxplot(apply(diff_Y_1, 2, mean), ylim=c(-1.3,0.5), 
        main="setting 1", ylab="Bias", xlab="Y(1)-Y(0)")
abline(h=0, col="red")
boxplot(apply(diff_Y_2, 2, mean), ylim=c(-1.3,0.5), 
        main="setting 2", ylab="Bias", xlab="Y(1)-Y(0)")
abline(h=0, col="red")
boxplot(apply(diff_Y_3, 2, mean), ylim=c(-1.3,0.5), 
        main="setting 3", ylab="Bias", xlab="Y(1)-Y(0)")
abline(h=0, col="red")
boxplot(apply(diff_Y_4, 2, mean), ylim=c(-1.3,0.5), 
        main="setting 4", ylab="Bias", xlab="Y(1)-Y(0)")
abline(h=0, col="red")
boxplot(apply(diff_Y_5, 2, mean), ylim=c(-1.3,0.5), 
        main="setting 5", ylab="Bias", xlab="Y(1)-Y(0)")
abline(h=0, col="red")

# MSE
par(mfrow=c(1,5))
boxplot(sqrt(apply(diff_Y_1^2, 2, mean)), ylim=c(0.7,4.2), 
        main="setting 1", ylab="MSE", xlab="Y(1)-Y(0)")
boxplot(sqrt(apply(diff_Y_2^2, 2, mean)), ylim=c(0.7,4.2), 
        main="setting 2", ylab="MSE", xlab="Y(1)-Y(0)")
boxplot(sqrt(apply(diff_Y_3^2, 2, mean)), ylim=c(0.7,4.2), 
        main="setting 3", ylab="MSE", xlab="Y(1)-Y(0)")
boxplot(sqrt(apply(diff_Y_4^2, 2, mean)), ylim=c(0.7,4.), 
        main="setting 4", ylab="MSE", xlab="Y(1)-Y(0)")
boxplot(sqrt(apply(diff_Y_5^2, 2, mean)), ylim=c(0.7,4.), 
        main="setting 5", ylab="MSE", xlab="Y(1)-Y(0)")


median_bias_Y<-c(median(apply(diff_Y_1, 2, mean)),
                 median(apply(diff_Y_2, 2, mean)),
                 median(apply(diff_Y_3, 2, mean)),
                 median(apply(diff_Y_4, 2, mean)),
                 median(apply(diff_Y_5, 2, mean)))
range_bias_Y<-c(IQR(apply(diff_Y_1, 2, mean)),
                IQR(apply(diff_Y_2, 2, mean)),
                IQR(apply(diff_Y_3, 2, mean)),
                IQR(apply(diff_Y_4, 2, mean)),
                IQR(apply(diff_Y_5, 2, mean)))

###################################################################################

# --- Y cond CLUSTER ---

identify_cluster<-function(part,Y_0,Y_1, cl_est){
  n_cl=table(part)
  n_cl=n_cl[order(n_cl,decreasing=TRUE)[1:cl_est]]
  est=sapply(as.numeric(names(n_cl)),function(l)
    mean(Y_1[part==l]-Y_0[part==l]))
  return(est[order(est)])
}

Y_cl_all_1=sapply(1:samples, function(c) 
  identify_cluster(part=model_1[[c]]$S_all_cluster,
                   Y_0=model_1[[c]]$post_Y_0_imp,
                   Y_1=model_1[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_all_1)))
Y_cl_all_2=sapply(1:samples, function(c) 
  identify_cluster(part=model_2[[c]]$S_all_cluster,
                   Y_0=model_2[[c]]$post_Y_0_imp,
                   Y_1=model_2[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_all_2)))
Y_cl_all_3=sapply(1:samples, function(c) 
  identify_cluster(part=model_3[[c]]$S_all_cluster,
                   Y_0=model_3[[c]]$post_Y_0_imp,
                   Y_1=model_3[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_all_3)))
Y_cl_all_4=sapply(1:samples, function(c) 
  identify_cluster(part=model_4[[c]]$S_all_cluster,
                   Y_0=model_4[[c]]$post_Y_0_imp,
                   Y_1=model_4[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_all_4)))
Y_cl_all_5=sapply(1:samples, function(c) 
  identify_cluster(part=model_5[[c]]$S_all_cluster,
                   Y_0=model_5[[c]]$post_Y_0_imp,
                   Y_1=model_5[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_all_5)))

Y_cl_strata_1=sapply(1:samples, function(c) 
  identify_cluster(part=model_1[[c]]$S_strata_cluster,
                   Y_0=model_1[[c]]$post_Y_0_imp,
                   Y_1=model_1[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_strata_1)))
Y_cl_strata_2=sapply(1:samples, function(c) 
  identify_cluster(part=model_2[[c]]$S_strata_cluster,
                   Y_0=model_2[[c]]$post_Y_0_imp,
                   Y_1=model_2[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_strata_2)))
Y_cl_strata_3=sapply(1:samples, function(c) 
  identify_cluster(part=model_3[[c]]$S_strata_cluster,
                   Y_0=model_3[[c]]$post_Y_0_imp,
                   Y_1=model_3[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_strata_3)))
Y_cl_strata_4=sapply(1:samples, function(c) 
  identify_cluster(part=model_4[[c]]$S_strata_cluster,
                   Y_0=model_4[[c]]$post_Y_0_imp,
                   Y_1=model_4[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_strata_4)))
Y_cl_strata_5=sapply(1:samples, function(c) 
  identify_cluster(part=model_5[[c]]$S_strata_cluster,
                   Y_0=model_5[[c]]$post_Y_0_imp,
                   Y_1=model_5[[c]]$post_Y_1_imp,
                   cl_est=median(n_cl_strata_5)))

# conditional causal effect 
par(mfrow=c(1,5))
boxplot(c(Y_cl_all_1)~rep(1:(median(n_cl_all_1)),samples),
        main="setting 1", ylab=" ", xlab="Y(1)-Y(0) | cluster")
#abline(h=c(0.5,9),col="red")
boxplot(c(Y_cl_all_2)~rep(1:(median(n_cl_all_2)),samples),
        main="setting 2", ylab=" ", xlab="Y(1)-Y(0) | cluster")
#abline(h=c(4,9,10.5),col="red")
boxplot(c(Y_cl_all_3)~rep(1:(median(n_cl_all_3)),samples),
        main="setting 3", ylab=" ", xlab="Y(1)-Y(0) | cluster")
#abline(h=c(4,6.5,6.875),col="red")
boxplot(c(Y_cl_all_4)~rep(1:(median(n_cl_all_4)),samples), 
        main="setting 4", ylab=" ", xlab="Y(1)-Y(0) | cluster")
#abline(h=c(4,6.5,6.875),col="red")
boxplot(c(Y_cl_all_5)~rep(1:(median(n_cl_all_5)),samples),
        main="setting 5", ylab=" ", xlab="Y(1)-Y(0) | cluster")
#abline(h=c(4,6.5,6.875),col="red")

par(mfrow=c(1,5))
boxplot(c(Y_cl_strata_1)~rep(1:(median(n_cl_strata_1)),samples),
        main="setting 1", ylab=" ", xlab="Y(1)-Y(0) | strata")
#abline(h=c(0.5,9),col="red")
boxplot(c(Y_cl_strata_2)~rep(1:(median(n_cl_strata_2)),samples), 
        main="setting 2", ylab=" ", xlab="Y(1)-Y(0) | strata")
#abline(h=c(6.5,10.5),col="red")
boxplot(c(Y_cl_strata_3)~rep(1:(median(n_cl_strata_3)),samples),
        main="setting 3", ylab=" ", xlab="Y(1)-Y(0) | strata")
#abline(h=c(6.5,10.5),col="red")
boxplot(c(Y_cl_strata_4)~rep(1:(median(n_cl_strata_4)),samples),
        main="setting 4", ylab=" ", xlab="Y(1)-Y(0) | strata")
#abline(h=c(5.25,6.875),col="red")
boxplot(c(Y_cl_strata_5)~rep(1:(median(n_cl_strata_5)),samples),
        main="setting 5", ylab=" ", xlab="Y(1)-Y(0) | strata")
#abline(h=c(5.25,6.875),col="red")


###################################################################################


