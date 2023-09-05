table(model_data_diff$S_all_cluster[,1])
table(model_data_diff$S_all_cluster[,2])
table(model_data_diff$S_strata_cluster[,1])
table(model_data_diff$S_strata_cluster[,2])

table(model_data_fu$S_all_cluster[,1])
table(model_data_fu$S_all_cluster[,2])
table(model_data_fu$S_strata_cluster[,1])
table(model_data_fu$S_strata_cluster[,2])

par(mfrow=c(2,3))
for (i in 1:n_cluster){
  plot(model_data_diff$post_eta[i,],type="l")
  abline(h=mean(model_data_diff$post_eta[i,]), col="red")
  legend("topright",legend=round(mean(model_data_diff$post_eta[i,]),3))
}
for (i in 1:n_cluster){
  plot(model_data_diff$post_var[i,],type="l")
  abline(h=mean(model_data_diff$post_var[i,]), col="red")
  legend("topright",legend=round(mean(model_data_diff$post_var[i,]),3))
}


par(mfrow=c(3,3))
for (i in 1:n_cluster){
  plot(model_data_fu$post_eta[i,],type="l")
  abline(h=mean(model_data_fu$post_eta[i,]), col="red")
  legend("topright",legend=round(mean(model_data_fu$post_eta[i,]),3))
}
for (i in 1:n_cluster){
  plot(model_data_fu$post_var[i,],type="l")
  abline(h=mean(model_data_fu$post_var[i,]), col="red")
  legend("topright",legend=round(mean(model_data_fu$post_var[i,]),3))
}

#par(mfrow=c(2,1))
for (i in 1:2){
  plot(model_data$post_lambda[i,],type="l")
  abline(h=mean(model_data$post_lambda[i,]), col="red")
  legend("topright",legend=round(mean(model_data$post_lambda[i,]),3))
}


apply(model_data$post_xi_0, 1, median)
apply(model_data$post_xi_1, 1, median)

##############################################################
# ---       observed PM     ----
##############################################################

par(mfrow=c(2,1))
hist(dataset_matched$pm_diff[dataset_matched$a==1],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="PM fu", main="treated")
hist(dataset_matched$pm_diff[dataset_matched$a==0],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="PM fu", main="control")


##############################################################
# ---       imputed PM     ----
##############################################################

par(mfrow=c(1,1))
P_0=apply(model_data_diff$post_P_0,1,median)
P_1=apply(model_data_diff$post_P_1,1,median)
par(mfrow=c(1,1))
hist(P_1-P_0, main="P(1)-P(0)", xlab="P(1)-P(0)", nclass=50)

par(mfrow=c(2,1))
hist(P_1,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="PM diff", main="treated")
hist(P_0,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="PM diff", main="control")


##############################################################

# --- boxplot ---

P_clusters=lapply(1:n_cluster, function(l)
  P_1[model_data_diff$S_all_cluster[,2]==l]-
    P_0[model_data_diff$S_all_cluster[,2]==l])

seq_cluster=lapply(1:n_cluster, function(l) 
  rep(l,sum(model_data_diff$S_all_cluster[,2]==l)))

P_strata=lapply(1:n_cluster, function(l)
  P_1[model_data_diff$S_strata_cluster[,2]==l]-
    P_0[model_data_diff$S_strata_cluster[,2]==l])

seq_strata=lapply(1:n_cluster, function(l) 
  rep(l,sum(model_data_diff$S_strata_cluster[,2]==l)))


par(mfrow=c(1,1))
boxplot(unlist(P_clusters)~unlist(seq_cluster), 
        main=" ", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")

par(mfrow=c(1,1))
boxplot(unlist(P_strata)~unlist(seq_strata), 
        main=" ", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")



Y_0=apply(model_data_diff$post_Y_0_imp,1,median)
Y_1=apply(model_data_diff$post_Y_1_imp,1,median)
par(mfrow=c(1,1))
hist(Y_1-Y_0, nclass=50)

##############################################################

# --- Y ---

par(mfrow=c(3,1))
hist((Y_1-Y_0)[model_data_diff$S_strata_cluster[,2]==1], 
     nclass=50, xlim=c(min(Y_1-Y_0),max(Y_1-Y_0)))
hist((Y_1-Y_0)[model_data_diff$S_strata_cluster[,2]==2], 
     nclass=50, xlim=c(min(Y_1-Y_0),max(Y_1-Y_0)))
hist((Y_1-Y_0)[model_data_diff$S_strata_cluster[,2]==3], 
     nclass=50, xlim=c(min(Y_1-Y_0),max(Y_1-Y_0)))
