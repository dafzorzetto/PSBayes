##############################################################
# ---       analysis regults     ----
##############################################################

#load results
load("C:/Users/dafne/Desktop/DOTTORATO/PROGETTO PhD _2/data/results_model_10vs12.RData")

##############################################################

apply(model_data_REGweights_10$S_all_cluster,2,table)
apply(model_data_REGweights_10$S_strata_cluster,2,table)
apply(model_data_REGweights_10$S_mix_cluster,2,table)

apply(model_data_REGweights_12$S_all_cluster,2,table)
apply(model_data_REGweights_12$S_strata_cluster,2,table)
apply(model_data_REGweights_12$S_mix_cluster,2,table)


##############################################################
# ---       observed PM     ----
##############################################################

par(mfrow=c(2,1))
hist(dataset_matched$pm_diff[dataset_matched$a==0],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="diff PM", main="control", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
hist(dataset_matched$pm_diff[dataset_matched$a==1],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="diff PM", main="treated", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")

##############################################################
# ---       imputed PM     ----
##############################################################

# same covariates of the weights

P_0_10=apply(model_data_REGweights_10$post_P_0,1,median)
P_1_10=apply(model_data_REGweights_10$post_P_1,1,median)

P_0_12=apply(model_data_REGweights_12$post_P_0,1,median)
P_1_12=apply(model_data_REGweights_12$post_P_1,1,median)

##############################################################
# ---        all models        ---
# ---    distributions of P    ---
##############################################################

par(mfrow=c(3,2))
hist(dataset_matched$pm_diff[dataset_matched$a==0],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="observed", main="control", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
hist(dataset_matched$pm_diff[dataset_matched$a==1],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="observed PM", main="treated", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")

hist(P_0_10,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="6 cov - 10", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
abline(v=mean(P_0_10), col="green")
hist(P_1_10,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="6 cov - 10", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")
abline(v=mean(P_1_10), col="green")

hist(P_0_12,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="6 cov - 12", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
abline(v=mean(P_0_12), col="blue")
hist(P_1_12,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="6 cov - 12", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")
abline(v=mean(P_1_12), col="blue")

##############################################################

# indirect effect: P(1) - P(0)

par(mfrow=c(2,1))
hist(P_1_10-P_0_10,nclass=50,xlim=c(-3,2),
     xlab="6 cov - 10", main="", freq=FALSE)
abline(v=mean(P_1_10-P_0_10), col="green")
hist(P_1_12-P_0_12,nclass=50,xlim=c(-3,2),
     xlab="6 cov - 12", main="", freq=FALSE)
abline(v=mean(P_1_12-P_0_12), col="blue")

##############################################################
# ---       groups distributions     ----
##############################################################

table_10_cl_b=table(model_data_REGweights_10$S_all_cluster[,2])
table_10_str_b=table(model_data_REGweights_10$S_strata_cluster[,2])
table_10_mix_b=table(model_data_REGweights_10$S_mix_cluster[,2])

table_10_cl_mode=table(model_data_REGweights_10$S_all_cluster[,3])
table_10_str_mode=table(model_data_REGweights_10$S_strata_cluster[,3])
table_10_mix_mode=table(model_data_REGweights_10$S_mix_cluster[,3])

table_12_cl_b=table(model_data_REGweights_12$S_all_cluster[,2])
table_12_str_b=table(model_data_REGweights_12$S_strata_cluster[,2])
table_12_mix_b=table(model_data_REGweights_12$S_mix_cluster[,2])

table_12_cl_mode=table(model_data_REGweights_12$S_all_cluster[,3])
table_12_str_mode=table(model_data_REGweights_12$S_strata_cluster[,3])
table_12_mix_mode=table(model_data_REGweights_12$S_mix_cluster[,3])


P_10_cl_b=lapply(names(table_10_cl_b[table_10_cl_b>2]), function(l)
  P_1_10[model_data_REGweights_10$S_all_cluster[,2]==l]-
    P_0_10[model_data_REGweights_10$S_all_cluster[,2]==l])
P_10_str_b=lapply(names(table_10_str_b[table_10_str_b>2]), function(l)
  P_1_10[model_data_REGweights_10$S_strata_cluster[,2]==l]-
    P_0_10[model_data_REGweights_10$S_strata_cluster[,2]==l])
P_10_mix_b=lapply(names(table_10_mix_b[table_10_mix_b>2]), function(l)
  P_1_10[model_data_REGweights_10$S_mix_cluster[,2]==l]-
    P_0_10[model_data_REGweights_10$S_mix_cluster[,2]==l])

P_10_cl_mode=lapply(names(table_10_cl_mode[table_10_cl_mode>2]), function(l)
  P_1_10[model_data_REGweights_10$S_all_cluster[,3]==l]-
    P_0_10[model_data_REGweights_10$S_all_cluster[,3]==l])
P_10_str_mode=lapply(names(table_10_str_mode[table_10_str_mode>2]), function(l)
  P_1_10[model_data_REGweights_10$S_strata_cluster[,3]==l]-
    P_0_10[model_data_REGweights_10$S_strata_cluster[,3]==l])
P_10_mix_mode=lapply(names(table_10_mix_mode[table_10_mix_mode>2]), function(l)
  P_1_10[model_data_REGweights_10$S_mix_cluster[,3]==l]-
    P_0_10[model_data_REGweights_10$S_mix_cluster[,3]==l])

P_12_cl_b=lapply(names(table_12_cl_b[table_12_cl_b>2]), function(l)
  P_1_12[model_data_REGweights_12$S_all_cluster[,2]==l]-
    P_0_12[model_data_REGweights_12$S_all_cluster[,2]==l])
P_12_str_b=lapply(names(table_12_str_b[table_12_str_b>2]), function(l)
  P_1_12[model_data_REGweights_12$S_strata_cluster[,2]==l]-
    P_0_12[model_data_REGweights_12$S_strata_cluster[,2]==l])
P_12_mix_b=lapply(names(table_12_mix_b[table_12_mix_b>2]), function(l)
  P_1_12[model_data_REGweights_12$S_mix_cluster[,2]==l]-
    P_0_12[model_data_REGweights_12$S_mix_cluster[,2]==l])

P_12_cl_mode=lapply(names(table_12_cl_mode[table_12_cl_mode>2]), function(l)
  P_1_12[model_data_REGweights_12$S_all_cluster[,3]==l]-
    P_0_12[model_data_REGweights_12$S_all_cluster[,3]==l])
P_12_str_mode=lapply(names(table_12_str_mode[table_12_str_mode>2]), function(l)
  P_1_12[model_data_REGweights_12$S_strata_cluster[,3]==l]-
    P_0_12[model_data_REGweights_12$S_strata_cluster[,3]==l])
P_12_mix_mode=lapply(names(table_12_mix_mode[table_12_mix_mode>2]), function(l)
  P_1_12[model_data_REGweights_12$S_mix_cluster[,3]==l]-
    P_0_12[model_data_REGweights_12$S_mix_cluster[,3]==l])

seq_10_cl_b=lapply(names(table_10_cl_b[table_10_cl_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_all_cluster[,2]==l)))
seq_10_str_b=lapply(names(table_10_str_b[table_10_str_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_strata_cluster[,2]==l)))
seq_10_mix_b=lapply(names(table_10_mix_b[table_10_mix_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_mix_cluster[,2]==l)))

seq_10_cl_mode=lapply(names(table_10_cl_mode[table_10_cl_mode>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_all_cluster[,3]==l)))
seq_10_str_mode=lapply(names(table_10_str_mode[table_10_str_mode>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_strata_cluster[,3]==l)))
seq_10_mix_mode=lapply(names(table_10_mix_mode[table_10_mix_mode>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_mix_cluster[,3]==l)))

seq_12_cl_b=lapply(names(table_12_cl_b[table_12_cl_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_12$S_all_cluster[,2]==l)))
seq_12_str_b=lapply(names(table_12_str_b[table_12_str_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_12$S_strata_cluster[,2]==l)))
seq_12_mix_b=lapply(names(table_12_mix_b[table_12_mix_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_12$S_mix_cluster[,2]==l)))

seq_12_cl_mode=lapply(names(table_12_cl_mode[table_12_cl_mode>2]), function(l) 
  rep(l,sum(model_data_REGweights_12$S_all_cluster[,3]==l)))
seq_12_str_mode=lapply(names(table_12_str_mode[table_12_str_mode>2]), function(l) 
  rep(l,sum(model_data_REGweights_12$S_strata_cluster[,3]==l)))
seq_12_mix_mode=lapply(names(table_12_mix_mode[table_12_mix_mode>2]), function(l) 
  rep(l,sum(model_data_REGweights_12$S_mix_cluster[,3]==l)))

# --- boxplot ---

par(mfrow=c(2,6))

boxplot(unlist(P_10_cl_b)~unlist(seq_10_cl_b), 
        main="10 - binder", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")
boxplot(unlist(P_10_str_b)~unlist(seq_10_str_b), 
        main="10 - binder", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")
boxplot(unlist(P_10_mix_b)~unlist(seq_10_mix_b), 
        main="10 - binder", ylab=" ", xlab="P(1)-P(0) | mix")
abline(h=0, col="red")

boxplot(unlist(P_10_cl_mode)~unlist(seq_10_cl_mode), 
        main="10 - mode", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")
boxplot(unlist(P_10_str_mode)~unlist(seq_10_str_mode), 
        main="10 - mode", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")
boxplot(unlist(P_10_mix_mode)~unlist(seq_10_mix_mode), 
        main="10 - mode", ylab=" ", xlab="P(1)-P(0) | mix")
abline(h=0, col="red")

boxplot(unlist(P_12_cl_b)~unlist(seq_12_cl_b), 
        main="12 - binder", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")
boxplot(unlist(P_12_str_b)~unlist(seq_12_str_b), 
        main="12 - binder", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")
boxplot(unlist(P_12_mix_b)~unlist(seq_12_mix_b), 
        main="12 - binder", ylab=" ", xlab="P(1)-P(0) | mix")
abline(h=0, col="red")

boxplot(unlist(P_12_cl_mode)~unlist(seq_12_cl_mode), 
        main="12 - mode", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")
boxplot(unlist(P_12_str_mode)~unlist(seq_12_str_mode), 
        main="12 - mode", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")
boxplot(unlist(P_12_mix_mode)~unlist(seq_12_mix_mode), 
        main="12 - mode", ylab=" ", xlab="P(1)-P(0) | mix")
abline(h=0, col="red")


# --- histogram only for some selected partitions ---

par(mfrow=c(3,1))
for(i in names(table_10_str_b[table_10_str_b>2])){
  hist((P_1_10-P_0_10)[model_data_REGweights_10$S_strata_cluster[,2]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | strata Binder 10", main=paste0("cluster ",i))
  abline(v=mean((P_1_10-P_0_10)[model_data_REGweights_10$S_strata_cluster[,2]==i]), 
         col="red", lwd=2)
}

par(mfrow=c(3,1))
for(i in names(table_10_str_mode[table_10_str_mode>2])){
  hist((P_1_10-P_0_10)[model_data_REGweights_10$S_strata_cluster[,3]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | strata mode 10", main=paste0("cluster ",i))
  abline(v=mean((P_1_10-P_0_10)[model_data_REGweights_10$S_strata_cluster[,3]==i]), 
         col="red", lwd=2)
}

par(mfrow=c(round(sum(table_12_str_b>2)/2),2))
for(i in names(table_12_str_b[table_12_str_b>2])){
  hist((P_1_12-P_0_12)[model_data_REGweights_12$S_strata_cluster[,2]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | strata Binder 12", main=paste0("cluster ",i))
  abline(v=mean((P_1_12-P_0_12)[model_data_REGweights_12$S_strata_cluster[,2]==i]), 
         col="red", lwd=2)
}

par(mfrow=c(3,1))
for(i in names(table_12_str_mode[table_12_str_mode>2])){
  hist((P_1_12-P_0_12)[model_data_REGweights_12$S_strata_cluster[,3]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | strata mode 12", main=paste0("cluster ",i))
  abline(v=mean((P_1_12-P_0_12)[model_data_REGweights_12$S_strata_cluster[,3]==i]), 
         col="red", lwd=2)
}

par(mfrow=c(round(sum(table_12_mix_mode>2)/2),2))
for(i in names(table_12_mix_mode[table_12_mix_mode>2])){
  hist((P_1_12-P_0_12)[model_data_REGweights_12$S_mix_cluster[,3]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | mix mode 12", main=paste0("cluster ",i))
  abline(v=mean((P_1_12-P_0_12)[model_data_REGweights_12$S_mix_cluster[,3]==i]), 
         col="red", lwd=2)
}

##############################################################
# ---       Y | groups distributions     ----
##############################################################

Y_0_10=apply(model_data_REGweights_10$post_Y_0_imp,1,median)
Y_0_12=apply(model_data_REGweights_12$post_Y_0_imp,1,median)

Y_1_10=apply(model_data_REGweights_10$post_Y_1_imp,1,median)
Y_1_12=apply(model_data_REGweights_12$post_Y_1_imp,1,median)

# --- histogram only for some selected partitions ---

par(mfrow=c(3,1))
for(i in names(table_10_str_b[table_10_str_b>2])){
  hist((Y_1_10-Y_0_10)[model_data_REGweights_10$S_strata_cluster[,2]==i], 
       nclass=15, xlim=c(-3,3), freq=FALSE,
       xlab="Y(1)-Y(0) | strata Binder 10", main=paste0("cluster ",i))
  abline(v=mean((Y_1_10-Y_0_10)[model_data_REGweights_10$S_strata_cluster[,2]==i]),
       col="red", lwd=2)
}

par(mfrow=c(3,1))
for(i in names(table_10_str_mode[table_10_str_mode>2])){
  hist((Y_1_10-Y_0_10)[model_data_REGweights_10$S_strata_cluster[,3]==i], 
       nclass=15, xlim=c(-3,3), freq=FALSE,
       xlab="Y(1)-Y(0) | strata mode 10", main=paste0("cluster ",i))
  abline(v=mean((Y_1_10-Y_0_10)[model_data_REGweights_10$S_strata_cluster[,3]==i]),
         col="red", lwd=2)
}

par(mfrow=c(3,1))
for(i in names(table_12_str_mode[table_12_str_mode>2])){
  hist((Y_1_12-Y_0_12)[model_data_REGweights_12$S_strata_cluster[,3]==i], 
       nclass=15, xlim=c(-3,3), freq=FALSE,
       xlab="Y(1)-Y(0) | strata mode 12", main=paste0("cluster ",i))
  abline(v=mean((Y_1_12-Y_0_12)[model_data_REGweights_12$S_strata_cluster[,3]==i]),
         col="red", lwd=2)
}

par(mfrow=c(4,1))
for(i in names(table_12_mix_mode[table_12_mix_mode>2])){
  hist((Y_1_12-Y_0_12)[model_data_REGweights_12$S_mix_cluster[,3]==i], 
       nclass=15, xlim=c(-3,3), freq=FALSE,
       xlab="Y(1)-Y(0) | mix mode 12", main=paste0("cluster ",i))
  abline(v=mean((Y_1_12-Y_0_12)[model_data_REGweights_12$S_mix_cluster[,3]==i]),
         col="red", lwd=2)
}

##############################################################
# ---     covariates analysis    ----
##############################################################

cov_10_str_b=sapply(names(table_10_str_b),function(c)
  apply(dataset_matched[which(model_data_REGweights_10$S_strata_cluster[,2]==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                          "PctFemale","PctPoor")],2,mean))
cov_10_str_mode=sapply(names(table_10_str_mode),function(c)
  apply(dataset_matched[which(model_data_REGweights_10$S_strata_cluster[,3]==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                          "PctFemale","PctPoor")],2,mean))
cov_12_str_mode=sapply(names(table_12_str_mode),function(c)
  apply(dataset_matched[which(model_data_REGweights_12$S_strata_cluster[,3]==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                          "PctFemale","PctPoor")],2,mean))
cov_12_mix_mode=sapply(names(table_12_mix_mode),function(c)
  apply(dataset_matched[which(model_data_REGweights_12$S_mix_cluster[,3]==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                          "PctFemale","PctPoor")],2,mean))

round(cov_10_str_b,3)
round(cov_10_str_mode,3)
round(cov_12_str_mode,3)
round(cov_12_mix_mode,3)


