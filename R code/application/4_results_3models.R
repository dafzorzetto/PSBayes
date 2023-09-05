##############################################################
# ---       analysis regults     ----
##############################################################

#load results
load("C:/Users/dafne/Desktop/DOTTORATO/PROGETTO PhD _2/data/results_3model.RData")

##############################################################

# partitions 

apply(model_data_NO$S_all_cluster,2,table)
apply(model_data_NO$S_strata_cluster,2,table)

apply(model_data_REGweights$S_all_cluster,2,table)
apply(model_data_REGweights$S_strata_cluster,2,table)
apply(model_data_REGweights$S_mix_cluster,2,table)

apply(model_data_REGall$S_all_cluster,2,table)
apply(model_data_REGall$S_strata_cluster,2,table)


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

# no covariates in Y-regression

P_0_NO=apply(model_data_NO$post_P_0,1,median)
P_1_NO=apply(model_data_NO$post_P_1,1,median)

# same covariates of the weights

P_0_w=apply(model_data_REGweights$post_P_0,1,median)
P_1_w=apply(model_data_REGweights$post_P_1,1,median)

# same covariates of the weights

P_0_all=apply(model_data_REGall$post_P_0,1,median)
P_1_all=apply(model_data_REGall$post_P_1,1,median)

##############################################################
# ---        all models        ---
# ---    distributions of P    ---
##############################################################

par(mfrow=c(4,2))
hist(dataset_matched$pm_diff[dataset_matched$a==0],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="observed", main="control", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
hist(dataset_matched$pm_diff[dataset_matched$a==1],nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="observed PM", main="treated", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")

hist(P_0_NO,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="no covariates", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
abline(v=mean(P_0_NO), col="blue")
hist(P_1_NO,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="no covariates", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")
abline(v=mean(P_1_NO), col="blue")

hist(P_0_w,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="6 covariates", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
abline(v=mean(P_0_w), col="green")
hist(P_1_w,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="6 covariates", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")
abline(v=mean(P_1_w), col="green")

hist(P_0_all,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="all covariates", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==0]),col="red")
abline(v=mean(P_0_all), col="pink")
hist(P_1_all,nclass=50,
     xlim=c(min(dataset_matched$pm_diff),max(dataset_matched$pm_diff)),
     xlab="all covariates", main="", freq=FALSE)
abline(v=mean(dataset_matched$pm_diff[dataset_matched$a==1]),col="red")
abline(v=mean(P_1_all), col="pink")

##############################################################

# indirect effect: P(1) - P(0)

par(mfrow=c(3,1))
hist(P_1_NO-P_0_NO,nclass=50,xlim=c(-3,2),
     xlab="no covariates", main="P(1) - P(0)", freq=FALSE)
abline(v=mean(P_1_NO-P_0_NO), col="blue")
hist(P_1_w-P_0_w,nclass=50,xlim=c(-3,2),
     xlab="6 covariates", main="", freq=FALSE)
abline(v=mean(P_1_w-P_0_w), col="green")
hist(P_1_all-P_0_all,nclass=50,,xlim=c(-3,2),
     xlab="all covariates", main="", freq=FALSE)
abline(v=mean(P_1_all-P_0_all), col="pink")


##############################################################
# ---       groups distributions     ----
##############################################################

table_NO_cl=table(model_data_NO$S_all_cluster[,2])
table_w_cl=table(model_data_REGweights$S_all_cluster[,2])
table_all_cl=table(model_data_REGall$S_all_cluster[,2])

table_NO_str=table(model_data_NO$S_strata_cluster[,2])
table_w_str=table(model_data_REGweights$S_strata_cluster[,2])
table_all_str=table(model_data_REGall$S_strata_cluster[,2])

table_w_mix=table(model_data_REGweights$S_mix_cluster[,2])

P_NO_cl=lapply(names(table_NO_cl[table_NO_cl>4]), function(l)
  P_1_NO[model_data_NO$S_all_cluster[,2]==l]-
    P_0_NO[model_data_NO$S_all_cluster[,2]==l])
P_w_cl=lapply(names(table_w_cl[table_w_cl>4]), function(l)
  P_1_w[model_data_REGweights$S_all_cluster[,2]==l]-
    P_0_w[model_data_REGweights$S_all_cluster[,2]==l])
P_all_cl=lapply(names(table_all_cl[table_all_cl>4]), function(l)
  P_1_all[model_data_REGall$S_all_cluster[,2]==l]-
    P_0_all[model_data_REGall$S_all_cluster[,2]==l])

P_NO_str=lapply(names(table_NO_str[table_NO_str>4]), function(l)
  P_1_NO[model_data_NO$S_strata_cluster[,2]==l]-
    P_0_NO[model_data_NO$S_strata_cluster[,2]==l])
P_w_str=lapply(names(table_w_str[table_w_str>4]), function(l)
  P_1_w[model_data_REGweights$S_strata_cluster[,2]==l]-
    P_0_w[model_data_REGweights$S_strata_cluster[,2]==l])
P_all_str=lapply(names(table_all_str[table_all_str>4]), function(l)
  P_1_all[model_data_REGall$S_strata_cluster[,2]==l]-
    P_0_all[model_data_REGall$S_strata_cluster[,2]==l])

P_w_mix=lapply(names(table_w_mix[table_w_mix>4]), function(l)
  P_1_w[model_data_REGweights$S_mix_cluster[,2]==l]-
    P_0_w[model_data_REGweights$S_mix_cluster[,2]==l])

seq_NO_cl=lapply(names(table_NO_cl[table_NO_cl>4]), function(l) 
  rep(l,sum(model_data_NO$S_all_cluster[,2]==l)))
seq_w_cl=lapply(names(table_w_cl[table_w_cl>4]), function(l) 
  rep(l,sum(model_data_REGweights$S_all_cluster[,2]==l)))
seq_all_cl=lapply(names(table_all_cl[table_all_cl>4]), function(l) 
  rep(l,sum(model_data_REGall$S_all_cluster[,2]==l)))

seq_NO_str=lapply(names(table_NO_str[table_NO_str>4]), function(l) 
  rep(l,sum(model_data_NO$S_strata_cluster[,2]==l)))
seq_w_str=lapply(names(table_w_str[table_w_str>4]), function(l) 
  rep(l,sum(model_data_REGweights$S_strata_cluster[,2]==l)))
seq_all_str=lapply(names(table_all_str[table_all_str>4]), function(l) 
  rep(l,sum(model_data_REGall$S_strata_cluster[,2]==l)))

seq_w_mix=lapply(names(table_w_mix[table_w_mix>4]), function(l) 
  rep(l,sum(model_data_REGweights$S_mix_cluster[,2]==l)))

# --- boxplot ---

par(mfrow=c(1,3))
boxplot(unlist(P_NO_cl)~unlist(seq_NO_cl), 
        main="no covariates", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")
boxplot(unlist(P_w_cl)~unlist(seq_w_cl), 
        main="6 covariates", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")
boxplot(unlist(P_all_cl)~unlist(seq_all_cl), 
        main="all caovariates", ylab=" ", xlab="P(1)-P(0) | cluster")
abline(h=0, col="red")

par(mfrow=c(1,3))
boxplot(unlist(P_NO_str)~unlist(seq_NO_str), 
        main="no covariates", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")
boxplot(unlist(P_w_str)~unlist(seq_w_str), 
        main="6 covariates", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")
boxplot(unlist(P_all_str)~unlist(seq_all_str), 
        main="all caovariates", ylab=" ", xlab="P(1)-P(0) | strata")
abline(h=0, col="red")

boxplot(unlist(P_w_mix)~unlist(seq_w_mix), 
        main="6 covariates", ylab=" ", xlab="P(1)-P(0) | mix")
abline(h=0, col="red")


# --- mode ---- 

table_mode_cl=table(model_data_REGweights$S_all_cluster[,3])
table_mode_str=table(model_data_REGweights$S_strata_cluster[,3])
table_mode_mix=table(model_data_REGweights$S_mix_cluster[,3])

P_mode_cl=lapply(names(table_mode_cl[table_mode_cl>4]), function(l)
  P_1_w[model_data_REGweights$S_all_cluster[,3]==l]-
    P_0_w[model_data_REGweights$S_all_cluster[,3]==l])
P_mode_str=lapply(names(table_mode_str[table_mode_str>4]), function(l)
  P_1_w[model_data_REGweights$S_strata_cluster[,3]==l]-
    P_0_w[model_data_REGweights$S_strata_cluster[,3]==l])
P_mode_mix=lapply(names(table_mode_mix[table_mode_mix>4]), function(l)
  P_1_w[model_data_REGweights$S_mix_cluster[,3]==l]-
    P_0_w[model_data_REGweights$S_mix_cluster[,3]==l])

seq_mode_cl=lapply(names(table_mode_cl[table_mode_cl>4]), function(l) 
  rep(l,sum(model_data_REGweights$S_all_cluster[,3]==l)))
seq_mode_str=lapply(names(table_mode_str[table_mode_str>4]), function(l) 
  rep(l,sum(model_data_REGweights$S_strata_cluster[,3]==l)))
seq_mode_mix=lapply(names(table_mode_mix[table_mode_mix>4]), function(l) 
  rep(l,sum(model_data_REGweights$S_mix_cluster[,3]==l)))

par(mfrow=c(1,3))
boxplot(unlist(P_mode_cl)~unlist(seq_mode_cl), 
        main="6 cov. MODE", ylab=" ", xlab="P(1)-P(0) | cluster ")
abline(h=0, col="red")
boxplot(unlist(P_mode_str)~unlist(seq_mode_str), 
        main="6 cov. MODE", ylab=" ", xlab="P(1)-P(0) | strata ")
abline(h=0, col="red")
boxplot(unlist(P_mode_mix)~unlist(seq_mode_mix), 
        main="6 cov. MODE", ylab=" ", xlab="P(1)-P(0) | mix ")
abline(h=0, col="red")

# histogram
par(mfrow=c(round(sum(table_w_str>4)/2),2))
for(i in names(table_w_str[table_w_str>4])){
  hist((P_1_w-P_0_w)[model_data_REGweights$S_strata_cluster[,2]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | cluster", main=paste0("cluster ",i))
  abline(v=mean((P_1_w-P_0_w)[model_data_REGweights$S_strata_cluster[,2]==i]), 
         col="red", lwd=2)
}

par(mfrow=c(round(sum(table_mode_str>4)/2),2))
for(i in names(table_mode_str[table_mode_str>4])){
  hist((P_1_w-P_0_w)[model_data_REGweights$S_strata_cluster[,3]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | cluster", main=paste0("cluster ",i))
  abline(v=mean((P_1_w-P_0_w)[model_data_REGweights$S_strata_cluster[,3]==i]), 
         col="red", lwd=2)
}
par(mfrow=c(round(sum(table_mode_mix>4)/2),2))
for(i in names(table_mode_mix[table_mode_mix>4])){
  hist((P_1_w-P_0_w)[model_data_REGweights$S_mix_cluster[,3]==i], 
       nclass=20, xlim=c(-3,2.5), freq=FALSE,
       xlab="P(1)-P(0) | mix", main=paste0("cluster ",i))
  abline(v=mean((P_1_w-P_0_w)[model_data_REGweights$S_mix_cluster[,3]==i]), 
         col="red", lwd=2)
}

##############################################################
# ---       Y | groups distributions     ----
##############################################################

Y_0_NO=apply(model_data_NO$post_Y_0_imp,1,median)
Y_0_w=apply(model_data_REGweights$post_Y_0_imp,1,median)
Y_0_all=apply(model_data_REGall$post_Y_0_imp,1,median)

Y_1_NO=apply(model_data_NO$post_Y_1_imp,1,median)
Y_1_w=apply(model_data_REGweights$post_Y_1_imp,1,median)
Y_1_all=apply(model_data_REGall$post_Y_1_imp,1,median)

# histograms

par(mfrow=c(round(sum(table_NO_str>4)/2),2))
for(i in names(table_NO_str[table_NO_str>4])){
  hist((Y_1_NO-Y_0_NO)[model_data_NO$S_strata_cluster[,2]==i], 
       nclass=10, xlim=c(-1,1), freq=FALSE,
       xlab="Y(1)-Y(0) | cluster", main=paste0("cluster ",i))
}
par(mfrow=c(round(sum(table_w_str>4)/2),2))
for(i in names(table_w_str[table_w_str>4])){
  hist((Y_1_w-Y_0_w)[model_data_REGweights$S_strata_cluster[,2]==i], 
       nclass=20, xlim=c(-2.5,2), freq=FALSE,
       xlab="Y(1)-Y(0) | cluster", main=paste0("cluster ",i))
  abline(v=mean((Y_1_w-Y_0_w)[model_data_REGweights$S_strata_cluster[,2]==i]), 
         col="red", lwd=2)
}
par(mfrow=c(round(sum(table_all_str>4)/2),2))
for(i in names(table_all_str[table_all_str>4])){
  hist((Y_1_all-Y_0_all)[model_data_REGall$S_strata_cluster[,2]==i], 
       nclass=10, xlim=c(-5,4), freq=FALSE,
       xlab="Y(1)-Y(0) | cluster", main=paste0("cluster ",i))
}

# mode - strata
par(mfrow=c(round(sum(table_mode_str>4)/2),2))
for(i in names(table_mode_str[table_mode_str>4])){
  hist((Y_1_w-Y_0_w)[model_data_REGweights$S_strata_cluster[,3]==i], 
       nclass=20, xlim=c(-2.5,2), freq=FALSE,
       xlab="Y(1)-Y(0) | cluster", main=paste0("cluster ",i))
  abline(v=mean((Y_1_w-Y_0_w)[model_data_REGweights$S_strata_cluster[,3]==i]), 
         col="red", lwd=2)
}

# clusters
par(mfrow=c(3,2))
for(i in names(table_w_cl[table_w_cl>4])){
  hist((Y_1_w-Y_0_w)[model_data_REGweights$S_all_cluster[,2]==i], 
       nclass=20, xlim=c(-2.5,2), freq=FALSE,
       xlab="Y(1)-Y(0) | cluster", main=paste0("cluster ",i))
}


##############################################################

cov_w_str=sapply(names(table_w_str),function(c)
  apply(dataset_matched[which(model_data_REGweights$S_strata_cluster[,2]==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban","PctFemale","PctPoor")],
        2,mean))

round(cov_w_str,3)
