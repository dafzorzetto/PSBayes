##############################################################
# ---       estimation model     ----
##############################################################


#load matched data
load("dataset_matched.RData")
#load("C:/Users/dafne/Desktop/DOTTORATO/PROGETTO PhD _2/data/dataset_matched.RData")


#library
library(mvtnorm)
library(CholWishart)
library(parallel)
library(truncnorm)
library(invgamma)
library(BNPmix)

###################################################################

# iterations
R=6000
R_burnin=3000

# max number clusters
n_cluster=10

###################################################################

library(corrgram)
#pdf(file="rorr_cov.pdf")
#corrgram(dataset_matched[,4:18], order=TRUE,
#         lower.panel=panel.shade, upper.panel=panel.conf)
#dev.off()

# prepar dataset --- common part

matrix_X=dataset_matched[,c("PctBlack","PctHisp","PctHighSchool",
                            "PctUrban","PctFemale","PctPoor")]
matrix_X=as.matrix(matrix_X)

T_var=dataset_matched$a
P_obs=dataset_matched$pm_diff
Y_obs=dataset_matched$all_causes*1000

###################################################################
###################################################################

# Provare:

# - tutte le covariate
# - diversi death causes


###################################################################
###################################################################

BNP_model_REG<-function(c){
  set.seed(c)
  n=length(T_var)
  
  #useful quantities
  n_X=dim(matrix_X)[2]
  T1=which(T_var==1)
  n1=length(T1)
  T0=which(T_var==0)
  n0=length(T0)
  P_mis=rep(0,n)
  P_obs_0=P_obs[T0]
  P_obs_1=P_obs[T1]
  Y_obs_0=Y_obs[T0]
  Y_obs_1=Y_obs[T1]
  
  # imputation
  P_mis=rep(0,n)
  P_imp=rep(0,n)
  Y_0_imp=rep(0,n)
  Y_1_imp=rep(0,n)
  
  #prior
  p_beta=c(0,5)
  p_sigma=c(6,0.1)
  #p_eta=c(12,5)
  p_eta=c(-3,5)
  p_theta=c(0,5)
  p_lambda=c(0,2,0,1)
  
  #inizializzazione
  beta_0=rep(0,n_X*(n_cluster-1))
  beta_1=rep(0,n_X*(n_cluster-1))
  sigma=rep(1,n_cluster)
  eta=seq(-3,-3+0.5*(n_cluster-1),0.5)
  theta_0=rep(0,2+(dim(matrix_REGR)[2]))
  theta_1=rep(0,4+(dim(matrix_REGR)[2]))
  lambda=rep(1,2)
  
  # some useful functions
  omega=function(X,beta){
    L=n_cluster-1
    alpha=c(sapply(1:L, function(l) pnorm(beta[(l*n_X-(n_X-1)):(l*n_X)]%*%X)),1)
    return(c(alpha[1],sapply(2:(L+1), function(l) alpha[l]*prod(1-alpha[1:(l-1)]))))
  }
  
  moda <- function(v) {
    as.numeric(names(which.max(table(v))))
  }
  
  partition_all_cluster<-function(xi_0,xi_1){
    grid_couple=expand.grid(1:n_cluster,1:n_cluster)
    
    #cartesian product of cluster allocation
    part=sapply(1:n, function(i)
      sapply(1:(R-R_burnin), function(r) 
        which(grid_couple$Var1==xi_0[i,r] & grid_couple$Var2==xi_1[i,r])))
    
    #point estimation
    WG=partition.BNPdens(list(clust=part),dist="VI")$partitions[1,]
    Binder=partition.BNPdens(list(clust=part),dist="Binder")$partitions[1,]
    tab=cbind(WG,Binder, apply(part,2, moda))
    dimnames(tab)=NULL
    return(tab)
  }
  
  partition_strata_cluster<-function(xi_0,xi_1, eta){
    #stratat allocation
    part=sapply(1:n, function(i)
      sapply(1:(R-R_burnin), function(r) 
        ifelse(xi_0[i,r]==xi_1[i,r], 0, 
               1*(eta[xi_1[i,r],r]>eta[xi_0[i,r],r])+(-1)*(eta[xi_0[i,r],r]>eta[xi_1[i,r],r]))))
    
    #point estimation
    WG=partition.BNPdens(list(clust=part),dist="VI")$partitions[1,]
    Binder=partition.BNPdens(list(clust=part),dist="Binder")$partitions[1,]
    tab=cbind(WG,Binder, apply(part,2, moda))
    dimnames(tab)=NULL
    return(tab)
  }
  
  partition_mix_cluster<-function(xi_0,xi_1,eta){
    grid_couple=expand.grid(1:n_cluster,1:n_cluster)
    
    #cartesian product of cluster allocation
    part=sapply(1:n, function(i)
      sapply(1:(R-R_burnin), function(r) 
        ifelse(xi_0[i,r]==xi_1[i,r], 0,
               which(grid_couple$Var1==xi_0[i,r] & grid_couple$Var2==xi_1[i,r]))))
    
    #point estimation
    WG=partition.BNPdens(list(clust=part),dist="VI")$partitions[1,]
    Binder=partition.BNPdens(list(clust=part),dist="Binder")$partitions[1,]
    tab=cbind(WG,Binder, apply(part,2, moda))
    dimnames(tab)=NULL
    return(tab)
  }
  
  
  xi_0i=sample(1:n_cluster,n0,replace=TRUE)
  xi_1i=sample(1:n_cluster,n1,replace=TRUE)
  t_xi_0=table(xi_0i)
  t_xi_1=table(xi_1i)
  eta_sogg_0=eta[xi_0i]
  eta_sogg_1=eta[xi_1i]
  Z_0=rep(list(rep(NA,n_cluster)),n0)
  Z_1=rep(list(rep(NA,n_cluster)),n1)
  
  # acceptation rate for the metropolis proposal step (variance Y-model)
  acc_0=0
  acc_1=0
  
  #chains:
  post_eta=matrix(NA,ncol=R-R_burnin,nrow=n_cluster)
  post_var=matrix(NA,ncol=R-R_burnin,nrow=n_cluster)
  post_xi_0=matrix(NA,ncol=R-R_burnin,nrow=n)     
  post_xi_1=matrix(NA,ncol=R-R_burnin,nrow=n)     
  post_beta=matrix(NA,ncol=R-R_burnin, nrow=length(beta_0)+length(beta_1))
  post_theta=matrix(NA,ncol=R-R_burnin,nrow=6+2*(dim(matrix_REGR)[2]))
  post_lambda=matrix(NA,ncol=R-R_burnin, nrow=2)
  
  post_P_0=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_P_1=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_Y_0_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_Y_1_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  
  for (r in 1:R){
    # ----------- Cluster Specific Parameters  ------------
    
    #ETA:
    for(l in 1:n_cluster){                                      
      v_eta_inv=(t_xi_0[l]+t_xi_1[l])/sigma[l]+1/p_eta[2]
      m_eta=(sum(P_obs_0[xi_0i==l])+sum(P_obs_1[xi_1i==l]))/sigma[l]+p_eta[1]/p_eta[2]
      eta[l]=rnorm(1,mean=1/v_eta_inv*(m_eta),
                   sd=sqrt(1/v_eta_inv))
    }
    
    eta_sogg_0=eta[xi_0i]
    eta_sogg_1=eta[xi_1i]
    
    #SIGMA:
    #sigma
    sigma=sapply(1:n_cluster, function(l) 
      rinvgamma(1,p_sigma[1]+(sum(xi_0i==l)+sum(xi_1i==l))/2,
                p_sigma[2]+(sum((P_obs_0[xi_0i==l]-eta[l])^2)+sum((P_obs_1[xi_1i==l]-eta[l])^2))/2))
    
    # ----------- Cluster Allocation  ------------
    
    #OMEGA:
    #omega_0
    omega_0=sapply(1:n0, function(i) omega(X=matrix_X[T0[i],],beta=beta_0))
    #omega_1
    omega_1=sapply(1:n1, function(i) omega(X=matrix_X[T1[i],],beta=beta_1))
    
    #Xi:
    #Xi_0
    dmn_0=sapply(1:n0, function(i) sapply(1:n_cluster, function(l) dnorm(P_obs_0[i], eta[l], sqrt(sigma[l]), log=TRUE)))+log(omega_0)
    dmn_0[which(is.nan(dmn_0))]=-100
    xi=sapply(1:n0, function(i) rmultinom(1,1,exp(dmn_0[,i])))
    xi_0i=sapply(1:n0, function(i) xi[,i]%*%(1:n_cluster))
    t_xi_0=apply(xi, 1, sum)
    
    #Xi_1
    dmn_1=sapply(1:n1, function(i) sapply(1:n_cluster, function(l) dnorm(P_obs_1[i], eta[l], sqrt(sigma[l]), log=TRUE)))+log(omega_1)
    dmn_1[which(is.nan(dmn_1))]=-100
    xi=sapply(1:n1, function(i) rmultinom(1,1,exp(dmn_1[,i])))
    xi_1i=sapply(1:n1, function(i) xi[,i]%*%(1:n_cluster))
    t_xi_1=apply(xi, 1, sum)
    
    # ----------- Augmentation Scheme  ------------
    
    # Z
    pesi_0=t(omega_0)
    mu_z=cbind((pesi_0[,1]),(pesi_0[,2]/(1-pesi_0[,1])))
    if (n_cluster>3){
      mu_z=cbind(mu_z,sapply(3:(n_cluster-1), function(l) (pesi_0[,l]/(1-apply(pesi_0[,1:(l-1)],1,sum)))))
    }
    mu_z[which(is.nan(mu_z))]=1
    mu_z[which(mu_z>1)]=1
    mu_z=mu_z-9.9e-15*(mu_z>(1-1e-16))
    #Z_0:
    for (i in 1:n0){
      for (l in 1:(min(xi_0i[i],n_cluster-1))) {
        if(l>1){
          if (l<xi_0i[i]){
            Z_0[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_0[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }else{
          if (l<xi_0i[i]){
            Z_0[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_0[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }
      }
    }
    
    #Z_1:
    pesi_1=t(omega_1)
    mu_z=cbind((pesi_1[,1]),(pesi_1[,2]/(1-pesi_1[,1])))
    if (n_cluster>3){
      mu_z=cbind(mu_z,sapply(3:(n_cluster-1), function(l) (pesi_1[,l]/(1-apply(pesi_1[,1:(l-1)],1,sum)))))
    }
    mu_z[which(is.nan(mu_z))]=1
    mu_z[which(mu_z>1)]=1
    mu_z=mu_z-9.9e-15*(mu_z>(1-1e-16))
    for (i in 1:n1){
      for (l in 1:(min(xi_1i[i],n_cluster-1))) {
        if(l>1){
          if (l<xi_1i[i]){
            Z_1[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_1[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }else{
          if (l<xi_1i[i]){
            Z_1[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_1[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }
      }
    }
    
    # ----------- Confounder-Dependent Weights  ------------
    
    #BETA:
    #beta_0
    gruppi=which(t_xi_0!=0)
    if (max(gruppi)==n_cluster)
      gruppi=gruppi[-length(gruppi)]
    for (l in gruppi){
      val=which(xi_0i>=l)
      z_tilde=unlist(sapply(val, function(i) Z_0[[i]][l]))
      x_tilde=matrix(matrix_X[T0[val],],ncol=n_X)
      V=solve(diag(n_X)/p_beta[2]+t(x_tilde)%*%x_tilde)
      beta_0[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,V%*%(1/p_beta[2]*(diag(n_X)%*%rep(p_beta[1],n_X))+t(x_tilde)%*%z_tilde),V)[,1:n_X]
    }
    vuoti=which(t_xi_0==0)
    if (length(vuoti)>0){
      if (max(vuoti)==n_cluster)
        vuoti=vuoti[-length(vuoti)]
    }
    for (l in vuoti){
      beta_0[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,rep(p_beta[1],n_X),diag(n_X)/p_beta[2])[,1:n_X]
    }
    #beta_1
    gruppi=which(t_xi_1!=0)
    if (max(gruppi)==n_cluster)
      gruppi=gruppi[-length(gruppi)]
    for (l in gruppi){
      val=which(xi_1i>=l)
      z_tilde=unlist(sapply(val, function(i) Z_1[[i]][l]))
      x_tilde=matrix(matrix_X[T1[val],],ncol=n_X)
      V=solve(diag(n_X)/p_beta[2]+t(x_tilde)%*%x_tilde)
      beta_1[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,V%*%(1/p_beta[2]*(diag(n_X)%*%rep(p_beta[1],n_X))+t(x_tilde)%*%z_tilde),V)[,1:n_X]
    }
    vuoti=which(t_xi_1==0)
    if (length(vuoti)>0){
      if (max(vuoti)==n_cluster)
        vuoti=vuoti[-length(vuoti)]
    }
    for (l in vuoti){
      beta_1[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,rep(p_beta[1],n_X),diag(n_X)/p_beta[2])[,1:n_X]
    }
    
    # --- imputing post-treatment missing variables ----
    
    #P_MIS:
    #t0->mis_P_1
    om_0=sapply(T0, function(i) omega(X=matrix_X[i,],beta=beta_1))
    om_0[which(is.nan(om_0))]=0
    p0=sapply(1:n0, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_0[,i])))
    P_mis[T0]=rnorm(T0,eta[p0],sigma[p0])
    #t1->mis_P_0
    om_1=sapply(T1, function(i) omega(X=matrix_X[i,],beta=beta_0))
    om_1[which(is.nan(om_1))]=0
    p1=sapply(1:n1, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_1[,i])))
    v_2=exp(lambda[1]+lambda[2]*P_obs_1)/((theta_1[3]+theta_1[4]*P_obs_1)^2)
    var_inv=1/sigma[p1]+1/v_2
    m_2=(Y_obs_1-theta_1[1]-theta_1[2]*P_obs_1)/(theta_1[3]+theta_1[4]*P_obs_1)
    P_mis[T1]=rnorm(T1,1/(var_inv)*(eta[p1]/sigma[p1]+m_2/v_2),1/(var_inv))
    
    # --- imputing post-treatment observed variables ----
    
    #P_IMP:
    #t=0
    #om_0=sapply(T0, function(i) omega(X=matrix_X[i,],beta=beta_0))
    #om_0[which(is.nan(om_0))]=0
    #p0__=sapply(1:n0, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_0[,i])))
    #P_imp[T0]=rnorm(T0,eta[p0__],sigma[p0__])
    #t=1
    #om_1=sapply(T1, function(i) omega(X=matrix_X[i,],beta=beta_1))
    #om_1[which(is.nan(om_1))]=0
    #p1__=sapply(1:n1, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_1[,i])))
    #P_imp[T1]=rnorm(T1,eta[p1__],sigma[p1__])
    
    # --- Y model ----
    
    #means regression
    P_tilde=cbind(1,P_obs_0,matrix_REGR[T0,])
    M=t(P_tilde)%*%diag(1/rep((exp(lambda[1]))^2,n0))%*%t(t(Y_obs_0))+
      p_theta[1]/(p_theta[2]^2)
    V=t(P_tilde)%*%diag(1/rep((exp(lambda[1]))^2,n0))%*%P_tilde+
      diag(rep(1/(p_theta[2]^2),2+(dim(matrix_REGR)[2])))
    V[lower.tri(V)] = t(V)[lower.tri(V)]
    theta_0=c(rmvnorm(1,solve(V)%*%M,solve(V)))
    
    P_tilde=cbind(1,P_obs_1,P_mis[T1],P_mis[T1]*P_obs_1,matrix_REGR[T1,])
    M=t(P_tilde)%*%diag(1/((exp(lambda[1]+lambda[2]*P_obs_1))^2))%*%t(t(Y_obs_1))+
      p_theta[1]/(p_theta[2]^2)
    V=t(P_tilde)%*%diag(1/((exp(lambda[1]+lambda[2]*P_obs_1))^2))%*%P_tilde+
      diag(rep(1/(p_theta[2]^2),4+(dim(matrix_REGR)[2])))
    V[lower.tri(V)] = t(V)[lower.tri(V)]
    theta_1=c(rmvnorm(1,solve(V)%*%M,solve(V)))
    
    # variances
    l0_star=rnorm(1,p_lambda[1],p_lambda[2])  
    l1_star=rnorm(1,p_lambda[3],p_lambda[4])  
    
    mean_M=c(cbind(1,P_obs_0,matrix_REGR[T0,])%*%theta_0,
             P_tilde%*%theta_1)
    var_0=c(rep(exp(lambda[1]),n0),exp(lambda[1]+lambda[2]*P_obs_1))
    var_star_0=c(rep(exp(l0_star),n0),exp(l0_star+lambda[2]*P_obs_1))
    
    a0=min(1, exp(sum(dnorm(c(Y_obs_0,Y_obs_1),mean_M,var_star_0,log=T))-
                    sum(dnorm(c(Y_obs_0,Y_obs_1),mean_M,var_0,log=T))))
    if (runif(1)<a0){
      lambda[1]=l0_star
      acc_0=acc_0+1
    } 
    
    var_star_1=exp(lambda[1]+l1_star*P_obs_1)
    
    a1=min(1, exp(sum(dnorm(Y_obs_1,mean_M[(n0+1):n],var_star_1,log=T))-
                    sum(dnorm(Y_obs_1,mean_M[(n0+1):n],var_0[(n0+1):n],log=T))))
    if (runif(1)<a1){
      lambda[2]=l1_star
      acc_1=acc_1+1
    } 
    
    #posteriors
    if(r>R_burnin){
      # --- imputing outcome ----
      P_0_r=rep(NA,n)
      P_0_r[T0]=P_obs_0
      P_0_r[T1]=P_mis[T1]
      Y_0_imp=rnorm(n,cbind(1,P_0_r,matrix_REGR)%*%t(t(theta_0)),
                    exp(lambda[1]))
      P_1_r=rep(NA,n)
      P_1_r[T1]=P_obs_1
      P_1_r[T0]=P_mis[T0]
      Y_1_imp=rnorm(n,cbind(1,P_1_r,P_0_r,P_0_r*P_1_r,matrix_REGR)%*%t(t(theta_1)),
                    exp(lambda[1]+lambda[2]*P_1_r))
      
      # P saved
      P_0_r[T0]=P_imp[T0]
      P_1_r[T1]=P_imp[T1]
      
      # --- save info ----
      post_eta[,r-R_burnin]=eta
      post_var[,r-R_burnin]=sigma
      post_beta[,r-R_burnin]=c(beta_0,beta_1)
      post_theta[,r-R_burnin]=c(theta_0,theta_1)
      post_lambda[,r-R_burnin]=lambda
      post_xi_0[T0,r-R_burnin]=xi_0i
      post_xi_0[T1,r-R_burnin]=p1
      post_xi_1[T1,r-R_burnin]=xi_1i
      post_xi_1[T0,r-R_burnin]=p0
      post_P_0[,r-R_burnin]=P_0_r
      post_P_1[,r-R_burnin]=P_1_r
      post_Y_0_imp[,r-R_burnin]=Y_0_imp
      post_Y_1_imp[,r-R_burnin]=Y_1_imp
    }
    
    #print(r)
    if (r%%500==0) print(r)
  }
  
  #partitions
  S_all_cluster=partition_all_cluster(xi_0=post_xi_0,xi_1=post_xi_1)
  S_strata_cluster=partition_strata_cluster(xi_0=post_xi_0,xi_1=post_xi_1, eta=post_eta)
  S_mix_cluster=partition_mix_cluster(xi_0=post_xi_0,xi_1=post_xi_1, eta=post_eta)
  
  print(paste0("campione n ",c))
  
  return(list(post_eta=post_eta, post_var=post_var,
              post_beta=post_beta,
              S_all_cluster=S_all_cluster, S_strata_cluster=S_strata_cluster, S_mix_cluster=S_mix_cluster,
              post_theta=post_theta, post_lambda=post_lambda,
              acc_rate_0=acc_0/R, acc_rate_1=acc_1/R,
              post_xi_0=post_xi_0, post_xi_1=post_xi_1,
              post_P_0=post_P_0, post_P_1=post_P_1,
              post_Y_0_imp=post_Y_0_imp, post_Y_1_imp=post_Y_1_imp
  ))
}  

BNP_model_NO<-function(c){
  set.seed(c)
  n=length(T_var)
  
  #useful quantities
  n_X=dim(matrix_X)[2]
  T1=which(T_var==1)
  n1=length(T1)
  T0=which(T_var==0)
  n0=length(T0)
  P_mis=rep(0,n)
  P_obs_0=P_obs[T0]
  P_obs_1=P_obs[T1]
  Y_obs_0=Y_obs[T0]
  Y_obs_1=Y_obs[T1]
  
  # imputation
  P_mis=rep(0,n)
  P_imp=rep(0,n)
  Y_0_imp=rep(0,n)
  Y_1_imp=rep(0,n)
  
  #prior
  p_beta=c(0,5)
  p_sigma=c(6,0.1)
  #p_eta=c(12,5)
  p_eta=c(-3,5)
  p_theta=c(0,5)
  p_lambda=c(0,2,0,1)
  
  #inizializzazione
  beta_0=rep(0,n_X*(n_cluster-1))
  beta_1=rep(0,n_X*(n_cluster-1))
  sigma=rep(1,n_cluster)
  eta=seq(-3,-3+0.5*(n_cluster-1),0.5)
  theta_0=rep(0,2)
  theta_1=rep(0,4)
  lambda=rep(1,2)
  
  # some useful functions
  omega=function(X,beta){
    L=n_cluster-1
    alpha=c(sapply(1:L, function(l) pnorm(beta[(l*n_X-(n_X-1)):(l*n_X)]%*%X)),1)
    return(c(alpha[1],sapply(2:(L+1), function(l) alpha[l]*prod(1-alpha[1:(l-1)]))))
  }
  
  partition_all_cluster<-function(xi_0,xi_1){
    grid_couple=expand.grid(1:n_cluster,1:n_cluster)
    
    #cartesian product of cluster allocation
    part=sapply(1:n, function(i)
      sapply(1:(R-R_burnin), function(r) 
        which(grid_couple$Var1==xi_0[i,r] & grid_couple$Var2==xi_1[i,r])))
    
    #point estimation
    WG=partition.BNPdens(list(clust=part),dist="VI")$partitions[1,]
    Binder=partition.BNPdens(list(clust=part),dist="Binder")$partitions[1,]
    tab=cbind(WG,Binder)
    dimnames(tab)=NULL
    return(tab)
  }
  
  partition_strata_cluster<-function(xi_0,xi_1, eta){
    #stratat allocation
    part=sapply(1:n, function(i)
      sapply(1:(R-R_burnin), function(r) 
        ifelse(xi_0[i,r]==xi_1[i,r], 0, 
               1*(eta[xi_1[i,r],r]>eta[xi_0[i,r],r])+(-1)*(eta[xi_0[i,r],r]>eta[xi_1[i,r],r]))))
    
    #point estimation
    WG=partition.BNPdens(list(clust=part),dist="VI")$partitions[1,]
    Binder=partition.BNPdens(list(clust=part),dist="Binder")$partitions[1,]
    tab=cbind(WG,Binder)
    dimnames(tab)=NULL
    return(tab)
    
  }
  
  
  xi_0i=sample(1:n_cluster,n0,replace=TRUE)
  xi_1i=sample(1:n_cluster,n1,replace=TRUE)
  t_xi_0=table(xi_0i)
  t_xi_1=table(xi_1i)
  eta_sogg_0=eta[xi_0i]
  eta_sogg_1=eta[xi_1i]
  Z_0=rep(list(rep(NA,n_cluster)),n0)
  Z_1=rep(list(rep(NA,n_cluster)),n1)
  
  # acceptation rate for the metropolis proposal step (variance Y-model)
  acc_0=0
  acc_1=0
  
  #chains:
  post_eta=matrix(NA,ncol=R-R_burnin,nrow=n_cluster)
  post_var=matrix(NA,ncol=R-R_burnin,nrow=n_cluster)
  post_xi_0=matrix(NA,ncol=R-R_burnin,nrow=n)     
  post_xi_1=matrix(NA,ncol=R-R_burnin,nrow=n)     
  post_beta=matrix(NA,ncol=R-R_burnin, nrow=length(beta_0)+length(beta_1))
  post_theta=matrix(NA,ncol=R-R_burnin,nrow=6)
  post_lambda=matrix(NA,ncol=R-R_burnin, nrow=2)
  
  post_P_0=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_P_1=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_Y_0_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_Y_1_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  
  for (r in 1:R){
    # ----------- Cluster Specific Parameters  ------------
    
    #ETA:
    for(l in 1:n_cluster){                                      
      v_eta_inv=(t_xi_0[l]+t_xi_1[l])/sigma[l]+1/p_eta[2]
      m_eta=(sum(P_obs_0[xi_0i==l])+sum(P_obs_1[xi_1i==l]))/sigma[l]+p_eta[1]/p_eta[2]
      eta[l]=rnorm(1,mean=1/v_eta_inv*(m_eta),
                   sd=sqrt(1/v_eta_inv))
    }
    
    eta_sogg_0=eta[xi_0i]
    eta_sogg_1=eta[xi_1i]
    
    #SIGMA:
    #sigma
    sigma=sapply(1:n_cluster, function(l) 
      rinvgamma(1,p_sigma[1]+(sum(xi_0i==l)+sum(xi_1i==l))/2,
                p_sigma[2]+(sum((P_obs_0[xi_0i==l]-eta[l])^2)+sum((P_obs_1[xi_1i==l]-eta[l])^2))/2))
    
    # ----------- Cluster Allocation  ------------
    
    #OMEGA:
    #omega_0
    omega_0=sapply(1:n0, function(i) omega(X=matrix_X[T0[i],],beta=beta_0))
    #omega_1
    omega_1=sapply(1:n1, function(i) omega(X=matrix_X[T1[i],],beta=beta_1))
    
    #Xi:
    #Xi_0
    dmn_0=sapply(1:n0, function(i) sapply(1:n_cluster, function(l) dnorm(P_obs_0[i], eta[l], sqrt(sigma[l]), log=TRUE)))+log(omega_0)
    dmn_0[which(is.nan(dmn_0))]=-100
    xi=sapply(1:n0, function(i) rmultinom(1,1,exp(dmn_0[,i])))
    xi_0i=sapply(1:n0, function(i) xi[,i]%*%(1:n_cluster))
    t_xi_0=apply(xi, 1, sum)
    
    #Xi_1
    dmn_1=sapply(1:n1, function(i) sapply(1:n_cluster, function(l) dnorm(P_obs_1[i], eta[l], sqrt(sigma[l]), log=TRUE)))+log(omega_1)
    dmn_1[which(is.nan(dmn_1))]=-100
    xi=sapply(1:n1, function(i) rmultinom(1,1,exp(dmn_1[,i])))
    xi_1i=sapply(1:n1, function(i) xi[,i]%*%(1:n_cluster))
    t_xi_1=apply(xi, 1, sum)
    
    # ----------- Augmentation Scheme  ------------
    
    # Z
    pesi_0=t(omega_0)
    mu_z=cbind((pesi_0[,1]),(pesi_0[,2]/(1-pesi_0[,1])))
    if (n_cluster>3){
      mu_z=cbind(mu_z,sapply(3:(n_cluster-1), function(l) (pesi_0[,l]/(1-apply(pesi_0[,1:(l-1)],1,sum)))))
    }
    mu_z[which(is.nan(mu_z))]=1
    mu_z[which(mu_z>1)]=1
    mu_z=mu_z-9.9e-15*(mu_z>(1-1e-16))
    #Z_0:
    for (i in 1:n0){
      for (l in 1:(min(xi_0i[i],n_cluster-1))) {
        if(l>1){
          if (l<xi_0i[i]){
            Z_0[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_0[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }else{
          if (l<xi_0i[i]){
            Z_0[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_0[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }
      }
    }
    
    #Z_1:
    pesi_1=t(omega_1)
    mu_z=cbind((pesi_1[,1]),(pesi_1[,2]/(1-pesi_1[,1])))
    if (n_cluster>3){
      mu_z=cbind(mu_z,sapply(3:(n_cluster-1), function(l) (pesi_1[,l]/(1-apply(pesi_1[,1:(l-1)],1,sum)))))
    }
    mu_z[which(is.nan(mu_z))]=1
    mu_z[which(mu_z>1)]=1
    mu_z=mu_z-9.9e-15*(mu_z>(1-1e-16))
    for (i in 1:n1){
      for (l in 1:(min(xi_1i[i],n_cluster-1))) {
        if(l>1){
          if (l<xi_1i[i]){
            Z_1[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_1[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }else{
          if (l<xi_1i[i]){
            Z_1[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_1[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }
      }
    }
    
    # ----------- Confounder-Dependent Weights  ------------
    
    #BETA:
    #beta_0
    gruppi=which(t_xi_0!=0)
    if (max(gruppi)==n_cluster)
      gruppi=gruppi[-length(gruppi)]
    for (l in gruppi){
      val=which(xi_0i>=l)
      z_tilde=unlist(sapply(val, function(i) Z_0[[i]][l]))
      x_tilde=matrix(matrix_X[T0[val],],ncol=n_X)
      V=solve(diag(n_X)/p_beta[2]+t(x_tilde)%*%x_tilde)
      beta_0[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,V%*%(1/p_beta[2]*(diag(n_X)%*%rep(p_beta[1],n_X))+t(x_tilde)%*%z_tilde),V)[,1:n_X]
    }
    vuoti=which(t_xi_0==0)
    if (length(vuoti)>0){
      if (max(vuoti)==n_cluster)
        vuoti=vuoti[-length(vuoti)]
    }
    for (l in vuoti){
      beta_0[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,rep(p_beta[1],n_X),diag(n_X)/p_beta[2])[,1:n_X]
    }
    #beta_1
    gruppi=which(t_xi_1!=0)
    if (max(gruppi)==n_cluster)
      gruppi=gruppi[-length(gruppi)]
    for (l in gruppi){
      val=which(xi_1i>=l)
      z_tilde=unlist(sapply(val, function(i) Z_1[[i]][l]))
      x_tilde=matrix(matrix_X[T1[val],],ncol=n_X)
      V=solve(diag(n_X)/p_beta[2]+t(x_tilde)%*%x_tilde)
      beta_1[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,V%*%(1/p_beta[2]*(diag(n_X)%*%rep(p_beta[1],n_X))+t(x_tilde)%*%z_tilde),V)[,1:n_X]
    }
    vuoti=which(t_xi_1==0)
    if (length(vuoti)>0){
      if (max(vuoti)==n_cluster)
        vuoti=vuoti[-length(vuoti)]
    }
    for (l in vuoti){
      beta_1[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,rep(p_beta[1],n_X),diag(n_X)/p_beta[2])[,1:n_X]
    }
    
    # --- imputing post-treatment missing variables ----
    
    #P_MIS:
    #t0->mis_P_1
    om_0=sapply(T0, function(i) omega(X=matrix_X[i,],beta=beta_1))
    om_0[which(is.nan(om_0))]=0
    p0=sapply(1:n0, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_0[,i])))
    P_mis[T0]=rnorm(T0,eta[p0],sigma[p0])
    #t1->mis_P_0
    om_1=sapply(T1, function(i) omega(X=matrix_X[i,],beta=beta_0))
    om_1[which(is.nan(om_1))]=0
    p1=sapply(1:n1, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_1[,i])))
    v_2=exp(lambda[1]+lambda[2]*P_obs_1)/((theta_1[3]+theta_1[4]*P_obs_1)^2)
    var_inv=1/sigma[p1]+1/v_2
    m_2=(Y_obs_1-theta_1[1]-theta_1[2]*P_obs_1)/(theta_1[3]+theta_1[4]*P_obs_1)
    P_mis[T1]=rnorm(T1,1/(var_inv)*(eta[p1]/sigma[p1]+m_2/v_2),1/(var_inv))
    
    # --- imputing post-treatment observed variables ----
    
    #P_IMP:
    #t=0
    #om_0=sapply(T0, function(i) omega(X=matrix_X[i,],beta=beta_0))
    #om_0[which(is.nan(om_0))]=0
    #p0__=sapply(1:n0, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_0[,i])))
    #P_imp[T0]=rnorm(T0,eta[p0__],sigma[p0__])
    #t=1
    #om_1=sapply(T1, function(i) omega(X=matrix_X[i,],beta=beta_1))
    #om_1[which(is.nan(om_1))]=0
    #p1__=sapply(1:n1, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_1[,i])))
    #P_imp[T1]=rnorm(T1,eta[p1__],sigma[p1__])
    
    # --- Y model ----
    
    #means regression
    P_tilde=cbind(1,P_obs_0)
    M=t(P_tilde)%*%diag(1/rep((exp(lambda[1]))^2,n0))%*%t(t(Y_obs_0))+
      p_theta[1]/(p_theta[2]^2)
    V=t(P_tilde)%*%diag(1/rep((exp(lambda[1]))^2,n0))%*%P_tilde+
      diag(rep(1/(p_theta[2]^2),2))
    V[lower.tri(V)] = t(V)[lower.tri(V)]
    theta_0=c(rmvnorm(1,solve(V)%*%M,solve(V)))
    
    P_tilde=cbind(1,P_obs_1,P_mis[T1],P_mis[T1]*P_obs_1)
    M=t(P_tilde)%*%diag(1/((exp(lambda[1]+lambda[2]*P_obs_1))^2))%*%t(t(Y_obs_1))+
      p_theta[1]/(p_theta[2]^2)
    V=t(P_tilde)%*%diag(1/((exp(lambda[1]+lambda[2]*P_obs_1))^2))%*%P_tilde+
      diag(rep(1/(p_theta[2]^2),4))
    V[lower.tri(V)] = t(V)[lower.tri(V)]
    theta_1=c(rmvnorm(1,solve(V)%*%M,solve(V)))
    
    # variances
    l0_star=rnorm(1,p_lambda[1],p_lambda[2])  
    l1_star=rnorm(1,p_lambda[3],p_lambda[4])  
    
    mean_M=c(cbind(1,P_obs_0)%*%theta_0,
             P_tilde%*%theta_1)
    var_0=c(rep(exp(lambda[1]),n0),exp(lambda[1]+lambda[2]*P_obs_1))
    var_star_0=c(rep(exp(l0_star),n0),exp(l0_star+lambda[2]*P_obs_1))
    
    a0=min(1, exp(sum(dnorm(c(Y_obs_0,Y_obs_1),mean_M,var_star_0,log=T))-
                    sum(dnorm(c(Y_obs_0,Y_obs_1),mean_M,var_0,log=T))))
    if (runif(1)<a0){
      lambda[1]=l0_star
      acc_0=acc_0+1
    } 
    
    var_star_1=exp(lambda[1]+l1_star*P_obs_1)
    
    a1=min(1, exp(sum(dnorm(Y_obs_1,mean_M[(n0+1):n],var_star_1,log=T))-
                    sum(dnorm(Y_obs_1,mean_M[(n0+1):n],var_0[(n0+1):n],log=T))))
    if (runif(1)<a1){
      lambda[2]=l1_star
      acc_1=acc_1+1
    } 
    
    #posteriors
    if(r>R_burnin){
      # --- imputing outcome ----
      P_0_r=rep(NA,n)
      P_0_r[T0]=P_obs_0
      P_0_r[T1]=P_mis[T1]
      Y_0_imp=rnorm(n,cbind(1,P_0_r)%*%t(t(theta_0)),
                    exp(lambda[1]))
      P_1_r=rep(NA,n)
      P_1_r[T1]=P_obs_1
      P_1_r[T0]=P_mis[T0]
      Y_1_imp=rnorm(n,cbind(1,P_1_r,P_0_r,P_0_r*P_1_r)%*%t(t(theta_1)),
                    exp(lambda[1]+lambda[2]*P_1_r))
      
      # P saved
      P_0_r[T0]=P_imp[T0]
      P_1_r[T1]=P_imp[T1]
      
      # --- save info ----
      post_eta[,r-R_burnin]=eta
      post_var[,r-R_burnin]=sigma
      post_beta[,r-R_burnin]=c(beta_0,beta_1)
      post_theta[,r-R_burnin]=c(theta_0,theta_1)
      post_lambda[,r-R_burnin]=lambda
      post_xi_0[T0,r-R_burnin]=xi_0i
      post_xi_0[T1,r-R_burnin]=p1
      post_xi_1[T1,r-R_burnin]=xi_1i
      post_xi_1[T0,r-R_burnin]=p0
      post_P_0[,r-R_burnin]=P_0_r
      post_P_1[,r-R_burnin]=P_1_r
      post_Y_0_imp[,r-R_burnin]=Y_0_imp
      post_Y_1_imp[,r-R_burnin]=Y_1_imp
    }
    
    #print(r)
    if (r%%500==0) print(r)
  }
  
  #partitions
  S_all_cluster=partition_all_cluster(xi_0=post_xi_0,xi_1=post_xi_1)
  S_strata_cluster=partition_strata_cluster(xi_0=post_xi_0,xi_1=post_xi_1, eta=post_eta)
  
  print(paste0("campione n ",c))
  
  return(list(post_eta=post_eta, post_var=post_var,
              post_beta=post_beta,
              S_all_cluster=S_all_cluster, S_strata_cluster=S_strata_cluster,
              post_theta=post_theta, post_lambda=post_lambda,
              acc_rate_0=acc_0/R, acc_rate_1=acc_1/R,
              post_xi_0=post_xi_0, post_xi_1=post_xi_1,
              post_P_0=post_P_0, post_P_1=post_P_1,
              post_Y_0_imp=post_Y_0_imp, post_Y_1_imp=post_Y_1_imp
  ))
}  

###################################################################

# Y-model with same covariates in the weights
matrix_REGR=matrix_X

# max number clusters
n_cluster=10

start.time <- Sys.time()
model_data_REGweights_10=BNP_model_REG(c=1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# max number clusters
n_cluster=12

start.time <- Sys.time()
model_data_REGweights_12=BNP_model_REG(c=1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

###################################################################

# Y-model with ALL covariates
matrix_REGR=dataset_matched[,c("PctBlack","PctHisp","PctHighSchool",
                               "PctUrban","PctFemale","PctPoor","MedianHHInc",
                               "PctOccupied","PctMovedIn5","MedianHValue",
                               "smokerate2000","avgdewpt","avgtemp",
                               "avgrelhum","Population")]
matrix_REGR=as.matrix(matrix_REGR)

start.time <- Sys.time()
model_data_REGall=BNP_model_REG(c=1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

###################################################################

# Y-model with NO covariates
matrix_REGR=NA

start.time <- Sys.time()
model_data_NO=BNP_model_NO(c=1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

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
