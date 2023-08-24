###################################################################################
# ---             ESTIMATION MODEL:             ---
# ---                5 settings                 ---
# ---  CONFONDER-DEPENDENT SHARED ATOMS MODEL   ---
# ---           WITH Fasano e Durante           ---
###################################################################################

#library
library(mvtnorm)
library(CholWishart)
library(parallel)
library(truncnorm)
library(invgamma)
library(BNPmix)

load("data_simulations.RData")

###################################################################################

# iterations
R=3000
R_burnin=1500

# max number clusters
n_cluster=7

###################################################################################

All_model<-function(c,sim){
  
  # set seed for riproducibility
  set.seed(c)
  
  # ------   prearing variables   ------
  
  # main variables
  matrix_X=cbind(rep(1,n),sim[[c]]$data$X)
  T_var=sim[[c]]$data$T
  P_obs=T_var*sim[[c]]$data$P_1+(1-T_var)*sim[[c]]$data$P_0
  Y_obs=T_var*sim[[c]]$data$Y_1+(1-T_var)*sim[[c]]$data$Y_0
  
  # number covariates (X)
  n_X=dim(matrix_X)[2]
  
  # dividing observation in the treatment levels
  # level t=1
  T1=which(T_var==1)
  n1=length(T1)
  P_obs_1=P_obs[T1]
  Y_obs_1=Y_obs[T1]
  # level t=0
  T0=which(T_var==0)
  n0=length(T0)
  P_obs_0=P_obs[T0]
  Y_obs_0=Y_obs[T0]
  
  
  # preparing vector for P imputation
  P_mis=rep(0,n)
  P_imp=rep(0,n)
  
  # ------   hyperparameters   -----
  p_beta=c(0,5)
  p_sigma=c(0.5,0.25)
  p_eta=c(0,5)
  p_sun_xi=0
  sd_sun=10
  q_sun=n_X*(n_cluster-1)
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
  
  partizione_function<-function(xi_0,xi_1){
    grid_couple=expand.grid(1:n_cluster,1:n_cluster)
    
    #cartesian product of cluster allocation
    part=sapply(1:n, function(i)
      sapply(1:(R-R_burnin), function(r) 
        which(grid_couple$Var1==xi_0[i,r] & grid_couple$Var2==xi_1[i,r])))
    
    #point estimation
    WG=partition.BNPdens(list(clust=part),dist = "VI")$partitions[1,]
    tab=cbind(WG)
    dimnames(tab)=NULL
    return(tab)
  }
  
  posterior_eta<-function(partizione){
    stime=list(stime_0=rep(NA,max(partizione)),stime_1=rep(NA,max(partizione)))
    for (g in 1:max(partizione)){
      
      #units in each group
      unita_0=intersect(T0,which(partizione==g))
      unita_1=intersect(T1,which(partizione==g))
      
      #estimation 
      stime$stime_0[g]=mean(sapply(1:length(unita_0), function(x) 
        mean(sapply(1:(R-R_burnin), function(j) 
          post_eta[post_xi_0[unita_0[x],j],R_burnin+j]), na.rm=FALSE)))
      
      stime$stime_1[g]=mean(sapply(1:length(unita_1), function(x) 
        mean(sapply(1:(R-R_burnin), function(j) 
          post_eta[post_xi_1[unita_1[x],j],R_burnin+j]), na.rm=FALSE)))
    }
    return(stime)
  }
  
  xi_0i=sample(1:n_cluster,n0,replace=TRUE)
  xi_1i=sample(1:n_cluster,n1,replace=TRUE)
  t_xi_0=table(xi_0i)
  t_xi_1=table(xi_1i)
  eta_sogg_0=eta[xi_0i]
  eta_sogg_1=eta[xi_1i]
  
  # acceptation rate for the metropolis proposal step (variance Y-model)
  acc_0=0
  acc_1=0
  
  #chains:
  post_eta=matrix(NA,ncol=R,nrow=n_cluster)
  post_var=matrix(NA,ncol=R,nrow=n_cluster)
  post_xi_0=matrix(NA,ncol=R-R_burnin,nrow=n)     
  post_xi_1=matrix(NA,ncol=R-R_burnin,nrow=n)     
  post_beta=matrix(NA,ncol=R, nrow=length(beta_0)+length(beta_1))
  post_theta=matrix(NA,ncol=R,nrow=6)
  post_lambda=matrix(NA,ncol=R, nrow=2)
  
  post_P_mis=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_P_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  
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
    
    # ----------- Fasano and Durante (2022)  ------------
    n_bar_0_star=xi_0i
    n_bar_0_star[n_bar_0_star==n_cluster]=n_cluster-1
    n_bar_0=sum(n_bar_0_star)
    n_bar_1_star=xi_1i
    n_bar_1_star[n_bar_1_star==n_cluster]=n_cluster-1
    n_bar_1=sum(n_bar_1_star)
    
    # variable for identify blocks
    cum_n_bar_0=cumsum(c(0,n_bar_0_star))+1
    cum_n_bar_1=cumsum(c(0,n_bar_1_star))+1
    
    bar_X_0=matrix(0, nrow=n_bar_0, ncol=q_sun)
    for(i in 1:n0){
      bar_X_0[cum_n_bar_0[i]:(cum_n_bar_0[i+1]-1),1:(n_bar_0_star[i]*n_X)]=
        diag(c(rep(-1,xi_0i[i]-1),rep(1,xi_0i[i]<n_cluster)))%x%t(matrix_X[T0[i],])
    }
    bar_X_1=matrix(0, nrow=n_bar_1, ncol=q_sun)
    for(i in 1:n1){
      bar_X_1[cum_n_bar_1[i]:(cum_n_bar_1[i+1]-1),1:(n_bar_1_star[i]*n_X)]=
        diag(c(rep(-1,xi_1i[i]-1),rep(1,xi_1i[i]<n_cluster)))%x%t(matrix_X[T1[i],])
    }
    
    d_post_0=diag(sqrt(sd_sun^2*diag(bar_X_0%*%t(bar_X_0))+1))
    d_post_1=diag(sqrt(sd_sun^2*diag(bar_X_1%*%t(bar_X_1))+1))
    
    Delta_0=sd_sun*t(bar_X_0)%*%solve(d_post_0)
    Delta_1=sd_sun*t(bar_X_1)%*%solve(d_post_1)
    
    Gamma_0=solve(d_post_0)%*%(sd_sun^2*bar_X_0%*%t(bar_X_0)+diag(n_bar_0))%*%solve(d_post_0)
    Gamma_1=solve(d_post_1)%*%(sd_sun^2*bar_X_1%*%t(bar_X_1)+diag(n_bar_1))%*%solve(d_post_1)
    
    gamma_pst_0=solve(d_post_0)%*%bar_X_0%*%rep(p_sun_xi,q_sun)
    gamma_pst_1=solve(d_post_1)%*%bar_X_1%*%rep(p_sun_xi,q_sun)
    
    # normal_{p+1}
    B_0pst_0=rmvnorm(1,rep(0,q_sun),round(diag(q_sun)-Delta_0%*%solve(Gamma_0)%*%t(Delta_0),8))
    B_0pst_1=rmvnorm(1,rep(0,q_sun),round(diag(q_sun)-Delta_1%*%solve(Gamma_1)%*%t(Delta_1),8))
    
    #  truncated normal multivariate
    # the Gamma matrix is positive defined but the deterinant is zero!
    # is.positive.definite(Gamma_0)  is.positive.definite(Gamma_1)     #TRUE
    # det(Gamma_0)                   det(Gamma_1)                   #0
    # then the rtmvnorm doesn't work
    # we have to add a very small variance in the diagonal so the determinant is different to zero
    
    while(det(Gamma_0)==0){
      Gamma_0=Gamma_0+diag(ncol(Gamma_0))*0.03   #abs(rnorm(ncol(Gamma_0),0.05,0.05))
    }
    while(det(Gamma_1)==0){
      Gamma_1=Gamma_1+diag(ncol(Gamma_1))*0.03   #abs(rnorm(ncol(Gamma_1),0.05,0.05))
    }
    
    #B_1pst_0=rtmvnorm(1,mean=rep(0,nrow(Gamma_0)),sigma=Gamma_0,lower=-c(gamma_pst_0),
    #                  algorithm="gibbsR")
    #B_1pst_1=rtmvnorm(1,mean=rep(0,nrow(Gamma_1)),sigma=Gamma_1,lower=-c(gamma_pst_1),
    #                  algorithm="gibbsR")
    
    # Extremelly tooo long!!!
    
    B_1pst_0=abs(rmvnorm(1,rep(0,nrow(Gamma_0)),Gamma_0))
    B_1pst_1=abs(rmvnorm(1,rep(0,nrow(Gamma_1)),Gamma_1))
    
    # beta - probit parameters
    beta_0=c(p_sun_xi+sd_sun*(t(B_0pst_0)+Delta_0%*%solve(Gamma_0)%*%t(B_1pst_0)))
    beta_1=c(p_sun_xi+sd_sun*(t(B_0pst_1)+Delta_1%*%solve(Gamma_1)%*%t(B_1pst_1)))
    
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
    P_mis[T1]=rnorm(T1,eta[p1],sigma[p1])
    
    # --- imputing post-treatment observed variables ----
    
    #P_IMP:
    #t=0
    om_0=sapply(T0, function(i) omega(X=matrix_X[i,],beta=beta_0))
    om_0[which(is.nan(om_0))]=0
    p0=sapply(1:n0, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_0[,i])))
    P_imp[T0]=rnorm(T0,eta[p0],sigma[p0])
    #t=1
    om_1=sapply(T1, function(i) omega(X=matrix_X[i,],beta=beta_1))
    om_1[which(is.nan(om_1))]=0
    p1=sapply(1:n1, function(i) (1:n_cluster)%*%(rmultinom(1,1,om_1[,i])))
    P_imp[T1]=rnorm(T1,eta[p1],sigma[p1])
    
    # --- Y model ----
    
    #p0=sapply(1:n0, function(i) which.max(om_0[,i]))
    #P_mis[T0]=rnorm(T0,eta[p0],sigma[p0])
    #p1=sapply(1:n1, function(i) which.max(om_1[,i]))
    #P_mis[T1]=rnorm(T1,eta[p1],sigma[p1])
    
    #P_mis=T_var*sim[[c]]$data$P_0+(1-T_var)*sim[[c]]$data$P_1
    
    #means regression
    P_tilde=cbind(1,P_obs_0)
    M=t(P_tilde)%*%diag(1/rep((exp(lambda[1]))^2,n0))%*%t(t(Y_obs_0))+p_theta[1]/(p_theta[2]^2)
    V=t(P_tilde)%*%diag(1/rep((exp(lambda[1]))^2,n0))%*%P_tilde+diag(rep(1/(p_theta[2]^2),2))
    theta_0=c(rmvnorm(1,solve(V)%*%M,solve(V)))
    
    P_tilde=cbind(1,P_obs_1,P_mis[T1],P_mis[T1]*P_obs_1)
    M=t(P_tilde)%*%diag(1/((exp(lambda[1]+lambda[2]*P_obs_1))^2))%*%t(t(Y_obs_1))+p_theta[1]/(p_theta[2]^2)
    V=t(P_tilde)%*%diag(1/((exp(lambda[1]+lambda[2]*P_obs_1))^2))%*%P_tilde+diag(rep(1/(p_theta[2]^2),4))
    V[lower.tri(V)] = t(V)[lower.tri(V)]
    theta_1=c(rmvnorm(1,solve(V)%*%M,solve(V)))
      
    # variances
    l0_star=rnorm(1,p_lambda[1],p_lambda[2])  
    l1_star=rnorm(1,p_lambda[3],p_lambda[4])  
    
    mean_M=c(cbind(1,P_obs_0)%*%theta_0,P_tilde%*%theta_1)
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
    
    # --- save info ----
    
    #posteriors
    post_eta[,r]=eta
    post_var[,r]=sigma
    post_beta[,r]=c(beta_0,beta_1)
    post_theta[,r]=c(theta_0,theta_1)
    post_lambda[,r]=lambda
    if(r>R_burnin){
      post_xi_0[T0,r-R_burnin]=xi_0i
      post_xi_0[T1,r-R_burnin]=p1
      post_xi_1[T1,r-R_burnin]=xi_1i
      post_xi_1[T0,r-R_burnin]=p0
      post_P_mis[,r-R_burnin]=P_mis
      post_P_imp[,r-R_burnin]=P_imp
    }
    
    #print(r)
    if (r%%100==0) print(r) 
  }
  
  #partizione e stime
  partizione=partizione_function(xi_0=post_xi_0,xi_1=post_xi_1)
  #stima_wg=posterior_eta(partizione)
  
  print(paste0("campione n ",c))
  
  return(list(post_eta=post_eta,
              post_beta=post_beta,
              post_var=post_var,
              partizione=partizione,
              post_theta=post_theta,
              post_lambda=post_lambda,
              acc_rate_0=acc_0/R,acc_rate_1=acc_1/R,
              post_xi_0=post_xi_0, post_xi_1=post_xi_1,
              post_P_mis=post_P_mis,post_P_imp=post_P_imp
  ))
}  


#########################################

start.time <- Sys.time()
#model_1=All_model(c=1, sim=scenario_1)
#model_1=lapply(1:3, function(i) All_model(c=i, sim=scenario_1))
model_1=mclapply(1:samples, All_model, sim=scenario_1, mc.cores = 5)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
#model_2=All_model(c=1, sim=scenario_2)
model_2=mclapply(1:samples, All_model, sim=scenario_2, mc.cores = 6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
#model_3=All_model(c=1, sim=scenario_3)
model_3=mclapply(1:samples, All_model, sim=scenario_3, mc.cores = 6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
#model_4=All_model(c=1, sim=scenario_4)
model_4=mclapply(1:samples, All_model, sim=scenario_4, mc.cores = 6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


save.image("results_allmodel.RData")