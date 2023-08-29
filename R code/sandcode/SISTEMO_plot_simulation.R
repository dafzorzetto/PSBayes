###################################################################################
# ---             SIMULATION STUDY:             ---
# ---                5 settings                 ---
# ---  CONFONDER-DEPENDENT SHARED ATOMS MODEL   ---
###################################################################################

#library
library(mvtnorm)
library(ggplot2)
library(ggExtra)

###################################################################################

# --- generation functions ---

n=500          # units for each samples
samples=100    # number of samples

# general function:
setup_sim<- function(seed,eta,sigma_p,allocation_0,allocation_1,beta_0,beta_1,sigma_y){
  
  # set seed for riproducibility
  set.seed(seed)
  
  # covariates
  X=cbind(rbinom(n,1,0.7),rbinom(n,1,0.6))
  
  # treatment
  reg_T=0+0.4*(X[,1])+0.3*(X[,2])
  logit_T=exp(reg_T)/(1+exp(reg_T))
  T=rbinom(n,1,logit_T)
  
  # cluster allocation
  S_dummy=cbind(apply(X,1,prod),
                (X[,1]==0),
                (X[,1]==1 & X[,2]==0))
  S_cl_0=S_dummy%*%allocation_0
  S_cl_1=S_dummy%*%allocation_1
  
  # outcomes simulation
  P=sapply(1:n, function(i) rmvnorm(1,c(eta[S_cl_0[i]],eta[S_cl_1[i]]),
                                    c(sigma_p[S_cl_0[i]],sigma_p[S_cl_1[i]])*diag(2)))
  
  Y=cbind(rnorm(n,beta_0[1]+beta_0[2]*P[1,],exp(sigma_y[1])),
          rnorm(n,beta_1[1]+beta_1[2]*P[2,]+beta_1[3]*P[1,]+beta_1[4]*P[1,]*P[2,],
                exp(sigma_y[1]+sigma_y[2]*P[2,])))
  
  # saving groups and strata allocation
  S_groups=S_cl_1*(S_cl_0+5)
  S_strata=1*(eta[S_cl_1]==eta[S_cl_0])+2*(eta[S_cl_1]>eta[S_cl_0])
  
  return(list(data=list(X=X,T=T,P_0=P[1,],P_1=P[2,],Y_0=Y[,1],Y_1=Y[,2]), 
              par_P=cbind(eta=eta, sigma_p=sigma_p, 
                          allocation_0=allocation_0, allocation_1=allocation_1),
              par_Y=list(beta_0=beta_0, beta_1=beta_1, sigma_y=sigma_y),
              clusters=list(S_groups=S_groups,S_strata=S_strata)))
}


# --- 5 settings ---

# 3 atoms - 1 DIS, 1 ASS+                                                                                                                                                          in common - allocation is shifted
scenario_1=lapply(1:samples, function(s) 
  prova=setup_sim(seed=s,
                  eta=c(1,2,3),
                  sigma_p=rep(0.05,3),
                  allocation_0=c(1,2,2),
                  allocation_1=c(1,3,3),
                  beta_0=c(1,2),
                  beta_1=c(1,2,-1,1.5),
                  sigma_y=c(-0.5,0.1)))

# 3 atoms - 1 DISS, 2 ASS+ with different effects
scenario_2=lapply(1:samples, function(s) 
  prova=setup_sim(seed=s,
                  eta=c(0.5,2,3),
                  sigma_p=rep(0.05,3),
                  allocation_0=c(1,2,3),
                  allocation_1=c(2,3,3),
                  beta_0=c(1,2),
                  beta_1=c(1,2,-1,1.5),
                  sigma_y=c(-0.5,0.1)))

# 3 atoms - 1 DISS, 1 ASS+, 1 ASS- 
scenario_3=lapply(1:samples, function(s) 
  prova=setup_sim(seed=s,
                  eta=c(1,2,3),
                  sigma_p=rep(0.05,3),
                  allocation_0=c(2,2,2),
                  allocation_1=c(1,2,3),
                  beta_0=c(1,2),
                  beta_1=c(1,2,-1,1.5),
                  sigma_y=c(-0.5,0.1)))

# 3 atoms - case2 closer
scenario_4=lapply(1:samples, function(s) 
  prova=setup_sim(seed=s,
                  eta=c(1.5,2,2.75),
                  sigma_p=c(0.2,0.1,0.15),
                  allocation_0=c(1,2,3),
                  allocation_1=c(2,3,3),
                  beta_0=c(1,2),
                  beta_1=c(1,2,-1,1.5),
                  sigma_y=c(-0.5,0.1)))

# 3 atoms - case3 closer
scenario_5=lapply(1:samples, function(s) 
  prova=setup_sim(seed=s,
                  eta=c(1.5,2,2.5),
                  sigma_p=c(0.2,0.1,0.15),
                  allocation_0=c(2,2,2),
                  allocation_1=c(1,2,3),
                  beta_0=c(1,2),
                  beta_1=c(1,2,-1,1.5),
                  sigma_y=c(-0.5,0.1)))

save.image("data_simulations.RData")

###################################################################################
# ---     plots: histograms of post-treatment or outcome distributions     ---
###################################################################################

# --- joint distribution of post treatment variable ----

plot_sim<-function(sim){
  
  sp <- ggplot(as.data.frame(sim$data), aes(x=P_0, y=P_1, colour=as.factor(T), shape = as.factor(T)))+
    geom_point(aes(shape=as.factor(T), color = as.factor(T)))+
    #scale_color_manual(values=c("#E69F00", "#56B4E9"))+
    scale_colour_manual(name = " ",
                        labels = c("Control", "Treatment"),
                        values = c("#E69F00", "#56B4E9")) +  
    #values = c("#B53636", "#2F9C2F")) +  
    scale_shape_manual(name = " ",
                       labels = c("Control", "Treatment"),
                       values = c(17, 19))+
    theme_bw()+
    geom_vline(xintercept=sim$par[unique(sim$par[,3]),1], color = "black", 
               size=0.4,linetype = "twodash")+
    geom_hline(yintercept=sim$par[unique(sim$par[,4]),1], color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "P(0)",y = "P(1)",col=" ", shape="Treatment")+
    theme(legend.position="bottom")
  
  ggMarginal(sp, type="histogram",xparams = list(  bins=60),yparams = list(  bins=60),
             color="white",fill="grey60", size=3)
  
}

plot_sim(sim=scenario_1[[1]])
plot_sim(sim=scenario_2[[1]])
plot_sim(sim=scenario_3[[1]])
plot_sim(sim=scenario_4[[1]])
plot_sim(sim=scenario_5[[1]])


# histograms

hist_sim<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$T,P=(data$T)*data$P_1+(1-data$T)*data$P_0))
  
  #setwd("C:/Users/dafne/Desktop/plot temporanei/princ_strata/CA")
  #pdf(file=paste0("sim_P_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(P, fill=as.factor(T))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n level",
                      labels = c("0", "1"),
                      values = c("#E69F00", "#56B4E9"))+
    theme_bw()+
    geom_vline(xintercept=data$par[,1], color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "Post-Treatment Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

hist_sim(data=scenario_1[[1]]$data,s_n=1)
hist_sim(data=scenario_2[[1]]$data,s_n=2)
hist_sim(data=scenario_3[[1]]$data,s_n=3)
hist_sim(data=scenario_4[[1]]$data,s_n=4)
hist_sim(data=scenario_5[[1]]$data,s_n=5)


# --- joint distribution of outcome variable ----

hist_sim_Y<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$T,Y=(data$T)*data$Y_1+(1-data$T)*data$Y_0))
  
  #setwd("C:/Users/dafne/Desktop/plot temporanei/princ_strata/CA")
  #pdf(file=paste0("sim_Y_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(Y, fill=as.factor(T))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n level",
                      labels = c("0", "1"),
                      values = c("#E69F00", "#56B4E9"))+
    theme_bw()+
    geom_vline(xintercept=data$par[,1], color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "Outcome Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

hist_sim_Y(data=scenario_1[[1]]$data,s_n=1)
hist_sim_Y(data=scenario_2[[1]]$data,s_n=2)
hist_sim_Y(data=scenario_3[[1]]$data,s_n=3)
hist_sim_Y(data=scenario_4[[1]]$data,s_n=4)
hist_sim_Y(data=scenario_5[[1]]$data,s_n=5)


# --- joint distribution of outcome variable -- cluster allocation ----

hist_sim_Y_cl<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                Y=(data$data$T)*data$data$Y_1+
                                  (1-data$data$T)*data$data$Y_0,
                                C=(data$data$T)*data$clusters_all[,1]*4+
                                  (1-data$data$T)*data$clusters_all[,1]))
  
  #setwd("C:/Users/dafne/Desktop/plot temporanei/princ_strata/CA")
  #pdf(file=paste0("sim_Y_cl_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(Y, fill=as.factor(C))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n - cluster",
                      labels = c("0-cl1","0-cl2", "1-cl1","1-cl2"),
                      values = c("#FF9B0D","#E69F00","#56B4E9","#0DA9FF"))+
    theme_bw()+
    geom_vline(xintercept=data$par[,1], color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "Outcome Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

hist_sim_Y_cl2<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                Y=(data$data$T)*data$data$Y_1+
                                  (1-data$data$T)*data$data$Y_0,
                                C=(data$data$T)*data$clusters_all[,1]*4+
                                  (1-data$data$T)*data$clusters_all[,1]))
  
  #pdf(file=paste0("sim_Y_cl_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(Y, fill=as.factor(C))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n - cluster",
                      labels = c("0-cl1","0-cl2","0-cl3","1-cl1","1-cl2","1-cl3"),
                      values = c("#FF4A00","#FF9B0D","#E69F00","#56B4E9","#0DA9FF","#3966B3"))+
    theme_bw()+
    geom_vline(xintercept=data$par[,1], color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "Outcome Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

hist_sim_Y_cl3<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                Y=(data$data$T)*data$data$Y_1+
                                  (1-data$data$T)*data$data$Y_0,
                                C=(data$data$T)*data$clusters_all[,2]*4+
                                  (1-data$data$T)*data$clusters_all[,2]))
  
  #pdf(file=paste0("sim_Y_cl_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(Y, fill=as.factor(C))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n - cluster",
                      labels = c("0-cl1","0-cl2","0-cl3","1-cl1","1-cl2","1-cl3"),
                      values = c("#FF4A00","#FF9B0D","#E69F00","#56B4E9","#0DA9FF","#3966B3"))+
    theme_bw()+
    geom_vline(xintercept=data$par[,1], color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "Outcome Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

hist_sim_Y_cl(data=sim1[[1]],s_n=1)
hist_sim_Y_cl2(data=scenario_2[[1]],s_n=2)
hist_sim_Y_cl3(data=sim3[[1]],s_n=3)
hist_sim_Y_cl2(data=sim4[[1]],s_n=4)
hist_sim_Y_cl3(data=sim5[[1]],s_n=5)


##############################################################
##############################################################
# ---    boxplot ----
##############################################################

#libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)

cbPalette <- c("#35D90B", "#F0D400", "#D90224")

boxplot_P<-function(data_frame,sim_n){
  
  bp <- ggplot(data_frame, aes(x=P, y=cl, fill=cl)) + 
    scale_fill_manual(values=cbPalette, name="")+
    geom_boxplot()+
    geom_vline(xintercept = 0, col="#00ace6", size=1) +
    theme(legend.position = "none",
          panel.background = element_rect(fill='white'),
          plot.background = element_rect(fill ="white"),
          #panel.grid.minor = element_line(color = "grey"),
          axis.title = element_text(size=23, face="bold"),
          #legend.text=element_text(size=14),
          plot.title = element_text(hjust = 0.5),
          title =element_text(size=15),
          #legend.background = element_rect(fill='transparent'),
          #panel.grid.major = element_line(color = "grey",size = 0.45)
          axis.text.y = element_text(size=15),
    )+
    #geom_hline(yintercept = vero, color = "#0BC2CF", size=0.65)+
    ylab(paste0("scenario ",sim_n,"\n")) +
    xlab("") +
    coord_cartesian(xlim = c(-2, 2)) +
  ggtitle("P(1)-P(0)|stratum")
  
  #pdf(file=paste0(sim_n, "_P.pdf"),width=10, height=5)
  print(bp)
  #dev.off()
  
}

data_frame=as.data.frame(cbind(P=scenario_1[[1]]$data$P_1-scenario_1[[1]]$data$P_0,
                               cl=c(scenario_1[[1]]$clusters$S_groups)))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==6]="0"
data_frame$cl=factor(data_frame$cl, levels=c("0", "+1"))

cbPalette <- c("#F0D400", "#D90224")
boxplot_P(data_frame,sim_n=1)

data_frame=as.data.frame(cbind(P=scenario_3[[1]]$data$P_1-scenario_3[[1]]$data$P_0,
                               cl=c(scenario_3[[1]]$clusters$S_groups)))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==14]="0"
data_frame$cl[data_frame$cl==7]="-1"
data_frame$cl=factor(data_frame$cl, levels=c("-1","0", "+1"))

cbPalette <- c( "#35D90B", "#F0D400", "#D90224")
boxplot_P(data_frame,sim_n=3)

data_frame=as.data.frame(cbind(P=scenario_4[[1]]$data$P_1-scenario_4[[1]]$data$P_0,
                               cl=c(scenario_4[[1]]$clusters$S_groups)))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==24]="0"
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==12]="+1."
data_frame$cl=factor(data_frame$cl, levels=c("0", "+1.", "+1"))

cbPalette <- c( "#F0D400", "#D90224", "#D90224")
boxplot_P(data_frame,sim_n=4)

data_frame=as.data.frame(cbind(P=scenario_5[[1]]$data$P_1-scenario_5[[1]]$data$P_0,
                               cl=c(scenario_5[[1]]$clusters$S_groups)))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==14]="0"
data_frame$cl[data_frame$cl==7]="-1"
data_frame$cl=factor(data_frame$cl, levels=c("-1","0", "+1"))

cbPalette <- c( "#35D90B", "#F0D400", "#D90224")
boxplot_P(data_frame,sim_n=5)

#########################################################

boxplot_Y<-function(data_frame,sim_n){
  
  bp <- ggplot(data_frame, aes(x=Y, y=cl, fill=cl)) + 
    scale_fill_manual(values=cbPalette, name="")+
    geom_boxplot()+
    geom_vline(xintercept = 0, col="#00ace6", size=1) +
    theme(legend.position = "none",
          panel.background = element_rect(fill='white'),
          plot.background = element_rect(fill ="white"),
          #panel.grid.minor = element_line(color = "grey"),
          axis.title = element_text(size=23, face="bold"),
          #legend.text=element_text(size=14),
          plot.title = element_text(hjust = 0.5),
          title =element_text(size=15),
          #legend.background = element_rect(fill='transparent'),
          #panel.grid.major = element_line(color = "grey",size = 0.45)
          axis.text.y = element_text(size=15),
    )+
    #geom_hline(yintercept = vero, color = "#0BC2CF", size=0.65)+
    ylab(paste0("scenario ",sim_n,"\n")) +
    xlab("") +
    #coord_cartesian(xlim = c(-2, 2)) +
    ggtitle("Y(1)-Y(0)|stratum")
  
  #pdf(file=paste0(sim_n, "_Y.pdf"),width=10, height=5)
  print(bp)
  #dev.off()
  
}

data_frame=as.data.frame(cbind(Y=scenario_1[[1]]$data$Y_1-scenario_1[[1]]$data$Y_0,
                               cl=c(scenario_1[[1]]$clusters$S_groups)))
data_frame$Y=as.numeric(data_frame$Y)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==6]="0"
data_frame$cl=factor(data_frame$cl, levels=c("0", "+1"))

cbPalette <- c("#F0D400", "#D90224")
boxplot_Y(data_frame,sim_n=1)

data_frame=as.data.frame(cbind(Y=scenario_2[[1]]$data$Y_1-scenario_2[[1]]$data$Y_0,
                               cl=c(scenario_2[[1]]$clusters$S_groups)))
data_frame$Y=as.numeric(data_frame$Y)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==24]="0"
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==12]="+1."
data_frame$cl=factor(data_frame$cl, levels=c("0", "+1", "+1."))

cbPalette <- c("#F0D400", "#D90224", "#D90224")
boxplot_Y(data_frame,sim_n=2)

data_frame=as.data.frame(cbind(Y=scenario_3[[1]]$data$Y_1-scenario_3[[1]]$data$Y_0,
                               cl=c(scenario_3[[1]]$clusters$S_groups)))
data_frame$Y=as.numeric(data_frame$Y)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==14]="0"
data_frame$cl[data_frame$cl==7]="-1"
data_frame$cl=factor(data_frame$cl, levels=c("-1","0", "+1"))

cbPalette <- c( "#35D90B", "#F0D400", "#D90224")
boxplot_Y(data_frame,sim_n=3)

data_frame=as.data.frame(cbind(Y=scenario_4[[1]]$data$Y_1-scenario_4[[1]]$data$Y_0,
                               cl=c(scenario_4[[1]]$clusters$S_groups)))
data_frame$Y=as.numeric(data_frame$Y)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==24]="0"
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==12]="+1."
data_frame$cl=factor(data_frame$cl, levels=c("0", "+1", "+1."))

cbPalette <- c( "#F0D400", "#D90224", "#D90224")
boxplot_Y(data_frame,sim_n=4)

data_frame=as.data.frame(cbind(Y=scenario_5[[1]]$data$Y_1-scenario_5[[1]]$data$Y_0,
                               cl=c(scenario_5[[1]]$clusters$S_groups)))
data_frame$Y=as.numeric(data_frame$Y)
data_frame$cl=as.character(data_frame$cl)
data_frame$cl[data_frame$cl==21]="+1"
data_frame$cl[data_frame$cl==14]="0"
data_frame$cl[data_frame$cl==7]="-1"
data_frame$cl=factor(data_frame$cl, levels=c("-1","0", "+1"))

cbPalette <- c( "#35D90B", "#F0D400", "#D90224")
boxplot_Y(data_frame,sim_n=5)


