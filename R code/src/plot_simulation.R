###################################################################################
# ---             SIMULATION STUDY:             ---
# ---                5 settings                 ---
# ---                  plots                    ---
###################################################################################

#library
library(mvtnorm)
library(ggplot2)
library(ggExtra)

###################################################################################

# POST_TREATMENT variable: divider by treatment level

hist_P_sim<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                P=(data$data$T)*data$data$P_1+
                                  (1-data$data$T)*data$data$P_0))
  
  #pdf(file=paste0("sim_P_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(P, fill=as.factor(T))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n level",
                      labels = c("0", "1"),
                      values = c("#E69F00", "#56B4E9"))+
    theme_bw()+
    geom_vline(xintercept=data$par_P$eta, color = "black", 
               size=0.4,linetype = "twodash")+
    labs(x = "Post-Treatment Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

# POST_TREATMENT variable: divider by treatment level

hist_diff_P_sim<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                P=data$data$P_1-data$data$P_0))
  
  #pdf(file=paste0("sim_P_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(P)) + 
    geom_histogram(alpha = 0.4, position="identity", fill="#56B4E9")+
    theme_bw()+
    geom_vline(xintercept=data$par_P$eta[data$par_P$allocation_1]-
                 data$par_P$eta[data$par_P$allocation_0], 
               color = "black", size=0.4,linetype = "twodash")+
    labs(x = "Post-Treatment Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

###################################################################################

# OUTCOME variable: divider by treatment level

hist_Y_sim<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                P=(data$data$T)*data$data$Y_1+
                                  (1-data$data$T)*data$data$Y_0))
  
  #pdf(file=paste0("sim_Y_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(P, fill=as.factor(T))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n level",
                      labels = c("0", "1"),
                      values = c("#E69F00", "#56B4E9"))+
    theme_bw()+
    labs(x = "Outcome Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}

# POST_TREATMENT variable: divider by treatment level

hist_diff_Y_sim<-function(data,s_n){
  
  data_hist=as.data.frame(cbind(T=data$data$T,
                                P=data$data$Y_1-data$data$Y_0))
  
  #pdf(file=paste0("sim_Y_",s_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(P)) + 
    geom_histogram(alpha = 0.4, position="identity", fill="#56B4E9")+
    theme_bw()+
    labs(x = "Post-Treatment Var.",y = " ")+
    ggtitle(paste0("Scenario ",s_n))
  print(g)
  #dev.off()
}
