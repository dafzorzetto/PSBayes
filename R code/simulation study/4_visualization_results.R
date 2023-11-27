#########################################################################
#     ---  MODEL COMPARISON   ----
#     ---  visualize results  ----
#     ---  tables and plots   ----
#########################################################################

# load data
load("estimands.RData")

#source
source("src/visualization_results_function.R")

#libraries
library(xtable)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)

#########################################################################
 
# adjusted RAND index
xtable(ARI, digits=4)

#########################################################################

# --- Boxplots for P: post-treatment ----

cbPalette <- c("#007399", "#80dfff")   #, "#00ace6", "#80dfff")

bias_P=a<-matrix(c(apply(bias_P_CASDMM1,2,mean),apply(bias_P_SLM1,2,mean),
                   apply(bias_P_CASDMM2,2,mean),apply(bias_P_SLM2,2,mean),
                   apply(bias_P_CASDMM3,2,mean),apply(bias_P_SLM3,2,mean),
                   apply(bias_P_CASDMM4,2,mean),apply(bias_P_SLM4,2,mean),
                   apply(bias_P_CASDMM5,2,mean),apply(bias_P_SLM5,2,mean)),
                 ncol=10)

round(apply(bias_P,2,median,na.rm = TRUE), 4)
round(apply(bias_P,2,IQR,na.rm = TRUE), 4)

bias_P_boxplot=as.data.frame(cbind(Xi=c(bias_P),
                                 Q=(rep(c(rep("CASBAH",samples),rep("SLM",samples)),5)),
                                 cov=paste0("scenario ",rep(1:5,each=samples*2))))
bias_P_boxplot$cov=as.character(bias_P_boxplot$cov)
bias_P_boxplot$Q=as.character(bias_P_boxplot$Q)
bias_P_boxplot$Xi=as.numeric(bias_P_boxplot$Xi)

pdf(file="bias_P_sim.pdf",width=10, height=5)
ggplot(bias_P_boxplot, aes(x=cov, y=Xi, fill=Q)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  geom_hline(yintercept = 0, col="#D90224", size=0.4) +
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=14),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("Bias Post-Treatment") +
  xlab("")
dev.off()

#cbPalette <- c("#80dfff")   #, "#00ace6", "#80dfff")

mse_P_boxplot=as.data.frame(cbind(Xi=c(apply(mse_P_CASDMM1,2,mean),apply(mse_P_SLM1,2,mean),
                                        apply(mse_P_CASDMM2,2,mean),apply(mse_P_SLM2,2,mean),
                                        apply(mse_P_CASDMM3,2,mean),apply(mse_P_SLM3,2,mean),
                                        apply(mse_P_CASDMM4,2,mean),apply(mse_P_SLM4,2,mean),
                                        apply(mse_P_CASDMM5,2,mean),apply(mse_P_SLM5,2,mean)),
                                  Q=(rep(c(rep("CASBAH",samples),rep("SLM",samples)),10)),
                                  cov=paste0("scenario ",rep(1:5,each=samples*2))))
mse_P_boxplot$cov=as.character(mse_P_boxplot$cov)
mse_P_boxplot$Q=as.character(mse_P_boxplot$Q)
mse_P_boxplot$Xi=as.numeric(mse_P_boxplot$Xi)

pdf(file="mse_P_sim.pdf",width=10, height=5)
ggplot(mse_P_boxplot, aes(x=cov, y=Xi, fill=Q)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=14),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("MSE Post-Treatment") +
  xlab("")
dev.off()


# --- Boxplots for Y: outcome ----

#cbPalette <- c("#007399")   #, "#00ace6", "#80dfff")

bias_Y=matrix(c(apply(bias_Y_CASDMM1,2,mean),apply(bias_Y_SLM1,2,mean),
                apply(bias_Y_CASDMM2,2,mean),apply(bias_Y_SLM2,2,mean),
                apply(bias_Y_CASDMM3,2,mean),apply(bias_Y_SLM3,2,mean),
                apply(bias_Y_CASDMM4,2,mean),apply(bias_Y_SLM4,2,mean),
                apply(bias_Y_CASDMM5,2,mean),apply(bias_Y_SLM5,2,mean, na.rm = FALSE)), ncol=10)

round(apply(bias_Y,2,median,na.rm = TRUE), 4)
round(apply(bias_Y,2,IQR,na.rm = TRUE), 4)

bias_Y_boxplot=as.data.frame(cbind(Xi=c(bias_Y),
                                   Q=(rep(c(rep("CASBAH",samples),rep("SLM",samples)),10)),
                                   cov=paste0("scenario ",rep(1:5,each=samples*2))))
bias_Y_boxplot$cov=as.character(bias_Y_boxplot$cov)
bias_Y_boxplot$Q=as.character(bias_Y_boxplot$Q)
bias_Y_boxplot$Xi=as.numeric(bias_Y_boxplot$Xi)

pdf(file="bias_Y_sim.pdf",width=10, height=5)
ggplot(bias_Y_boxplot, aes(x=cov, y=Xi, fill=Q)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  geom_hline(yintercept = 0, col="#D90224", size=0.4) +
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=14),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("Bias Outcome") +
  xlab("")
dev.off()

#cbPalette <- c("#80dfff")   #, "#00ace6", "#80dfff")

mse_Y_boxplot=as.data.frame(cbind(Xi=c(apply(mse_Y_CASDMM1,2,mean),apply(mse_Y_SLM1,2,mean),
                                       apply(mse_Y_CASDMM2,2,mean),apply(mse_Y_SLM2,2,mean),
                                       apply(mse_Y_CASDMM3,2,mean),apply(mse_Y_SLM3,2,mean),
                                       apply(mse_Y_CASDMM4,2,mean),apply(mse_Y_SLM4,2,mean),
                                       apply(mse_Y_CASDMM5,2,mean),apply(mse_Y_SLM5,2,mean)),
                                  Q=(rep(c(rep("CASBAH",samples),rep("SLM",samples)),10)),
                                  cov=paste0("scenario ",rep(1:5,each=samples*2))))
mse_Y_boxplot$cov=as.character(mse_Y_boxplot$cov)
mse_Y_boxplot$Q=as.character(mse_Y_boxplot$Q)
mse_Y_boxplot$Xi=as.numeric(mse_Y_boxplot$Xi)

pdf(file="mse_Y_sim.pdf",width=10, height=5)
ggplot(mse_Y_boxplot, aes(x=cov, y=Xi, fill=Q)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=14),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("MSE Outcome") +
  xlab("") +
  ylim(c(0,120))
dev.off()

# --- Boxplots for principal strata ----

Princ_strata_1=lapply(1:samples,function(i) order_strata(CASDMM_1[[i]]))
matrix_PS_1=matrix(unlist(Princ_strata_1), ncol=3, byrow=TRUE)
Princ_strata_2=lapply(1:samples,function(i) order_strata(CASDMM_2[[i]]))
matrix_PS_2=matrix(unlist(Princ_strata_2), ncol=3, byrow=TRUE)
Princ_strata_3=lapply(1:samples,function(i) order_strata(CASDMM_3[[i]]))
matrix_PS_3=matrix(unlist(Princ_strata_3), ncol=3, byrow=TRUE)
Princ_strata_4=lapply(1:samples,function(i) order_strata(CASDMM_4[[i]]))
matrix_PS_4=matrix(unlist(Princ_strata_4), ncol=3, byrow=TRUE)
Princ_strata_5=lapply(1:samples,function(i) order_strata(CASDMM_5[[i]]))
matrix_PS_5=matrix(unlist(Princ_strata_5), ncol=3, byrow=TRUE)

boxplots_PS(matrix_PS=matrix_PS_1,scenario_n=1)
boxplots_PS(matrix_PS=matrix_PS_2,scenario_n=2)
boxplots_PS(matrix_PS=matrix_PS_3,scenario_n=3)
boxplots_PS(matrix_PS=matrix_PS_4,scenario_n=4)
boxplots_PS(matrix_PS=matrix_PS_5,scenario_n=5)

boxplots_PS_est(matrix_PS=matrix_PS_1,scenario_n=1)
boxplots_PS_est(matrix_PS=matrix_PS_2,scenario_n=2)
boxplots_PS_est(matrix_PS=matrix_PS_3,scenario_n=3)
boxplots_PS_est(matrix_PS=matrix_PS_4,scenario_n=4)
boxplots_PS_est(matrix_PS=matrix_PS_5,scenario_n=5)
