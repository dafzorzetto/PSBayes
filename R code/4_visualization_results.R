#########################################################################
#     ---  MODEL COMPARISON   ----
#     ---  visualize results  ----
#     ---  tables and plots   ----
#########################################################################

# load data
load("estimands.RData")

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

cbPalette <- c("#007399")   #, "#00ace6", "#80dfff")

bias_P_boxplot=as.data.frame(cbind(Xi=c(apply(bias_P_CASDMM1,2,mean),
                                        apply(bias_P_CASDMM2,2,mean),
                                        apply(bias_P_CASDMM3,2,mean),
                                        apply(bias_P_CASDMM4,2,mean),
                                        apply(bias_P_CASDMM5,2,mean)),
                                 Q=(rep(rep("CASDMM",samples),5)),
                                 cov=paste0("scenario ",rep(1:5,each=samples))))
bias_P_boxplot$cov=as.character(bias_P_boxplot$cov)
bias_P_boxplot$Q=as.character(bias_P_boxplot$Q)
bias_P_boxplot$Xi=as.numeric(bias_P_boxplot$Xi)

#pdf(file="bias_P_sim.pdf",width=10, height=5)
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
#dev.off()

cbPalette <- c("#80dfff")   #, "#00ace6", "#80dfff")

mse_P_boxplot=as.data.frame(cbind(Xi=c(apply(mse_P_CASDMM1,2,mean),
                                        apply(mse_P_CASDMM2,2,mean),
                                        apply(mse_P_CASDMM3,2,mean),
                                        apply(mse_P_CASDMM4,2,mean),
                                        apply(mse_P_CASDMM5,2,mean)),
                                   Q=(rep(rep("CASDMM",samples),5)),
                                   cov=paste0("scenario ",rep(1:5,each=samples))))
mse_P_boxplot$cov=as.character(mse_P_boxplot$cov)
mse_P_boxplot$Q=as.character(mse_P_boxplot$Q)
mse_P_boxplot$Xi=as.numeric(mse_P_boxplot$Xi)

#pdf(file="mse_P_sim.pdf",width=10, height=5)
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
#dev.off()


# --- Boxplots for Y: outcome ----

cbPalette <- c("#007399")   #, "#00ace6", "#80dfff")

bias_Y_boxplot=as.data.frame(cbind(Xi=c(apply(bias_Y_CASDMM1,2,mean),
                                        apply(bias_Y_CASDMM2,2,mean),
                                        apply(bias_Y_CASDMM3,2,mean),
                                        apply(bias_Y_CASDMM4,2,mean),
                                        apply(bias_Y_CASDMM5,2,mean)),
                                   Q=(rep(rep("CASDMM",samples),5)),
                                   cov=paste0("scenario ",rep(1:5,each=samples))))
bias_Y_boxplot$cov=as.character(bias_Y_boxplot$cov)
bias_Y_boxplot$Q=as.character(bias_Y_boxplot$Q)
bias_Y_boxplot$Xi=as.numeric(bias_Y_boxplot$Xi)

#pdf(file="bias_Y_sim.pdf",width=10, height=5)
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
#dev.off()

cbPalette <- c("#80dfff")   #, "#00ace6", "#80dfff")

mse_Y_boxplot=as.data.frame(cbind(Xi=c(apply(mse_Y_CASDMM1,2,mean),
                                       apply(mse_Y_CASDMM2,2,mean),
                                       apply(mse_Y_CASDMM3,2,mean),
                                       apply(mse_Y_CASDMM4,2,mean),
                                       apply(mse_Y_CASDMM5,2,mean)),
                                  Q=(rep(rep("CASDMM",samples),5)),
                                  cov=paste0("scenario ",rep(1:5,each=samples))))
mse_Y_boxplot$cov=as.character(mse_Y_boxplot$cov)
mse_Y_boxplot$Q=as.character(mse_Y_boxplot$Q)
mse_Y_boxplot$Xi=as.numeric(mse_Y_boxplot$Xi)

#pdf(file="bias_Y_sim.pdf",width=10, height=5)
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
  xlab("")
#dev.off()
