##############################################################
# ---     GRAPH FOR PAPER    ----
##############################################################

#load results
load("results_application.RData")

#libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)

##############################################################

##############################################################

table(results$S_all_cluster)
table(results$S_strata_cluster)

par(mfrow=c(2,2))
hist(P_obs[T_var==0], nclass=40, xlim=c(-8,2), main="P(0) obs")
hist(results$post_P_0_imp[T_var==1], nclass=40, xlim=c(-8,2), main="P(0) imput")
hist(P_obs[T_var==1], nclass=40, xlim=c(-8,2), main="P(1) obs")
hist(results$post_P_1_imp[T_var==0], nclass=40, xlim=c(-8,2), main="P(1) imput")

par(mfrow=c(2,2))
hist(Y_obs[T_var==0], nclass=40, xlim=c(-250,50), main="Y(0) obs")
hist(results$post_Y_0_imp[T_var==1], nclass=40, xlim=c(-250,50), main="Y(0) imput")
hist(Y_obs[T_var==1], nclass=40, xlim=c(-250,50), main="Y(1) obs")
hist(results$post_Y_1_imp[T_var==0], nclass=40, xlim=c(-250,50), main="Y(1) imput")

par(mfrow=c(2,1))
hist(results$post_P_1_imp-results$post_P_0_imp,nclass=50, main="P(1)-P(0)")
hist(results$post_Y_1_imp-results$post_Y_0_imp,nclass=50, main="Y(1)-Y(0)")

par(mfrow=c(3,1))
diff_P=(results$post_P_1_imp-results$post_P_0_imp)
hist(diff_P[results$S_strata_cluster=="-1"],nclass=50, main="P(1)-P(0)|strata=-1", xlim=c(min(diff_P),max(diff_P)))
hist(diff_P[results$S_strata_cluster=="0"],nclass=50, main="P(1)-P(0)|strata=0", xlim=c(min(diff_P),max(diff_P)))
hist(diff_P[results$S_strata_cluster=="1"],nclass=50, main="P(1)-P(0)|strata=+1", xlim=c(min(diff_P),max(diff_P)))

par(mfrow=c(3,1))
diff_Y=(results$post_Y_1_imp-results$post_Y_0_imp)
hist(diff_Y[results$S_strata_cluster=="-1"],nclass=50, main="Y(1)-Y(0)|strata=-1", xlim=c(min(diff_Y),max(diff_Y)))
hist(diff_Y[results$S_strata_cluster=="0"],nclass=50, main="Y(1)-Y(0)|strata=0", xlim=c(min(diff_Y),max(diff_Y)))
hist(diff_Y[results$S_strata_cluster=="1"],nclass=50, main="Y(1)-Y(0)|strata=+1", xlim=c(min(diff_Y),max(diff_Y)))

covariates=sapply(c("-1","0","1"),function(c)
  apply(matrix_X[which(results$S_strata_cluster==c),],2,mean))

##############################################################
# chains

par(mfrow=c(2,3))
for (i in 1:12){
  plot(results$post_eta[i,], type="l")
  abline(h=mean(results$post_eta[i,]), col="red")
}

par(mfrow=c(3,4))
for (i in 1:36){
  plot(results$chains_theta[i,], type="l")
  abline(h=mean(results$chains_theta[i,]), col="red")
}

##############################################################
# ---    boxplot ----
##############################################################

library(latex2exp)

cbPalette <- c("#35D90B", "#F0D400", "#D90224")

data_frame=as.data.frame(cbind(P=results$post_P_1_imp-results$post_P_0_imp,
                               cl=results$S_strata_cluster))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=factor(data_frame$cl)

#quantiles
quantile(data_frame$P[data_frame$cl=="-1"], prob=c(0.05,0.25,0.5,0.75,0.95))
quantile(data_frame$P[data_frame$cl=="0"], prob=c(0.05,0.25,0.5,0.75,0.95))
quantile(data_frame$P[data_frame$cl=="1"], prob=c(0.05,0.25,0.5,0.75,0.95))

IQR(data_frame$P[data_frame$cl=="-1"])
IQR(data_frame$P[data_frame$cl=="0"])
IQR(data_frame$P[data_frame$cl=="1"])

#pdf(file="dati_P.pdf",width=9, height=6)
ggplot(data_frame, aes(x=P, y=cl, fill=cl)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot()+
  geom_vline(xintercept = 0, col="#00ace6", size=1) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=16),
        #legend.text=element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        title =element_text(size=20),
        #legend.background = element_rect(fill='transparent'),
        #panel.grid.major = element_line(color = "grey",size = 0.45)
        axis.text.y = element_text(size=14),
  )+
  #geom_hline(yintercept = vero, color = "#0BC2CF", size=0.65)+
  ylab("") +
  xlab(" ") +
  ggtitle(expression(paste("E [ ", P[i](1)-P[i](0), " | stratum ]")))
#dev.off()

data_frame=as.data.frame(cbind(P=results$post_Y_1_imp-results$post_Y_0_imp,
                               cl=results$S_strata_cluster))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=factor(data_frame$cl)

#quantiles
quantile(data_frame$P[data_frame$cl=="-1"], prob=c(0.05,0.25,0.5,0.75,0.95))
quantile(data_frame$P[data_frame$cl=="0"], prob=c(0.05,0.25,0.5,0.75,0.95))
quantile(data_frame$P[data_frame$cl=="1"], prob=c(0.05,0.25,0.5,0.75,0.95))

IQR(data_frame$P[data_frame$cl=="-1"])
IQR(data_frame$P[data_frame$cl=="0"])
IQR(data_frame$P[data_frame$cl=="1"])


#pdf(file="dati_Y.pdf",width=9, height=6)
ggplot(data_frame, aes(x=P, y=cl, fill=cl)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot()+
  geom_vline(xintercept = 0, col="#00ace6", size=1) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=16),
        #legend.text=element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        title =element_text(size=20),
        #legend.background = element_rect(fill='transparent'),
        #panel.grid.major = element_line(color = "grey",size = 0.45)
        axis.text.y = element_text(size=14),
  )+
  #geom_hline(yintercept = vero, color = "#0BC2CF", size=0.65)+
  ylab("") +
  xlab(" ") +
  #coord_cartesian(xlim = c(-2, 0.5)) +
  ggtitle(expression(paste("E [ ", Y[i](1)-Y[i](0), " | stratum ]")))
#dev.off()

##############################################################
##############################################################
# ---    spiderplot ----
##############################################################

data_spider=as.data.frame( rbind(max=rep(1,dim(matrix_X)[2]),min=rep(-1,dim(matrix_X)[2]),
                                 t(covariates),
                                 mean=apply(matrix_X,2,mean)))
#colnames(data_spider)=c("Black","Hispanic","Education","Urban","Female","Poor")
colnames(data_spider)=c("Urban","Black","Hispanic","Education","Income","Poor",
                        "Female","Occupied","Moved", "House Value", "Population", "Smoke Rate",
                        "Dew Point", "Temperature", "Humidity")
rownames(data_spider)=c("max","min","strata: -1","strata: 0","strata: +1","mean pop.")

# a spider plot for each group:
par(mfrow=c(1,1))
for(cl in 1:3){
  #pdf(file=paste0("spiderplot_",cl,".pdf"))
  data_cl= data_spider[c(1,2,6,cl+2),]
  radarchart( data_cl  , axistype=1 , 
              pcol=c("grey",cbPalette[cl]), plwd=4 ,plty=1,
              cglcol="grey", cglty=1, axislabcol="black", caxislabels=rep("",5), cglwd=0.8, 
              vlcex=0.8 
  )
  legend(x=0.97, y=-1.05, legend=rownames(data_spider[c(cl+2,6),]), bty = "n", pch=20 , 
         col=c(cbPalette[cl],"grey") ,  cex=0.7, pt.cex=3)
  #dev.off()
}
