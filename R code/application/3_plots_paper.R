##############################################################
# ---     GRAPHS FOR PAPER    ----
##############################################################

#load results
load("results_application.RData")

#libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(fmsb)

##############################################################

covariates=sapply(c("-1","0","1"),function(c)
  apply(matrix_COV[which(results$S_strata_cluster==c),],2,mean))


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
data_frame$P[which(data_frame$P<(-100))]<-0
data_frame$P[which(data_frame$P>(100))]<-0

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
# ---    spiderplot ----
##############################################################

data_spider=as.data.frame( rbind(max=rep(1,dim(matrix_COV)[2]),min=rep(-1,dim(matrix_COV)[2]),
                                 t(covariates),
                                 mean=apply(matrix_COV,2,mean)))
#colnames(data_spider)=c("Black","Hispanic","Education","Urban","Female","Poor")
colnames(data_spider)=c("Urban","Black","Hispanic","Education","Income","Poor",
                        "Female","Occupied","Moved", "House Value", "Population", "Smoke Rate",
                        "Dew Point", "Temperature", "Humidity")
rownames(data_spider)=c("max","min","strata: -1","strata: 0","strata: +1","mean pop.")

# a spider plot for each group:
par(mfrow=c(1,1))
for(cl in 1:3){
  pdf(file=paste0("spiderplot_",cl,"_new.pdf"))
  data_cl= data_spider[c(1,2,6,cl+2),]
  radarchart( data_cl  , axistype=1 , 
              pcol=c("grey",cbPalette[cl]), plwd=4 ,plty=1,
              cglcol="grey", cglty=1, axislabcol="black", caxislabels=rep("",5), cglwd=0.8, 
              vlcex=0.8 
  )
  legend(x=0.97, y=-1.05, legend=rownames(data_spider[c(cl+2,6),]), bty = "n", pch=20 , 
         col=c(cbPalette[cl],"grey") ,  cex=0.7, pt.cex=3)
  dev.off()
}


#################################################################
CE_Y <- results$all_Y1-results$all_Y0

PCE_Y_chains2 <- t(sapply(1:5000, function(r)
  sapply(c(-1,0,1), function(v) 
    median(CE_Y[which(results$chians_strata[r,]==v),r]))))

CE_median = apply(PCE_Y_chains2,2,median)
CE_mean =apply(PCE_Y_chains2,2,mean)
CE_quanitles =apply(PCE_Y_chains2,2,quantile, prob=c(0.025,0.05,0.95,0.975))

par(mfrow=c(3,1))
plot(PCE_Y_chains2[,1], type = "l") #, ylim=c(-lim_val,lim_val))
abline(h=c(CE_quanitles[2,1],CE_quanitles[3,1]), col=2)
abline(h=0, col=3)
plot(PCE_Y_chains2[,2], type = "l") #, ylim=c(-lim_val,lim_val))
abline(h=c(CE_quanitles[2,2],CE_quanitles[3,2]), col=2)
abline(h=0, col=3)
plot(PCE_Y_chains2[,3], type = "l") #, ylim=c(-lim_val,lim_val))
abline(h=c(CE_quanitles[2,3],CE_quanitles[3,3]), col=2)
abline(h=0, col=3)


cbPalette <- c("#35D90B", "#F0D400", "#D90224")

par(mfrow=c(1,1))
plot(0,0, type = "n", frame.plot = FALSE, axes = FALSE,
     ylim = c(-1, 1), xlim = c(-45, 40), xlab = " ", ylab = " ",
     main = expression("E[" ~ Y[i] ~ "(1)" ~  - Y[i] ~ "(0) | " ~ i %in% ~ "stratum ]"))
points(CE_mean, c(-1,0,1), 
       col=cbPalette, pch=16, cex=2)
segments(CE_quanitles[2,1],-1,CE_quanitles[3,1],-1, col=cbPalette[1], lwd = 3)
segments(CE_quanitles[2,2],0,CE_quanitles[3,2],0, col=cbPalette[2], lwd = 3)
segments(CE_quanitles[2,3],1,CE_quanitles[3,3],1, col=cbPalette[3], lwd = 3)
abline(v=0, col="#00ace6")
axis(1, at = seq(-30, 40, by = 10))
text(-40, -1.1, labels = "associative \n negative", pos = 3, cex = 0.8, col = "black")
text(-40, 0, labels = "dissociative", pos = 3, cex = 0.8, col = "black")
text(-40, 0.8, labels = "associative \n positive", pos = 3, cex = 0.8, col = "black")

ce_P <- results$post_P_1_imp - results$post_P_0_imp

CE_P_median = sapply(c("-1","0","1"),function(c) 
  median(ce_P[which(results$S_strata_cluster==c)]))
CE_P_mean = sapply(c("-1","0","1"),function(c) 
  mean(ce_P[which(results$S_strata_cluster==c)]))
CE_P_quanitles = sapply(c("-1","0","1"),function(c) 
  quantile(ce_P[which(results$S_strata_cluster==c)], prob=c(0.025,0.05,0.95,0.975)))

augmented_ce_P <- sapply(1:5000, function(c)
  sapply(c("-1","0","1"), function(v) 
    median(ce_P[which(results$chians_strata[c,]==v)])))

par(mfrow=c(3,1))
plot(augmented_ce_P[1,], type = "l") #, ylim=c(-lim_val,lim_val))
abline(h=0, col=2)
plot(augmented_ce_P[2,], type = "l") #, ylim=c(-lim_val,lim_val))
abline(h=0, col=2)
plot(augmented_ce_P[3,], type = "l") #, ylim=c(-lim_val,lim_val))
abline(h=0, col=2)

CE_Pa_median = apply(augmented_ce_P[,1:5000],1,median)
CE_Pa_mean =apply(augmented_ce_P[,1:5000],1,mean)
CE_Pa_quanitles =apply(augmented_ce_P[,1:5000],1,quantile, prob=c(0.025,0.05,0.95,0.975,0.25,0.75))

par(mfrow=c(1,1))
plot(0,0, type = "n", frame.plot = FALSE, axes = FALSE,
     ylim = c(-1.2, 1.2), xlim = c(-2, 1.2), xlab = " ", ylab = " ",
     main = expression("E[" ~ P[i] ~ "(1)" ~  - P[i] ~ "(0) | " ~ i %in% ~ "stratum ]"))
points(CE_Pa_mean, c(-1,0,1), 
       col=cbPalette, pch=16, cex=2)
segments(CE_Pa_quanitles[2,1],-1,CE_Pa_quanitles[3,1],-1, col=cbPalette[1], lwd = 3)
segments(CE_Pa_quanitles[2,2],0,CE_Pa_quanitles[3,2],0, col=cbPalette[2], lwd = 3)
segments(CE_Pa_quanitles[2,3],1,CE_Pa_quanitles[3,3],1, col=cbPalette[3], lwd = 3)
abline(v=0, col="#00ace6")
axis(1, at = round(seq(-1.4, 1, by = 0.2),1))
text(-1.85, -1.1, labels = "associative \n negative", pos = 3, cex = 0.8, col = "black")
text(-1.85, 0, labels = "dissociative", pos = 3, cex = 0.8, col = "black")
text(-1.85, 0.8, labels = "associative \n positive", pos = 3, cex = 0.8, col = "black")


PCE_P_mean <- sapply(c("-1","0","1"),function(c) 
  mean(ce_P[which(results$S_strata_cluster==c)]))
PCE_P_quantiles <- sapply(c("-1","0","1"),function(c) 
  quantile(ce_P[which(results$S_strata_cluster==c)], prob=c(0.025,0.05,0.95,0.975,0.25,0.75,0.5)))

par(mfrow=c(1,1))
plot(0,0, type = "n", frame.plot = FALSE, axes = FALSE,
     ylim = c(-1, 1), xlim = c(-1.5, 1), xlab = " ", ylab = " ")
points(PCE_P_mean, c(-1,0,1), 
       col=cbPalette, pch=16, cex=2)
segments(PCE_P_quantiles[2,1],-1,PCE_P_quantiles[3,1],-1, col=cbPalette[1], lwd = 3)
segments(PCE_P_quantiles[2,2],0,PCE_P_quantiles[3,2],0, col=cbPalette[2], lwd = 3)
segments(PCE_P_quantiles[2,3],1,PCE_P_quantiles[3,3],1, col=cbPalette[3], lwd = 3)
abline(v=0, col="#00ace6")
axis(1, at = seq(-1, 1, by = 0.2))
text(-1.35, -1.1, labels = "associative \n negative", pos = 3, cex = 0.8, col = "black")
text(-1.35, 0, labels = "dissociative", pos = 3, cex = 0.8, col = "black")
text(-1.35, 0.8, labels = "associative \n positive", pos = 3, cex = 0.8, col = "black")

par(mfrow=c(1,1))
plot(0,0, type = "n", frame.plot = FALSE, axes = FALSE,
     ylim = c(-1, 1), xlim = c(-2.5, 1.5), xlab = " ", ylab = " ")
points(PCE_P_mean, c(-1,0,1), 
       col=cbPalette, pch=16, cex=2)
segments(PCE_P_quantiles[5,1],-1,PCE_P_quantiles[6,1],-1, col=cbPalette[1], lwd = 3)
segments(PCE_P_quantiles[5,2],0,PCE_P_quantiles[6,2],0, col=cbPalette[2], lwd = 3)
segments(PCE_P_quantiles[5,3],1,PCE_P_quantiles[6,3],1, col=cbPalette[3], lwd = 3)
abline(v=0, col="#00ace6")
axis(1, at = seq(-1.5, 1.5, by = 0.5))
text(-2.35, -1.1, labels = "associative \n negative", pos = 3, cex = 0.8, col = "black")
text(-2.35, 0, labels = "dissociative", pos = 3, cex = 0.8, col = "black")
text(-2.35, 0.8, labels = "associative \n positive", pos = 3, cex = 0.8, col = "black")
