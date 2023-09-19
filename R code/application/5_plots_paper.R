##############################################################
# ---     GRAPH FOR PAPER    ----
##############################################################

#load results
load("C:/Users/dafne/Desktop/DOTTORATO/PROGETTO PhD _2/data/results_model_10vs12.RData")

#libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)

##############################################################

results=Gibbs_CASDMM
results=Gibbs_CASDMM_cov

##############################################################

table(results$S_all_cluster)
table(results$S_strata_cluster)

par(mfrow=c(2,2))
hist(P_obs[T_var==0], nclass=40, xlim=c(-8,2), main="P(0) obs")
hist(results$post_P_0_imp[T_var==1], nclass=40, xlim=c(-8,2), main="P(0) imput")
hist(P_obs[T_var==1], nclass=40, xlim=c(-6,0), main="P(1) obs")
hist(results$post_P_1_imp[T_var==0], nclass=40, xlim=c(-6,0), main="P(1) imput")

par(mfrow=c(2,2))
hist(Y_obs[T_var==0], nclass=40, xlim=c(-200,50), main="Y(0) obs")
hist(results$post_Y_0_imp[T_var==1], nclass=40, xlim=c(-200,50), main="Y(0) imput")
hist(Y_obs[T_var==1], nclass=40, xlim=c(-250,50), main="Y(1) obs")
hist(results$post_Y_1_imp[T_var==0], nclass=40, xlim=c(-250,50), main="Y(1) imput")

par(mfrow=c(2,1))
hist(results$post_P_1_imp-results$post_P_0_imp,nclass=50, main="P(1)-P(0)")
hist(results$post_Y_1_imp-results$post_Y_0_imp,nclass=50, main="Y(1)-Y(0)")

par(mfrow=c(3,1))
hist((results$post_P_1_imp-results$post_P_0_imp)[results$S_strata_cluster==1],nclass=50, main="P(1)-P(0)|strata=1")
hist((results$post_P_1_imp-results$post_P_0_imp)[results$S_strata_cluster==2],nclass=50, main="P(1)-P(0)|strata=2")
hist((results$post_P_1_imp-results$post_P_0_imp)[results$S_strata_cluster==3],nclass=50, main="P(1)-P(0)|strata=3")

par(mfrow=c(3,1))
hist((results$post_P_1_imp-results$post_P_0_imp)[results$S_strata_cluster=="-1"],nclass=50, main="P(1)-P(0)|strata=-1")
hist((results$post_P_1_imp-results$post_P_0_imp)[results$S_strata_cluster=="0"],nclass=50, main="P(1)-P(0)|strata=0")
hist((results$post_P_1_imp-results$post_P_0_imp)[results$S_strata_cluster=="1"],nclass=50, main="P(1)-P(0)|strata=+1")

par(mfrow=c(3,1))
hist((results$post_Y_1_imp-results$post_Y_0_imp)[results$S_strata_cluster=="-1"],nclass=50, main="Y(1)-Y(0)|strata=-1")
hist((results$post_Y_1_imp-results$post_Y_0_imp)[results$S_strata_cluster=="0"],nclass=50, main="Y(1)-Y(0)|strata=0")
hist((results$post_Y_1_imp-results$post_Y_0_imp)[results$S_strata_cluster=="1"],nclass=50, main="Y(1)-Y(0)|strata=+1")


mean((results$post_Y_1_imp-results$post_Y_0_imp)[results$S_strata_cluster=="-1"])
mean((results$post_Y_1_imp-results$post_Y_0_imp)[results$S_strata_cluster=="0"])
mean((results$post_Y_1_imp-results$post_Y_0_imp)[results$S_strata_cluster=="1"])

covariates=sapply(c("-1","0","1"),function(c)
  apply(dataset_matched[which(results$S_strata_cluster==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                          "PctFemale","PctPoor")],2,mean))

##############################################################
# ---    boxplot ----
##############################################################


cbPalette <- c("#35D90B", "#F0D400", "#D90224")

data_frame=as.data.frame(cbind(P=results$post_P_1_imp-results$post_P_0_imp,
                               cl=results$S_strata_cluster))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=factor(data_frame$cl)

pdf(file="dati_P.pdf",width=9, height=6)
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
  ggtitle("P(1)-P(0)")
dev.off()

data_frame=as.data.frame(cbind(P=results$post_Y_1_imp-results$post_Y_0_imp,
                               cl=results$S_strata_cluster))
data_frame$P=as.numeric(data_frame$P)
data_frame$cl=factor(data_frame$cl)

pdf(file="dati_Y.pdf",width=9, height=6)
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
  ggtitle("Y(1)-Y(0)")
dev.off()

##############################################################
##############################################################
# ---    spiderplot ----
##############################################################

data_spider=as.data.frame( rbind(max=rep(1,6),min=rep(-0.15,6),
                                 t(covariates),
                                 mean=apply(dataset_matched[,c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                                                               "PctFemale","PctPoor")],2,mean)))
colnames(data_spider)=c("Black","Hispanic","Education","Urban","Female","Poor")
rownames(data_spider)=c("max","min","strata: -1","strata: 0","strata: +1","mean pop.")

# a spider plot for each group:
par(mfrow=c(1,1))
for(cl in 1:3){
  pdf(file=paste0("spiderplot_",cl,".pdf"))
  data_cl= data_spider[c(1,2,6,cl+2),]
  radarchart( data_cl  , axistype=1 , 
              pcol=c("grey",cbPalette[cl]), plwd=4 ,plty=1,
              cglcol="grey", cglty=1, axislabcol="black", caxislabels=rep("",5), cglwd=0.8, 
              vlcex=0.8 
  )
  legend(x=0.8, y=-0.8, legend=rownames(data_spider[c(cl+2,6),]), bty = "n", pch=20 , 
         col=c(cbPalette[cl],"grey") ,  cex=0.7, pt.cex=3)
  dev.off()
}