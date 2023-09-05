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

table_10_str_b=table(model_data_REGweights_10$S_strata_cluster[,2])

##############################################################
# ---       imputed PM and Mortality   ----
##############################################################

P_0_10=apply(model_data_REGweights_10$post_P_0,1,median)
P_1_10=apply(model_data_REGweights_10$post_P_1,1,median)

Y_0_10=apply(model_data_REGweights_10$post_Y_0_imp,1,median)
Y_1_10=apply(model_data_REGweights_10$post_Y_1_imp,1,median)

##############################################################
# ---       groups distributions     ----
##############################################################

P_10_str_b=lapply(names(table_10_str_b[table_10_str_b>2]), function(l)
  P_1_10[model_data_REGweights_10$S_strata_cluster[,2]==l]-
    P_0_10[model_data_REGweights_10$S_strata_cluster[,2]==l])

Y_10_str_b=lapply(names(table_10_str_b[table_10_str_b>2]), function(l)
  Y_1_10[model_data_REGweights_10$S_strata_cluster[,2]==l]-
    Y_0_10[model_data_REGweights_10$S_strata_cluster[,2]==l])

seq_10_str_b=lapply(names(table_10_str_b[table_10_str_b>2]), function(l) 
  rep(l,sum(model_data_REGweights_10$S_strata_cluster[,2]==l)))

##############################################################
# ---     covariates analysis    ----
##############################################################

cov_10_str_mode=sapply(names(table_10_str_mode),function(c)
  apply(dataset_matched[which(model_data_REGweights_10$S_strata_cluster[,3]==c),
                        c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                          "PctFemale","PctPoor")],2,mean))

##############################################################
##############################################################


##############################################################
##############################################################
# ---    boxplot ----
##############################################################


cbPalette <- c("#35D90B", "#F0D400", "#D90224")

boxplot_P<-function(values, seq_cluster){
  
  data_frame=as.data.frame(cbind(P=unlist(values),cl=unlist(seq_cluster)))
  data_frame$P=as.numeric(data_frame$P)
  data_frame$cl[data_frame$cl==2]="-1"
  data_frame$cl[data_frame$cl==1]="0"
  data_frame$cl[data_frame$cl==3]="+1"
  data_frame$cl=factor(data_frame$cl, levels=c("-1", "0", "+1"))
  
  setwd("C:/Users/dafne/Desktop/plot temporanei/princ_strata/CA")
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
    coord_cartesian(xlim = c(-2, 0.5)) +
    ggtitle("P(1)-P(0)")
  dev.off()
  
}

boxplot_P(values=P_10_str_b, seq_cluster=seq_10_str_b)

boxplot_Y<-function(values, seq_cluster){
  
  data_frame=as.data.frame(cbind(P=unlist(values),cl=unlist(seq_cluster)))
  data_frame$P=as.numeric(data_frame$P)
  data_frame$cl[data_frame$cl==2]="-1"
  data_frame$cl[data_frame$cl==1]="0"
  data_frame$cl[data_frame$cl==3]="+1"
  data_frame$cl=factor(data_frame$cl, levels=c("-1", "0", "+1"))
  
  setwd("C:/Users/dafne/Desktop/plot temporanei/princ_strata/CA")
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
  
}

boxplot_Y(values=Y_10_str_b, seq_cluster=seq_10_str_b)


##############################################################
##############################################################
# ---    spiderplot ----
##############################################################

data_spider=as.data.frame( rbind(max=rep(1,6),min=rep(-0.15,6),
                  t(cov_10_str_mode),
                 mean=apply(dataset_matched[,c("PctBlack","PctHisp","PctHighSchool","PctUrban",
                                          "PctFemale","PctPoor")],2,mean)))
colnames(data_spider)=c("Black","Hispanic","Education","Urban","Female","Poor")
rownames(data_spider)=c("max","min","strata: -1","strata: 0","strata: +1","mean pop.")

# a spiderplot for each group:

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

