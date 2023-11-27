###################################################################################
# ---     results - plots: principal strata and causal effects     ---
###################################################################################

#create the dataset for the boxplots
order_strata<-function(sample_est){
  
  strata=sample_est$S_strata_cluster
  P_diff=sample_est$post_P_1_imp-sample_est$post_P_0_imp
  Y_diff=sample_est$post_Y_1_imp-sample_est$post_Y_0_imp
  
  n_strata=unique(strata)
  Princ_strata=sapply(n_strata, function(s) mean(P_diff[strata==s]))
  order_strata=order(Princ_strata, decreasing = FALSE)
  Princ_strata=Princ_strata[order_strata]
  Est_CE=sapply(n_strata[order_strata], function(s) mean(Y_diff[strata==s]))
  
  return(t(matrix(c(Princ_strata,Est_CE,n_strata), ncol=3)))
}

#libraries
library(ggplot2)

cbPalette <- c("#35D90B", "#F0D400", "#D90224")

boxplots_PS<-function(matrix_PS,scenario_n){
  
  delate= which(table(matrix_PS[,3])<(samples/2))
  if (length(delate)>0){
    matrix_PS=matrix_PS[-which(matrix_PS[,3]==delate),]
  }
  
  
  boxplot_df=as.data.frame(cbind(Xi=matrix_PS[,1],
                                 Q=matrix_PS[,3]))
  boxplot_df$Q=as.character(boxplot_df$Q)
  boxplot_df$Q[boxplot_df$Q==1]="V=-1"
  boxplot_df$Q[boxplot_df$Q==2]="V=0"
  boxplot_df$Q[boxplot_df$Q==3]="V=1"
  
  g<-ggplot(boxplot_df, aes(x=Q, y=Xi, fill=Q)) + 
    scale_fill_manual(values=cbPalette, name="strata")+
    geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
    #geom_hline(yintercept = sim_par$eta_1-sim_par$eta_0, col="#00ace6", size=0.5, linetype = "dotdash") +
    theme(panel.background = element_rect(fill='white'),
          plot.background = element_rect(fill ="white"),
          axis.title = element_text(size=14),
          legend.text=element_text(size=10),
          plot.title = element_text(hjust = 0.05),
          title =element_text(size=12),
          legend.background = element_rect(fill='transparent'),
          panel.grid.major = element_line(color = "grey",size = 0.1))+
    ylab("GATEs") +
    xlab("")+
    ggtitle(paste0("    Scenario ",scenario_n)) 
  
  pdf(file=paste0("PS_strata_",scenario_n,".pdf"),width=7, height=5)
  print(g)
  dev.off()
  
}

boxplots_PS_est<-function(matrix_PS,scenario_n){
  
  delate= which(table(matrix_PS[,3])<(samples/2))
  if (length(delate)>0){
    matrix_PS=matrix_PS[-which(matrix_PS[,3]==delate),]
  }
  
  
  boxplot_df=as.data.frame(cbind(Xi=matrix_PS[,2],
                                 Q=matrix_PS[,3]))
  boxplot_df$Q=as.character(boxplot_df$Q)
  boxplot_df$Q[boxplot_df$Q==1]="V=-1"
  boxplot_df$Q[boxplot_df$Q==2]="V=0"
  boxplot_df$Q[boxplot_df$Q==3]="V=1"
  
  g<-ggplot(boxplot_df, aes(x=Q, y=Xi, fill=Q)) + 
    scale_fill_manual(values=cbPalette, name="Princ.CE")+
    geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
    #geom_hline(yintercept = sim_par$eta_1-sim_par$eta_0, col="#00ace6", size=0.5, linetype = "dotdash") +
    theme(panel.background = element_rect(fill='white'),
          plot.background = element_rect(fill ="white"),
          axis.title = element_text(size=14),
          legend.text=element_text(size=10),
          plot.title = element_text(hjust = 0.05),
          title =element_text(size=12),
          legend.background = element_rect(fill='transparent'),
          panel.grid.major = element_line(color = "grey",size = 0.1))+
    ylab("GATEs") +
    xlab("")+
    ggtitle(paste0("    Scenario ",scenario_n)) 
  
  pdf(file=paste0("PS_CE_",scenario_n,".pdf"),width=7, height=5)
  print(g)
  dev.off()
  
}
