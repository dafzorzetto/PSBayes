

prob_diss<-function(x,a,L){

  alpha=a[1]+a[2]*x
  rho_1=(pnorm(alpha))
  rho_2=(pnorm(alpha)^2)
  
  mix_rho= rho_2-2*rho_1
    
  prob=rho_2*(1-(1+mix_rho)^L)/(-mix_rho)
  
  return(prob)
}

a=matrix(c(0,0,
           0,0.2,
           0,0.5,
           0,1,
           0,2,
           0,-0.5,
           0,-1,
           0.5,-0.5,
           1,-0.5), ncol=2, byrow=TRUE)

funct_1<-function(x){
  prob_diss(x,a=a[1,],100)
}
funct_2<-function(x){
  prob_diss(x,a=a[2,],100)
}
funct_3<-function(x){
  prob_diss(x,a=a[3,],100)
}
funct_4<-function(x){
  prob_diss(x,a=a[4,],100)
}
funct_5<-function(x){
  prob_diss(x,a=a[5,],100)
}
funct_6<-function(x){
  prob_diss(x,a=a[6,],100)
}
funct_7<-function(x){
  prob_diss(x,a=a[7,],100)
}
funct_8<-function(x){
  prob_diss(x,a=a[8,],100)
}
funct_9<-function(x){
  prob_diss(x,a=a[9,],100)
}

legend_label<-sapply(1:9, function(i) bquote(beta[0]==.(a[i,1]) ~ beta[1]==.(a[i,2])))
legend_color<-c("black","#0C00F0" ,"#0080F0" ,"#00D8F0" ,"#00F0B8","#E01E00" ,"#E08F00","#E000B7","#E094D2" )



curve(funct_1,
      from=-5,to=7.5, 
      ylim=c(-0.3,1.3), xlim=c(-5,11.5),
      ylab=' ',
      col=legend_color[1],
      frame.plot = FALSE, yaxt = "n", xaxt = "n")
title(ylab=expression(paste("Pr[",S(0)==S(1),"]")), 
      mgp=c(2,0,1))
axis(2, at = seq(0, 1, by = 0.25))
axis(1, at = seq(-5, 7.5, by = 2.5))
curve(funct_2, add=TRUE,from=-5,to=7.5,  col=legend_color[2])
curve(funct_3, add=TRUE,from=-5,to=7.5,  col=legend_color[3])
curve(funct_4, add=TRUE,from=-5,to=7.5,  col=legend_color[4])
curve(funct_5, add=TRUE,from=-5,to=7.5,  col=legend_color[5])
curve(funct_6, add=TRUE,from=-5,to=7.5,  col=legend_color[6])
curve(funct_7, add=TRUE,from=-5,to=7.5,  col=legend_color[7])
curve(funct_8, add=TRUE,from=-5,to=7.5,  col=legend_color[8])
curve(funct_9, add=TRUE,from=-5,to=7.5,  col=legend_color[9])
legend("right", legend=legend_label, 
       col=legend_color, lty=1)
