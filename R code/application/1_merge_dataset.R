
#library
library(data.table)
library(WeightIt) 
library(cobalt) 

#load datasets
load("Dataset_Policies.RData")
load("Dataset_Causes_Death_2000.2005.RData")
load("Dataset_Causes_Death_2010.2016.RData")

diff_deaths_rate<-cbind(Causes_death_2000.2005[,1:3],
                        Causes_death_2010.2016[,4:15]-Causes_death_2000.2005[,4:15])

par(mfrow=c(1,1))
hist(diff_deaths_rate[,4], nclass=100)
hist(diff_deaths_rate[,5], nclass=100)

summary(diff_deaths_rate[,4])
summary(diff_deaths_rate[,5])

# keep only all causes
All_causes=na.omit(diff_deaths_rate[,2:5])

# clean zigler's data
data$FIPS=as.integer(data$FIPS)

counties<-intersect(data$FIPS,All_causes$County.Code)
length(unique(data$FIPS))
counties_where<-sapply(counties, function(i) which(data$FIPS==i))

data_reduced<-data[1:2,]
for (i in 3:length(counties_where)){
  if (length(counties_where[[i]])==1){
    data_reduced=rbind(data_reduced,data[counties_where[[i]],])
  }else{
    summary_FIP<-data[c(counties_where[[i]]),]
    unique_fip<-c(unique(summary_FIP$FIPS),
                  unique(summary_FIP$a),
                  apply(summary_FIP[,3:6],2,mean),
                  apply(summary_FIP[,7:16]*summary_FIP$logpop/sum(summary_FIP$logpop),2,sum),
                  log(sum(exp(summary_FIP$logpop))),
                  sum(summary_FIP$smokerate2000*summary_FIP$logpop/sum(summary_FIP$logpop)),
                  apply(summary_FIP[,19:22],2,mean))
    unique_fip<-as.data.frame(t(unique_fip))
    colnames(unique_fip)<-colnames(data_reduced)
    data_reduced=rbind(data_reduced,unique_fip)
  }
}

FIPS_codes<-unique(data_reduced$FIPS)
counties_deaths<-sapply(FIPS_codes, function(i) which(All_causes$County.Code==i))

setdiff(counties,data$FIPS)

# --- final dataset: merge of death_causes + pollution policies ---

data_merged<-cbind(data_reduced,
                   All_causes[counties_deaths,])


plot(data_merged$pm_diff, data_merged$all_causes, pch=19)
plot(data_merged$pm_diff, data_merged$all_causes_ADJ, pch=19)

save(data_merged, file = "Dataset_merged.RData")
