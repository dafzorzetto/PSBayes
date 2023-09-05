
#library
library(data.table)
library(WeightIt) 
library(cobalt) 

load("C:/Users/dafne/Downloads/causes_death_counties_/Dataset_Policies.RData")
load("C:/Users/dafne/Dropbox/DafneZorzetto/2_BNPCausal/2_PrincipalStratification/Data/Dataset_Causes_Death.RData")


head(data)
head(Dataset_all_causes)

Dataset_all_causes=Dataset_all_causes[!is.na(Dataset_all_causes$County.Code),]

head(data$FIPS)
tail(data$FIPS)
data$FIPS=as.integer(data$FIPS)

summary(Dataset_all_causes)

head(Dataset_all_causes$County.Code)
tail(Dataset_all_causes$County.Code)

counties<-intersect(data$FIPS,Dataset_all_causes$County.Code)
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
counties_deaths<-sapply(FIPS_codes, function(i) which(Dataset_all_causes$County.Code==i))


# --- final dataset: merge of death_causes + pollution policies ---

data_merged<-cbind(data_reduced,
                   Dataset_all_causes[counties_deaths,3:9])
