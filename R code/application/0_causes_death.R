
# upload all death causes by county level
# all the datasets has the same counties 

All_causes_death <- read.delim("All_causes_death.txt")
All_causes_death=All_causes_death[,-c(1,6)]
colnames(All_causes_death)[3]="all_causes"

Mental_and_behavioural_disorders <- read.delim("Mental_and_behavioural_disorders.txt")
All_causes_death$mental_dis=Mental_and_behavioural_disorders[,4]
All_causes_death$mental_dis[All_causes_death$mental_dis=="Suppressed"]=NaN
All_causes_death$mental_dis=as.integer(All_causes_death$mental_dis)

Diseases_nervous_system <- read.delim("Diseases_nervous_system.txt")
All_causes_death$nervous_dis=Diseases_nervous_system[,4]
All_causes_death$nervous_dis[All_causes_death$nervous_dis=="Suppressed"]=NaN
All_causes_death$nervous_dis=as.integer(All_causes_death$nervous_dis)

Diseases_circulatory_system <- read.delim("Diseases_circulatory_system.txt")
All_causes_death$circulatory_dis=Diseases_circulatory_system[,4]
All_causes_death$circulatory_dis[All_causes_death$circulatory_dis=="Suppressed"]=NaN
All_causes_death$circulatory_dis=as.integer(All_causes_death$circulatory_dis)

Diseases_respiratory_system <- read.delim("Diseases_respiratory_system.txt")
All_causes_death$respiratory_dis=Diseases_respiratory_system[,4]
All_causes_death$respiratory_dis[All_causes_death$respiratory_dis=="Suppressed"]=NaN
All_causes_death$respiratory_dis=as.integer(All_causes_death$respiratory_dis)

Congenital_malformations <- read.delim2("Congenital_malformations.txt")
All_causes_death$malformations=Congenital_malformations[,4]
All_causes_death$malformations[All_causes_death$malformations=="Suppressed"]=NaN
All_causes_death$malformations=as.integer(All_causes_death$malformations)

Dataset_all_causes=All_causes_death[, c(1,2,4,3,5:9)]

# summary dataset
summary(Dataset_all_causes)

par(mfrow=c(3,2))
for (i in c(4:9)){
  hist(Dataset_all_causes[,i], nclass=150, main=colnames(Dataset_all_causes)[i],
       xlab=" ", xlim=c(0, quantile(Dataset_all_causes[,i],0.99, na.rm=TRUE)))
}

save.image("Dataset_Causes_Death.RData")
