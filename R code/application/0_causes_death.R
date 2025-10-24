#################################################################
# --- clean dataset for death 2000-2005
#################################################################

All_causes_death <- read.delim("All_causes_death_2000.2005.txt")
All_causes_death=All_causes_death[,c(2,3,5,4,7)]
All_causes_death=All_causes_death[!is.na(All_causes_death$County.Code),]
colnames(All_causes_death)[4]="all_causes"
colnames(All_causes_death)[5]="all_causes_ADJ"
All_causes_death$Population=as.integer(All_causes_death$Population)
All_causes_death$all_causes=as.integer(All_causes_death$all_causes)
All_causes_death$all_causes_ADJ=as.integer(All_causes_death$all_causes_ADJ)

Mental_and_behavioural_disorders <- read.delim("Mental_and_behavioural_disorders_2000.2005.txt")
Mental_and_behavioural_disorders=Mental_and_behavioural_disorders[!is.na(Mental_and_behavioural_disorders$County.Code),]
Mental_and_behavioural_disorders[Mental_and_behavioural_disorders$Deaths=="Suppressed",4]=NaN
All_causes_death$mental_dis=as.integer(Mental_and_behavioural_disorders[,4])
Mental_and_behavioural_disorders[Mental_and_behavioural_disorders[,6]=="Suppressed",6]=NaN
All_causes_death$mental_dis_ADJ=as.integer(Mental_and_behavioural_disorders[,6])

Diseases_nervous_system <- read.delim("Diseases_nervous_system_2000.2005.txt")
Diseases_nervous_system=Diseases_nervous_system[!is.na(Diseases_nervous_system$County.Code),]
Diseases_nervous_system[Diseases_nervous_system$Deaths=="Suppressed",4]=NaN
All_causes_death$nervous_dis=as.integer(Diseases_nervous_system[,4])
Diseases_nervous_system[Diseases_nervous_system[,6]=="Suppressed",6]=NaN
All_causes_death$nervous_dis_ADJ=as.integer(Diseases_nervous_system[,6])

Diseases_circulatory_system <- read.delim("Diseases_circulatory_system_2000.2005.txt")
Diseases_circulatory_system=Diseases_circulatory_system[!is.na(Diseases_circulatory_system$County.Code),]
Diseases_circulatory_system[Diseases_circulatory_system$Deaths=="Suppressed",4]=NaN
All_causes_death$circulatory_dis=as.integer(Diseases_circulatory_system[,4])
Diseases_circulatory_system[Diseases_circulatory_system[,6]=="Suppressed",6]=NaN
All_causes_death$circulatory_dis_ADJ=as.integer(Diseases_circulatory_system[,6])

Diseases_respiratory_system <- read.delim("Diseases_respiratory_system_2000.2005.txt")
Diseases_respiratory_system=Diseases_respiratory_system[!is.na(Diseases_respiratory_system$County.Code),]
Diseases_respiratory_system[Diseases_respiratory_system$Deaths=="Suppressed",4]=NaN
All_causes_death$respiratory_dis=as.integer(Diseases_respiratory_system[,4])
Diseases_respiratory_system[Diseases_respiratory_system[,6]=="Suppressed",6]=NaN
All_causes_death$respiratory_dis_ADJ=as.integer(Diseases_respiratory_system[,6])

Congenital_malformations <- read.delim2("Congenital_malformations_2000.2005.txt")
Congenital_malformations=Congenital_malformations[!is.na(Congenital_malformations$County.Code),]
Congenital_malformations[Congenital_malformations$Deaths=="Suppressed",4]=NaN
All_causes_death$malformations=as.integer(Congenital_malformations[,4])
Congenital_malformations[Congenital_malformations[,6]=="Suppressed",6]=NaN
All_causes_death$malformations_ADJ=as.integer(Congenital_malformations[,6])

Causes_death_2000.2005=All_causes_death

# summary dataset
summary(Causes_death_2000.2005)

par(mfrow=c(3,2))
for (i in c(4:15)){
  hist(Causes_death_2000.2005[,i], nclass=150, main=colnames(Causes_death_2000.2005)[i],
       xlab=" ", xlim=c(0, quantile(Causes_death_2000.2005[,i],0.99, na.rm=TRUE)))
}

save.image("Dataset_Causes_Death_2000.2005.RData")

#################################################################
# --- clean dataset for death 2010-2016
#################################################################

All_causes_death <- read.delim("All_causes_death_2010.2016.txt")
All_causes_death=All_causes_death[,c(2,3,5,4,7)]
All_causes_death=All_causes_death[!is.na(All_causes_death$County.Code),]
colnames(All_causes_death)[4]="all_causes"
colnames(All_causes_death)[5]="all_causes_ADJ"
All_causes_death$Population=as.integer(All_causes_death$Population)
All_causes_death$all_causes=as.integer(All_causes_death$all_causes)
All_causes_death$all_causes_ADJ=as.integer(All_causes_death$all_causes_ADJ)

Mental_and_behavioural_disorders <- read.delim("Mental_and_behavioural_disorders_2010.2016.txt")
Mental_and_behavioural_disorders=Mental_and_behavioural_disorders[!is.na(Mental_and_behavioural_disorders$County.Code),]
Mental_and_behavioural_disorders[Mental_and_behavioural_disorders$Deaths=="Suppressed",4]=NaN
All_causes_death$mental_dis=as.integer(Mental_and_behavioural_disorders[,4])
Mental_and_behavioural_disorders[Mental_and_behavioural_disorders[,6]=="Suppressed",6]=NaN
All_causes_death$mental_dis_ADJ=as.integer(Mental_and_behavioural_disorders[,6])

Diseases_nervous_system <- read.delim("Diseases_nervous_system_2010.2016.txt")
Diseases_nervous_system=Diseases_nervous_system[!is.na(Diseases_nervous_system$County.Code),]
Diseases_nervous_system[Diseases_nervous_system$Deaths=="Suppressed",4]=NaN
All_causes_death$nervous_dis=as.integer(Diseases_nervous_system[,4])
Diseases_nervous_system[Diseases_nervous_system[,6]=="Suppressed",6]=NaN
All_causes_death$nervous_dis_ADJ=as.integer(Diseases_nervous_system[,6])

Diseases_circulatory_system <- read.delim("Diseases_circulatory_system_2010.2016.txt")
Diseases_circulatory_system=Diseases_circulatory_system[!is.na(Diseases_circulatory_system$County.Code),]
Diseases_circulatory_system[Diseases_circulatory_system$Deaths=="Suppressed",4]=NaN
All_causes_death$circulatory_dis=as.integer(Diseases_circulatory_system[,4])
Diseases_circulatory_system[Diseases_circulatory_system[,6]=="Suppressed",6]=NaN
All_causes_death$circulatory_dis_ADJ=as.integer(Diseases_circulatory_system[,6])

Diseases_respiratory_system <- read.delim("Diseases_respiratory_system_2010.2016.txt")
Diseases_respiratory_system=Diseases_respiratory_system[!is.na(Diseases_respiratory_system$County.Code),]
Diseases_respiratory_system[Diseases_respiratory_system$Deaths=="Suppressed",4]=NaN
All_causes_death$respiratory_dis=as.integer(Diseases_respiratory_system[,4])
Diseases_respiratory_system[Diseases_respiratory_system[,6]=="Suppressed",6]=NaN
All_causes_death$respiratory_dis_ADJ=as.integer(Diseases_respiratory_system[,6])

Congenital_malformations <- read.delim2("Congenital_malformations_2010.2016.txt")
Congenital_malformations=Congenital_malformations[!is.na(Congenital_malformations$County.Code),]
Congenital_malformations[Congenital_malformations$Deaths=="Suppressed",4]=NaN
All_causes_death$malformations=as.integer(Congenital_malformations[,4])
Congenital_malformations[Congenital_malformations[,6]=="Suppressed",6]=NaN
All_causes_death$malformations_ADJ=as.integer(Congenital_malformations[,6])

Causes_death_2010.2016=All_causes_death

# summary dataset
summary(Causes_death_2010.2016)

par(mfrow=c(3,2))
for (i in c(4:15)){
  hist(Causes_death_2010.2016[,i], nclass=150, main=colnames(Causes_death_2010.2016)[i],
       xlab=" ", xlim=c(0, quantile(Causes_death_2010.2016[,i],0.99, na.rm=TRUE)))
}

save.image("Dataset_Causes_Death_2010.2016.RData")

