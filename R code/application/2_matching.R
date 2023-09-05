##############################################################
# ---       matching dataset     ----
##############################################################

#library
library(data.table)
library(WeightIt) 
library(cobalt) 
library(optmatch)
library(MatchIt)
library(tidyverse)

# load dataset:
load("C:/Users/dafne/Dropbox/DafneZorzetto/2_BNPCausal/2_PrincipalStratification/Data/final_dataset.RData")

##############################################################
# ---       preparing data     ----
##############################################################

# prepare data for matching
data_temp <- data_merged %>% 
  select(a,                                          #treatment
         pmbase2002_2004, pmfu, pm_diff,             # post-treatment
         all_causes,                                 #outcome
         PctUrban, PctBlack, PctHisp, PctHighSchool,
         MedianHHInc, PctPoor, PctFemale,
         PctOccupied, PctMovedIn5, MedianHValue,
         smokerate2000, avgdewpt, avgtemp, avgrelhum, Population,
         FIPS)
data_temp_rnames <- row.names(data_temp)

data_temp$all_causes=data_temp$all_causes/data_temp$Population
data_temp$PctUrban=1*(data_temp$PctUrban>median(data_temp$PctUrban))
data_temp$PctBlack=1*(data_temp$PctBlack>median(data_temp$PctBlack))
data_temp$PctHisp=1*(data_temp$PctHisp>median(data_temp$PctHisp))
data_temp$PctHighSchool=1*(data_temp$PctHighSchool>median(data_temp$PctHighSchool))
data_temp$MedianHHInc=1*(data_temp$MedianHHInc>median(data_temp$MedianHHInc))
data_temp$PctPoor=1*(data_temp$PctPoor>median(data_temp$PctPoor))
data_temp$PctFemale=1*(data_temp$PctFemale>median(data_temp$PctFemale))
data_temp$PctOccupied=1*(data_temp$PctOccupied>median(data_temp$PctOccupied))
data_temp$PctMovedIn5=1*(data_temp$PctMovedIn5>median(data_temp$PctMovedIn5))
data_temp$MedianHValue=1*(data_temp$MedianHValue>median(data_temp$MedianHValue))
data_temp$smokerate2000=1*(data_temp$smokerate2000>median(data_temp$smokerate2000))
data_temp$avgdewpt=1*(data_temp$avgdewpt>median(data_temp$avgdewpt))
data_temp$avgtemp=1*(data_temp$avgtemp>median(data_temp$avgtemp))
data_temp$avgrelhum=1*(data_temp$avgrelhum>median(data_temp$avgrelhum))
data_temp$Population=1*(data_temp$Population>median(data_temp$Population))


##############################################################
# ---       matching      ----
##############################################################


m_temp <- matchit(a ~ PctUrban + PctBlack + PctHisp + PctHighSchool +
                  MedianHHInc + PctPoor + PctFemale +
                  PctOccupied + PctMovedIn5 + MedianHValue +
                  smokerate2000 + avgdewpt + avgtemp + avgrelhum + Population, 
                  data = data_temp, method = "nearest", distance = "glm", replace = FALSE)

matched_data <- match.data(m_temp, data = data_temp)
matched_rnames <- row.names(matched_data)

# loveplot
data_love <- matched_data %>%
  select(a, 
         PctUrban, PctBlack, PctHisp, PctHighSchool,
         MedianHHInc, PctPoor, PctFemale,
         PctOccupied, PctMovedIn5, MedianHValue,
         smokerate2000, avgdewpt, avgtemp, avgrelhum, Population)
w_out <- weightit(a ~ ., data = data_love, estimand = "ATE", method = "ps")
love.plot(bal.tab(w_out), binary = 'std', thresholds = c(m = .1))

pdf(file="love_plot.pdf")
love.plot(bal.tab(w_out), binary = 'std', thresholds = c(m = .1))
dev.off()

###########################################################################
# ----   get_matches()     ----
###########################################################################

match_ <- get_matches(m_temp)
dataset_matched <- match_[,c(4:24)]

#Save
save(dataset_matched, file = "dataset_matched.RData")

