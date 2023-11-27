
# dataset
load("Dataset_matched.RData")
load("Dataset_Policies.RData")

# data -- ZIP CODES
library(ggplot2)
library(sf)
library(data.table)
library(tidyverse)

zip_2010 <- st_read("/n/home_fasse/dzorzetto/code_bnp_pm25/zip_code/ESRI10USZIP5_POLY_WGS84.shp")

unique(zip_2010$STATE)

zip_2010 <- zip_2010 %>%
  filter(STATE != "AK")

zip_2010 <- zip_2010 %>%
  filter(STATE != "HI")

zip_2010 <- zip_2010 %>%
  filter(STATE != "PR")

ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  geom_point(data = dataset_matched, aes(x = Longitude, y = Latitude),
             col="black", size=2)+
  theme_void() +
  theme(legend.position = "right")

dataset_matched$Attainment=as.factor(dataset_matched$a)
data_merged$Attainment=as.factor(data_merged$a)

pdf(file="attainment_standard.pdf", height=4)
#jpeg(file="map_attainment_standard.jpeg", quality = 150, height=300)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Attainment of the Air Pollution Standard") + 
  #geom_point(data = dataset_matched, 
  geom_point(data = data_merged,
             aes(x = Longitude, y = Latitude, color=Attainment, fill=Attainment),
              size=2)+
  theme_void() +
  scale_color_manual(values=c("#00E01A", "#E01E00")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.92,.2))
  #theme(legend.position = "none")
dev.off()

pdf(file="pollution.pdf", height=4)
#jpeg(file="map_pollution.jpeg", quality = 150, height=300)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Variation of PM2.5") + 
  #geom_point(data = dataset_matched, 
  geom_point(data = data_merged,
             aes(x = Longitude, y = Latitude, color=pm_diff),
             size=2)+
  theme_void() +
  scale_color_gradient(low="#FADE00", high="#FA3700") +
  labs(color = " ") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()

pdf(file="mortality.pdf", height=4)
#jpeg(file="map_mortality.jpeg", quality = 150, height=300)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Variation of Mortality Rate") + 
  #geom_point(data = dataset_matched, 
  geom_point(data = data_merged,
             aes(x = Longitude, y = Latitude, color=all_causes_ADJ),
             size=2)+
  theme_void() +
  scale_color_gradient(low="#00FA5E", high="#300DFF") +
  labs(color = " ") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()


################################################################################

load("Dataset_merged.RData")
load("results_application.RData")

data_merged$prob_1pos=apply(results$chians_strata==1,2,mean)
data_merged$prob_1neg=apply(results$chians_strata==(-1),2,mean)
data_merged$prob_0=apply(results$chians_strata==0,2,mean)

pdf(file="prob_1pos.pdf", height=4)
#jpeg(file="map_mortality.jpeg", quality = 150, height=300)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Probability of associative negative stratum") + 
  geom_point(data = data_merged, 
             aes(x = Longitude, y = Latitude, color=prob_1pos),
             size=2)+
  theme_void() +
  scale_color_gradient(low="white", high="#D90224") +
  labs(color = "prob") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()


pdf(file="prob_1neg.pdf", height=4)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Probability of associative positive stratum") + 
  geom_point(data = data_merged, 
             aes(x = Longitude, y = Latitude, color=prob_1neg),
             size=2)+
  theme_void() +
  scale_color_gradient(low="white", high="#35D90B") +
  labs(color = "prob") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()

pdf(file="prob_0.pdf", height=4)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Probability of dissociative stratum") + 
  geom_point(data = data_merged, 
             aes(x = Longitude, y = Latitude, color=prob_0),
             size=2)+
  theme_void() +
  scale_color_gradient(low="white", high="#F0D400") +
  labs(color = "prob") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()

data_merged$strata=as.factor(results$S_strata_cluster)

pdf(file="strata.pdf", height=4)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Strata allocation") + 
  #geom_point(data = dataset_matched, 
  geom_point(data = data_merged,
             aes(x = Longitude, y = Latitude, color=strata, fill=strata),
             size=2)+
  theme_void() +
  scale_color_manual(values=c("#35D90B", "#F0D400","#D90224")) +
  labs(color = "strata") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.92,.2))
#theme(legend.position = "none")
dev.off()