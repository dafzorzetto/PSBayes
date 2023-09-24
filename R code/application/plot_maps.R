
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

pdf(file="attainment_standard.pdf", height=4)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Attainment of the Air Pollution Standard") + 
  geom_point(data = dataset_matched, 
             aes(x = Longitude, y = Latitude, color=Attainment, fill=Attainment),
              size=2)+
  theme_void() +
  scale_color_manual(values=c("#00E01A", "#E01E00")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.92,.2))
  #theme(legend.position = "none")
dev.off()

pdf(file="pollution.pdf", height=4)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Variation of PM2.5") + 
  geom_point(data = dataset_matched, 
             aes(x = Longitude, y = Latitude, color=pm_diff),
             size=2)+
  theme_void() +
  scale_color_gradient(low="#FADE00", high="#FA3700") +
  labs(color = " ") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()

pdf(file="mortality.pdf", height=4)
ggplot()+ 
  geom_sf(data = zip_2010$geometry,  col = NA)  +
  ggtitle("Variation of Mortality Rate") + 
  geom_point(data = dataset_matched, 
             aes(x = Longitude, y = Latitude, color=all_causes_ADJ),
             size=2)+
  theme_void() +
  scale_color_gradient(low="#00FA5E", high="#300DFF") +
  labs(color = " ") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.95,.3))
dev.off()
