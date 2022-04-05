#read in new csv
#get average price
#convert price per lb to price per gram
#https://www.fisheries.noaa.gov/inport/item/29742
prices<-read.csv("ESD_RetailMonitoring2016_v01.csv")

library(tidyverse)

#all reef to see species
prices_sum <- prices %>%
  filter(GROUP == "REEF") %>%
  group_by(SPECIES) %>%
  summarise(AvgPrice = mean(na.omit(PRICE)))

#get herbivores
#https://dlnr.hawaii.gov/marine30x30/herbivoremanagement/participate-in-the-process/

herbivores <- c("A'AWA", "KOLE", "MANINI", "NAENAE", "NENUE (CHUB)", "PALANI", #"TABLE BOSS (HOGFISH - BLACKSPOT WRASSE)",
                "UHU (BLUE - MALE PARROTFISH)", "UHU (FEMALE - UNSPECIFIED)", "UHU (PARROTFISH)", "UKU ULIULI (SPECTACLED PARROTFISH)")

#average price per lb for all reef herbivores
prices_sum <- prices %>%
  filter( SPECIES %in% herbivores) %>%
  summarise(AvgPricePerLb = mean(na.omit(PRICE)))

prices_sum$AvgPricePerG<-prices_sum$AvgPricePerLb / 453.592 #number of grams in a lb

#then use this inflation calculator to get 2021 value: $6.57/lb
prices_sum$PricePerLB_2021 <- 6.57
prices_sum$PricePerG_2021 <- prices_sum$PricePerLB_2021 / 453.592



