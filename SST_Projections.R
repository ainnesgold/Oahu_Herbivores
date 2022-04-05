
# Code from: https://rpubs.com/markpayne/358146
#Dataset was too large to upload to github but can be downloaded here
#CMIP6: monthly, single levels, SST, CESM2 (USA), whole temporal range, whole available region
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cmip6?tab=form

library(tidyverse)
library(raster)
library(ncdf4)

#2015-2064
##SSP2-4.5
fname <- "Temp options/CMIP6 trials/4.5/tos_Omon_CESM2_ssp245_r4i1p1f1_gn_201501-206412_v20200528.nc"

temp_file <- brick(fname)

lon.pts <- seq(202,202.99,by=0.01) #~long of Oahu edges (calculated as 360 - lon to get degrees east). 202-202.99 by 0.01
lat.pts <- seq(201.5, 202.49, by=0.01)  #degrees north. equator is 180 + oahu latittude edges (21) = 201, 201.99
extract.pts <- cbind(lon.pts,lat.pts)

##LOOP TO CALCULATE AVG TEMP FOR EACH YEAR 
number_years <- (nlayers(temp_file) / 12.0) - 1
years <- array(2015:2064)
mean_temps <- c()

for (t in 0:number_years) {
  year_step <- t * 12
  b <- temp_file[[t+1:t+12]]
  ann.extract <- extract(b,extract.pts,method="bilinear")
  mean_temp <- mean(ann.extract)
  mean_temps <- append(mean_temps, mean_temp)
}

mean_temps<- cbind(years, mean_temps)


#2065-2100
fname2 <- "Temp options/CMIP6 trials/4.5/tos_Omon_CESM2_ssp245_r4i1p1f1_gn_206501-210012_v20200528.nc"
temp_file2 <- brick(fname2)

##LOOP TO CALCULATE AVG TEMP FOR EACH YEAR 
number_years <- (nlayers(temp_file2) / 12.0) - 1
years <- array(2065:2100)
mean_temps2 <- c()

for (t in 0:number_years) {
  year_step <- t * 12
  b <- temp_file2[[t+1:t+12]]
  ann.extract <- extract(b,extract.pts,method="bilinear")
  mean_temp <- mean(ann.extract)
  mean_temps2 <- append(mean_temps2, mean_temp)
}

mean_temps2<- cbind(years, mean_temps2)

mean_temps<-rbind(mean_temps, mean_temps2)
mean_temps <- as.data.frame(mean_temps)



#Calculate future years deviation from present
#historical data, too big for github but can be downloaded from:
#https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html HadISST netCDF at bottom of page
lon.pts <- seq(157.5, 158.5, by = 1)
lat.pts <- seq(21, 22, by = 1)
extract.pts <- cbind(lon.pts,lat.pts)

histname <- "Temp options/CMIP6 trials/HadISST_sst.nc"
hist_file <- brick(histname)

##LOOP TO CALCULATE AVG TEMP FOR EACH YEAR 
hist_number_years <- (nlayers(hist_file) / 12.0)
years <- array(1870:2021)
hist_mean_temps <- c()

for (t in 0:hist_number_years) {
  year_step <- t * 12
  b <- hist_file[[t+1:t+12]]
  ann.extract <- extract(b,extract.pts,method="bilinear")
  hist_mean_temp <- mean(ann.extract)
  hist_mean_temps <- append(hist_mean_temps, hist_mean_temp)
}

hist_mean_temps<- cbind(years, hist_mean_temps)

hist_mean_temps<-as.data.frame(hist_mean_temps)

current_temp <- hist_mean_temps %>%
  filter(years >= 1991,
         years <= 2021) %>%
  summarize(avgtemp = mean(hist_mean_temps))

#create dataframe of future values only, starting at 2021. I want 50 years right now to match the model runs.
mean_temps_crop <- mean_temps %>%
  filter(years > 2020, years < 2071)

projections <- mean_temps %>%
  filter(years > 2020, years < 2071) %>%
  summarise(anomaly = mean_temps - current_temp$avgtemp)

rcp4.5<-mean_temps

#supplemental temp figure
ggplot(data=rcp4.5 %>% filter (years > 2019, years < 2071), aes(x=years, y=mean_temps)) +
  geom_point() +
  geom_line() +
  labs(x="Year", y="°C") +
  theme_minimal()

