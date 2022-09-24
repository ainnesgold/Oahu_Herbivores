#####################################POPULATION GROWTH########################################################################
##Population growth - Ricker model - no temp
calculate_population_growth <- function(population, intrinsic_growth_rate, carrying_capacity) {
  population <-population * exp(intrinsic_growth_rate * (1 - population / carrying_capacity))
  return(population)
}

#Fraction of stock harvested
calculate_fraction_harvested <- function(fishing_effort, catchability) {
  fraction_harvested <- 1 - exp(-fishing_effort * catchability)
  return(fraction_harvested)
}
#Fisheries Harvest
calculate_fisheries_harvest <- function(population, fraction_harvested, patch_area_m2) {
  harvest <- population * fraction_harvested * patch_area_m2
  return(harvest)
}

#Escaped stock biomass
calculate_escaped_stock_biomass <- function(population, fraction_harvested) {
  escaped_stock_biomass <- population * (1 - fraction_harvested)
  return(escaped_stock_biomass)
}

###################################THERMAL EFFECTS############################################################################
#Climate 1: temp dependent r
a = 0.3
b = 0
c = -0.0037 

calculate_r_temp <- function(SST_dev, a, b, c) {
  r_temp <- a + b*SST_dev + c*SST_dev^2
  return(r_temp)
}

calculate_population_growth_temp_r <- function(population, r_temp, carrying_capacity) {
  population <-population * exp(r_temp * (1 - population / carrying_capacity))
  return(population)
}

#Climate 2: temp dependent K
#m = -10.76933
#d = 635.5131
m=-7.234
d=541.932
calculate_CC_temp <- function(mean_temp_series) {
  CC_temp <- (m * mean_temp_series) + d
  return(CC_temp)
}

calculate_population_growth_temp_CC <- function(population, intrinsic_growth_rate, CC_temp) {
  population <-population * exp(intrinsic_growth_rate * (1 - population / CC_temp))
  return(population)
}

#final equation using both R and K
calculate_population_growth_temp_r_CC <- function(population, r_temp, CC_temp) {
  population <-population * exp(r_temp * (1 - population / CC_temp))
  return(population)
}
#####################################################################################################################
