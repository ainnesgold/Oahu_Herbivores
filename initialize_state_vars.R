# State variables
population <- array(NA, dim = c(timesteps, number_patches))
recruits <- array(NA, dim = c(timesteps, number_patches)) 
recruits_dispersal <- array(NA, dim = c(timesteps, number_patches)) 

# set initial starting population in each patch
starting_population <- c(14,14) #14g/m2 median herbivore biomass value for Oahu - from Donovan et al. in prep
population[1, ] <- starting_population

#fisheries
fraction_harvested <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))

#escaped biomass
escapement <- array(NA, dim = c(timesteps, number_patches))

#temp dependent r
r_temp <- array(NA, dim = c(timesteps, number_patches)) 

#temp dependent CC
#temps_percent_dev <- array(NA, dim = c(timesteps, number_patches)) 
CC_temp <- array(NA, dim = c(timesteps, number_patches)) 
