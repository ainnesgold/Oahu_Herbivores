##Time
timesteps <- 100
##Area
patch_area <- c(0.1, 0.9)
number_patches <- length(patch_area)

fishing_effort <- c(0.1, 0)
catchability <- 1

intrinsic_growth_rate <- 0.3
carrying_capacity <- c(350, 350)

source('seafood_prices.R')
price <- prices_sum$PricePerG_2021
value_added_ratio <- 0.9  #from Grafeld et al. 2017

S <- 0.8
S_recruit <- 0.1
source('dispersal_fn.R')

##Grid of all different possible combinations


## Population Dynamics Functions ------------------------------------------------------------------------------------------------

source('pop_dynamics_fns.R')

## State Variables -----------------------------------------------------------------------------------------------


tmp1 <- dispersal(patch_area, number_patches, S)
tmp2 <- dispersal(patch_area, number_patches, S_recruit)

source('initialize_state_vars.R')

## Run Population Model ---------------------------------------------------------------------------------------------------------

for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
  for(i in 1:number_patches){
    population[t, i] <- calculate_population_growth(population[t-1, i], intrinsic_growth_rate, 
                                                    carrying_capacity[i])
    recruits[t, i] <- population[t, i] - population[t-1, i]
    adults[t, i] <- population[t, i] - recruits[t, i]
  }
  
  recruits[t,] <- recruits[t,] %*% tmp2
  
  for(i in 1:number_patches){
    # harvest
    fraction_harvested[t, i] <- calculate_fraction_harvested(fishing_effort[i],
                                                             catchability)
    
    harvest[t, i] <- calculate_fisheries_harvest(adults[t, i], fraction_harvested[t, i], 
                                                 (patch_area[i]*5.04e+8)) 
    
    # escapement
    escapement[t, i] <- calculate_escaped_stock_biomass(adults[t, i], fraction_harvested[t, i])
  }  
  
  
  adults[t,] <- escapement[t,] %*% tmp1
  
  revenue <- (price * harvest) * value_added_ratio
  for(i in 1:number_patches){
    population[t, i] <- recruits[t, i] + adults[t, i]
    if (population[t, i] > carrying_capacity[i]) {
      population[t, i] = carrying_capacity[i]
    }
  }
}  

par(mfrow=c(1,2))
plot(population[,1])
plot(population[,2])