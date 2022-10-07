
library(purrr)

##Time
timesteps <- 50
##Area
patch_area_sequences <- list(seq(0, 1, by = 0.1))
patch_area_grid <- do.call(expand.grid, patch_area_sequences)
patch_area_grid$Var2 <- 1 - patch_area_grid$Var1
patch_area_list <- split(patch_area_grid, 1:nrow(patch_area_grid))
number_patches <- ncol(patch_area_grid)

fishing_effort_sequences <- list(seq(0, 1, by = 0.1), 0)
fishing_effort_grid <- do.call(expand.grid, fishing_effort_sequences)
fishing_effort_list <- split(fishing_effort_grid, 1:nrow(fishing_effort_grid))
catchability <- 1

intrinsic_growth_rate <- 0.3
carrying_capacity <- 350

source('seafood_prices.R')
price <- prices_sum$PricePerG_2021
value_added_ratio <- 0.9  #from Grafeld et al. 2017

S <- seq(0, 1, by = 0.1)
S_recruit <- seq(0, 1, by = 0.1)
source('dispersal_fn.R')

##Grid of all different possible combinations
parameter_grid <- expand.grid(S = S, S_recruit = S_recruit,
                              patch_area = patch_area_list,
                              fishing_effort = fishing_effort_list)


## Population Dynamics Functions ------------------------------------------------------------------------------------------------

source('pop_dynamics_fns.R')

## State Variables -----------------------------------------------------------------------------------------------

outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))


for (iter in 1:nrow(parameter_grid)) {
  
  tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, parameter_grid[iter, "S"])
  tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, parameter_grid[iter, "S_recruit"])
  
  source('initialize_state_vars.R')
  
  ## Run Population Model ---------------------------------------------------------------------------------------------------------
  
  for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
    for(i in 1:number_patches){
      population[t, i] <- calculate_population_growth(population[t-1, i], intrinsic_growth_rate, 
                                                      carrying_capacity)
      recruits[t, i] <- population[t, i] - population[t-1, i]
      adults[t, i] <- population[t, i] - recruits[t, i]
    }
    recruits_dispersal[t,] <- recruits[t,] %*% tmp2
    
    for(i in 1:number_patches){
      # harvest
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      
      harvest[t, i] <- calculate_fisheries_harvest(adults[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8)) 
      
      # escapement
      escapement[t, i] <- calculate_escaped_stock_biomass(adults[t, i], fraction_harvested[t, i])
    }  
    
    # Dispersal
    adults <- escapement[,] %*% tmp1
    
    for (i in 1:number_patches){
      population[t, i] <- adults[t, i] + recruits_dispersal[t, i]
      if (population[t, i] > carrying_capacity) {
        population[t, i] == carrying_capacity
      }
    }
    #calculate market equivalent revenue
    revenue <- (price * harvest) * value_added_ratio
  }  
  
  outcome_biomass[iter, ] <- population[t, ]
  outcome_harvest[iter, ] <- harvest[t, ]
  
}


# cbind outcome with parameter grid
colnames(outcome_biomass) <- c("patch_1_population", "patch_2_population")
colnames(outcome_harvest) <- c("patch_1_harvest", "patch_2_harvest")

outcome <- cbind(parameter_grid, outcome_biomass, outcome_harvest)

outcome$fishing_effort_patch1 <- map_dbl(outcome$fishing_effort, 1)
outcome$area_patch1 <- map_dbl(outcome$patch_area, 1)
outcome$area_patch2 <- map_dbl(outcome$patch_area, 2)


