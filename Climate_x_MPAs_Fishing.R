library(tidyverse)
library(viridis)
library(ggpubr)
library(grid)
library(purrr)

#CLIMATE 
source('SST_Projections.R')
SST_dev <- projections[['anomaly']]
mean_temps_series <- mean_temps_crop[['extracted_mean_temps']]

##Time
timesteps <- 50

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
carrying_capacity <- c(342.4, 342.4)

source('seafood_prices.R')
price <- prices_sum$PricePerG_2021
value_added_ratio <- 0.9  #from Grafeld et al. 2017

S <- 0.7
S_recruit <- 0.1
source('dispersal_fn.R')

##Grid of all different possible combinations
parameter_grid <- expand.grid(patch_area = patch_area_list,
                              fishing_effort = fishing_effort_list)

## Population Dynamics Functions ------------------------------------------------------------------------------------------------

source('pop_dynamics_fns.R')

## Parameters and State Variables -----------------------------------------------------------------------------------------------
outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_revenue <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))



for (iter in 1:nrow(parameter_grid)) {
  
  tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S_recruit)
  
  source('initialize_state_vars.R')
  
  ## Run Population Model ---------------------------------------------------------------------------------------------------------
  ##fish population dynamics
  
  for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
    for(i in 1:number_patches){
      population[t, i] <- calculate_population_growth(population[t-1, i], intrinsic_growth_rate, 
                                                      carrying_capacity[i])
      recruits[t, i] <- population[t, i] - population[t-1, i]
    }
    recruits_dispersal[t,] <- recruits[t,] %*% tmp2
    
    for(i in 1:number_patches){
      population[t, i] <- population[t, i] - recruits[t, i] + recruits_dispersal[t, i]
      # harvest
      # use function to calculate harvest
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      # We will have to figure out a way to scale fishing effort because the negative exponential of a high number like 10000, or even
      # just 10, gives a pretty small number. 
      # try running this code -- this shows that at fishing effort = 0, the survival is 1
      # x <- seq(0, 10, length.out = 100)
      # y <- exp(-x)
      # plot(x,y)
      
      harvest[t, i] <- calculate_fisheries_harvest(population[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8)) 
      
      # escapement
      # subtract the harvest at time t and patch i from the population[t, i]
      escapement[t, i] <- calculate_escaped_stock_biomass(population[t, i], fraction_harvested[t, i])
    }  
    
    # After the above population dynamics run for a timestep, we want the fish to disperse between patches.
    # At time t, have the dispersal function determine:
    # 1. the number (or biomass) of fish in patch i=1 that stay in patch i=1 and the number that move to patch i=2
    # 2. the number in patch i=2 that stay in i=2 and the number that move to i=1.
    population <- escapement[,] %*% tmp1
    #calculate market equivalent revenue
    revenue <- (price * harvest) * value_added_ratio
    # after the above population dynamics run for time t and each patch, the next timestep (t+1) will proceed using values from time t.
  }  
  
  outcome_biomass[iter, ] <- population[t, ]
  outcome_harvest[iter, ] <- harvest[t, ]
  outcome_revenue[iter, ] <- revenue[t, ]
}


# cbind outcome with parameter grid
colnames(outcome_biomass) <- c("patch_1_population", "patch_2_population")
colnames(outcome_harvest) <- c("patch_1_harvest", "patch_2_harvest")

outcome <- cbind(parameter_grid, outcome_biomass, outcome_harvest)

outcome$fishing_effort_patch1 <- map_dbl(outcome$fishing_effort, 1)
outcome$area_patch1 <- map_dbl(outcome$patch_area, 1)
outcome$area_patch2 <- map_dbl(outcome$patch_area, 2)

outcome$area_patch1_m <- outcome$area_patch1 * 5.04e+8
outcome$area_patch2_m <- outcome$area_patch2 * 5.04e+8
outcome$total_pop_p1 <- (outcome$patch_1_population * outcome$area_patch1_m) /1000
outcome$total_pop_p2 <- (outcome$patch_2_population * outcome$area_patch2_m) /1000
outcome$total_pop <- outcome$total_pop_p1 + outcome$total_pop_p2
outcome$CPUE <- (outcome$patch_1_harvest/1000) / outcome$fishing_effort_patch1
outcome$CPUE[is.na(outcome$CPUE)] <- 0

outcome <- outcome[c("patch_1_harvest","total_pop_p1","total_pop_p2", "total_pop", "CPUE", 
                     "area_patch2", "fishing_effort_patch1")]
outcome$patch_1_harvest <- outcome$patch_1_harvest/1000

outcome_long <- outcome %>%
  pivot_longer(cols = patch_1_harvest:CPUE, names_to = "patch", values_to = "population")
outcome_long$Patch2_Percent <- outcome_long$area_patch2 * 100



########################################################first climate run
outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_revenue <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))



for (iter in 1:nrow(parameter_grid)) {
  
  tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S_recruit)
  
  source('initialize_state_vars.R')
  
  ## Run Population Model ---------------------------------------------------------------------------------------------------------
  ##fish population dynamics
  
  for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
    for(i in 1:number_patches){
      r_temp[t, i] <- calculate_r_temp(SST_dev[t], a, b, c)
      population[t, i] <- calculate_population_growth_temp_r(population[t-1, i], r_temp[t, i], 
                                                             carrying_capacity[i])
      recruits[t, i] <- population[t, i] - population[t-1, i]
    }
    recruits_dispersal[t,] <- recruits[t,] %*% tmp2
    
    for(i in 1:number_patches){
      population[t, i] <- population[t, i] - recruits[t, i] + recruits_dispersal[t, i]
      # harvest
      # use function to calculate harvest
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      # We will have to figure out a way to scale fishing effort because the negative exponential of a high number like 10000, or even
      # just 10, gives a pretty small number. 
      # try running this code -- this shows that at fishing effort = 0, the survival is 1
      # x <- seq(0, 10, length.out = 100)
      # y <- exp(-x)
      # plot(x,y)
      
      harvest[t, i] <- calculate_fisheries_harvest(population[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8)) 
      
      # escapement
      # subtract the harvest at time t and patch i from the population[t, i]
      escapement[t, i] <- calculate_escaped_stock_biomass(population[t, i], fraction_harvested[t, i])
    }  
    
    # After the above population dynamics run for a timestep, we want the fish to disperse between patches.
    # At time t, have the dispersal function determine:
    # 1. the number (or biomass) of fish in patch i=1 that stay in patch i=1 and the number that move to patch i=2
    # 2. the number in patch i=2 that stay in i=2 and the number that move to i=1.
    population <- escapement[,] %*% tmp1
    #calculate market equivalent revenue
    revenue <- (price * harvest) * value_added_ratio
    # after the above population dynamics run for time t and each patch, the next timestep (t+1) will proceed using values from time t.
  }  
  
  outcome_biomass[iter, ] <- population[t, ]
  outcome_harvest[iter, ] <- harvest[t, ]
  outcome_revenue[iter, ] <- revenue[t, ]
}

##
colnames(outcome_biomass) <- c("patch_1_population", "patch_2_population")
colnames(outcome_harvest) <- c("patch_1_harvest", "patch_2_harvest")

outcome_climate1 <- cbind(parameter_grid, outcome_biomass, outcome_harvest)
outcome_climate1$fishing_effort_patch1 <- map_dbl(outcome_climate1$fishing_effort, 1)
outcome_climate1$area_patch1 <- map_dbl(outcome_climate1$patch_area, 1)
outcome_climate1$area_patch2 <- map_dbl(outcome_climate1$patch_area, 2)

outcome_climate1$area_patch1_m <- outcome_climate1$area_patch1 * 5.04e+8
outcome_climate1$area_patch2_m <- outcome_climate1$area_patch2 * 5.04e+8
outcome_climate1$total_pop_p1 <- (outcome_climate1$patch_1_population * outcome_climate1$area_patch1_m) /1000
outcome_climate1$total_pop_p2 <- (outcome_climate1$patch_2_population * outcome_climate1$area_patch2_m) /1000
outcome_climate1$total_pop <- outcome_climate1$total_pop_p1 + outcome_climate1$total_pop_p2
outcome_climate1$CPUE <- (outcome_climate1$patch_1_harvest/1000) / outcome_climate1$fishing_effort_patch1
outcome_climate1$CPUE[is.na(outcome_climate1$CPUE)] <- 0


outcome_climate1 <- outcome_climate1[c("patch_1_harvest","total_pop_p1","total_pop_p2", "total_pop", "CPUE",
                                       "area_patch2", "fishing_effort_patch1")]
outcome_climate1$patch_1_harvest <- outcome_climate1$patch_1_harvest/1000
outcome_climate1_long <- outcome_climate1 %>%
  pivot_longer(cols = patch_1_harvest:CPUE, names_to = "patch", values_to = "population")
outcome_climate1_long$Patch2_Percent <- outcome_climate1_long$area_patch2 * 100






########################################################second climate run
outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_revenue <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))



for (iter in 1:nrow(parameter_grid)) {
  
  tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S_recruit)
  
  source('initialize_state_vars.R')
  
  ## Run Population Model ---------------------------------------------------------------------------------------------------------
  ##fish population dynamics
  
  for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
    for(i in 1:number_patches){
      CC_temp[t, i] <- calculate_CC_temp(mean_temps_series[t])
      population[t, i] <- calculate_population_growth_temp_CC(population[t-1, i], intrinsic_growth_rate, 
                                                              CC_temp[t, i])
      recruits[t, i] <- population[t, i] - population[t-1, i]
    }
    recruits_dispersal[t,] <- recruits[t,] %*% tmp2
    
    for(i in 1:number_patches){
      population[t, i] <- population[t, i] - recruits[t, i] + recruits_dispersal[t, i]
      # harvest
      # use function to calculate harvest
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      # We will have to figure out a way to scale fishing effort because the negative exponential of a high number like 10000, or even
      # just 10, gives a pretty small number. 
      # try running this code -- this shows that at fishing effort = 0, the survival is 1
      # x <- seq(0, 10, length.out = 100)
      # y <- exp(-x)
      # plot(x,y)
      
      harvest[t, i] <- calculate_fisheries_harvest(population[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8)) 
      
      # escapement
      # subtract the harvest at time t and patch i from the population[t, i]
      escapement[t, i] <- calculate_escaped_stock_biomass(population[t, i], fraction_harvested[t, i])
    }  
    
    # After the above population dynamics run for a timestep, we want the fish to disperse between patches.
    # At time t, have the dispersal function determine:
    # 1. the number (or biomass) of fish in patch i=1 that stay in patch i=1 and the number that move to patch i=2
    # 2. the number in patch i=2 that stay in i=2 and the number that move to i=1.
    population <- escapement[,] %*% tmp1
    #calculate market equivalent revenue
    revenue <- (price * harvest) * value_added_ratio
    # after the above population dynamics run for time t and each patch, the next timestep (t+1) will proceed using values from time t.
  }  
  
  outcome_biomass[iter, ] <- population[t, ]
  outcome_harvest[iter, ] <- harvest[t, ]
  outcome_revenue[iter, ] <- revenue[t, ]
}

##
colnames(outcome_biomass) <- c("patch_1_population", "patch_2_population")
colnames(outcome_harvest) <- c("patch_1_harvest", "patch_2_harvest")

outcome_climate2 <- cbind(parameter_grid, outcome_biomass, outcome_harvest)
outcome_climate2$fishing_effort_patch1 <- map_dbl(outcome_climate2$fishing_effort, 1)
outcome_climate2$area_patch1 <- map_dbl(outcome_climate2$patch_area, 1)
outcome_climate2$area_patch2 <- map_dbl(outcome_climate2$patch_area, 2)

outcome_climate2$area_patch1_m <- outcome_climate2$area_patch1 * 5.04e+8
outcome_climate2$area_patch2_m <- outcome_climate2$area_patch2 * 5.04e+8
outcome_climate2$total_pop_p1 <- (outcome_climate2$patch_1_population * outcome_climate2$area_patch1_m) /1000
outcome_climate2$total_pop_p2 <- (outcome_climate2$patch_2_population * outcome_climate2$area_patch2_m) /1000
outcome_climate2$total_pop <- outcome_climate2$total_pop_p1 + outcome_climate2$total_pop_p2
outcome_climate2$CPUE <- (outcome_climate2$patch_1_harvest/1000) / outcome_climate2$fishing_effort_patch1
outcome_climate2$CPUE[is.na(outcome_climate2$CPUE)] <- 0

outcome_climate2 <- outcome_climate2[c("patch_1_harvest","total_pop_p1","total_pop_p2", "total_pop", "CPUE", 
                                       "area_patch2", "fishing_effort_patch1")]
outcome_climate2$patch_1_harvest <- outcome_climate2$patch_1_harvest/1000
outcome_climate2_long <- outcome_climate2 %>%
  pivot_longer(cols = patch_1_harvest:CPUE, names_to = "patch", values_to = "population")
outcome_climate2_long$Patch2_Percent <- outcome_climate2_long$area_patch2 * 100




########################################################third climate run
outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_revenue <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))



for (iter in 1:nrow(parameter_grid)) {
  
  tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S_recruit)
  
  source('initialize_state_vars.R')
  
  ## Run Population Model ---------------------------------------------------------------------------------------------------------
  ##fish population dynamics
  
  for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
    for(i in 1:number_patches){
      r_temp[t, i] <- calculate_r_temp(SST_dev[t], a, b, c)
      CC_temp[t, i] <- calculate_CC_temp(mean_temps_series[t])
      population[t, i] <- calculate_population_growth_temp_r_CC(population[t-1, i], r_temp[t, i], 
                                                                CC_temp[t, i])
      recruits[t, i] <- population[t, i] - population[t-1, i]
    }
    recruits_dispersal[t,] <- recruits[t,] %*% tmp2
    
    for(i in 1:number_patches){
      population[t, i] <- population[t, i] - recruits[t, i] + recruits_dispersal[t, i]
      # harvest
      # use function to calculate harvest
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      # We will have to figure out a way to scale fishing effort because the negative exponential of a high number like 10000, or even
      # just 10, gives a pretty small number. 
      # try running this code -- this shows that at fishing effort = 0, the survival is 1
      # x <- seq(0, 10, length.out = 100)
      # y <- exp(-x)
      # plot(x,y)
      
      harvest[t, i] <- calculate_fisheries_harvest(population[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8)) 
      
      # escapement
      # subtract the harvest at time t and patch i from the population[t, i]
      escapement[t, i] <- calculate_escaped_stock_biomass(population[t, i], fraction_harvested[t, i])
    }  
    
    # After the above population dynamics run for a timestep, we want the fish to disperse between patches.
    # At time t, have the dispersal function determine:
    # 1. the number (or biomass) of fish in patch i=1 that stay in patch i=1 and the number that move to patch i=2
    # 2. the number in patch i=2 that stay in i=2 and the number that move to i=1.
    population <- escapement[,] %*% tmp1
    #calculate market equivalent revenue
    revenue <- (price * harvest) * value_added_ratio
    # after the above population dynamics run for time t and each patch, the next timestep (t+1) will proceed using values from time t.
  }  
  
  outcome_biomass[iter, ] <- population[t, ]
  outcome_harvest[iter, ] <- harvest[t, ]
  outcome_revenue[iter, ] <- revenue[t, ]
}

##
colnames(outcome_biomass) <- c("patch_1_population", "patch_2_population")
colnames(outcome_harvest) <- c("patch_1_harvest", "patch_2_harvest")

outcome_climate3 <- cbind(parameter_grid, outcome_biomass, outcome_harvest)
outcome_climate3$fishing_effort_patch1 <- map_dbl(outcome_climate3$fishing_effort, 1)
outcome_climate3$area_patch1 <- map_dbl(outcome_climate3$patch_area, 1)
outcome_climate3$area_patch2 <- map_dbl(outcome_climate3$patch_area, 2)

outcome_climate3$area_patch1_m <- outcome_climate3$area_patch1 * 5.04e+8
outcome_climate3$area_patch2_m <- outcome_climate3$area_patch2 * 5.04e+8
outcome_climate3$total_pop_p1 <- (outcome_climate3$patch_1_population * outcome_climate3$area_patch1_m) /1000
outcome_climate3$total_pop_p2 <- (outcome_climate3$patch_2_population * outcome_climate3$area_patch2_m) /1000
outcome_climate3$total_pop <- outcome_climate3$total_pop_p1 + outcome_climate3$total_pop_p2
outcome_climate3$CPUE <- (outcome_climate3$patch_1_harvest/1000) / outcome_climate3$fishing_effort_patch1
outcome_climate3$CPUE[is.na(outcome_climate3$CPUE)] <- 0

outcome_climate3 <- outcome_climate3[c("patch_1_harvest","total_pop_p1","total_pop_p2", "total_pop", "CPUE",
                                       "area_patch2", "fishing_effort_patch1")]
outcome_climate3$patch_1_harvest <- outcome_climate3$patch_1_harvest/1000
outcome_climate3_long <- outcome_climate3 %>%
  pivot_longer(cols = patch_1_harvest:CPUE, names_to = "patch", values_to = "population")
outcome_climate3_long$Patch2_Percent <- outcome_climate3_long$area_patch2 * 100



patch_population_compare <- cbind(outcome_long$Patch2_Percent, outcome_long$fishing_effort_patch1,
                                  outcome_long$patch, outcome_long$population, 
                                  outcome_climate1_long$population, outcome_climate2_long$population,
                                  outcome_climate3_long$population)



patch_population_compare <- as.data.frame(patch_population_compare)

colnames(patch_population_compare)[1] <- "Patch2Percent"
colnames(patch_population_compare)[2] <- "FishingEffort"
colnames(patch_population_compare)[3] <- "Patch"
colnames(patch_population_compare)[4] <- "NoClimate"
colnames(patch_population_compare)[5] <- "Climate1"
colnames(patch_population_compare)[6] <- "Climate2"
colnames(patch_population_compare)[7] <- "Climate3"



patch_population_compare <- patch_population_compare %>%
  pivot_longer(cols = NoClimate:Climate3, names_to = "Scenario", values_to = "Biomass")

patch_population_compare$Patch <- factor(patch_population_compare$Patch, levels = c("total_pop_p1", "total_pop_p2",
                                                                                    "total_pop", "patch_1_harvest", "CPUE"),
                                         labels = c("Fishable Biomass", "Protected Biomass", "Total Biomass", "Total Harvest", "CPUE"))

patch_population_compare$Scenario <- factor(patch_population_compare$Scenario, levels = c("NoClimate", "Climate1", "Climate2", "Climate3"),
                                            labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", 
                                                       "Temp-dependent r + K"))

#Main text plot
#biomass
p1 <- ggplot(data=patch_population_compare %>% filter(Patch=="Total Biomass",
                                                      FishingEffort == 0.1 | FishingEffort ==0.3 | FishingEffort == 0.7),
             aes(x = as.numeric(Patch2Percent), y = as.numeric(Biomass), color = Scenario)) +
  geom_line(aes(size = as.factor(Scenario))) +
  facet_wrap(~FishingEffort) +
  scale_color_viridis_d(option="turbo", name="", 
                        breaks = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K"),
                        labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K")) +
  labs(x="% Coastal Waters Closed to Fishing", y="Total Fish Biomass (kg)") +
  scale_size_manual(name = "", values = c(1, 1, 1, 1)) +
  theme_minimal()+
  theme(legend.position="bottom", text = element_text(size=20), plot.title = element_text(size = 20, hjust = 0.5))+
  ggtitle("Fishing Effort")

#harvest
p2 <- ggplot(data=patch_population_compare %>% filter(Patch=="Total Harvest",
                                                      FishingEffort == 0.1 | FishingEffort ==0.3 | FishingEffort == 0.7),
             aes(x = as.numeric(Patch2Percent), y = as.numeric(Biomass), color = Scenario)) +
  geom_line(aes(size = as.factor(Scenario))) +
  facet_wrap(~FishingEffort) +
  scale_color_viridis_d(option="turbo", name="", 
                        breaks = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K"),
                        labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K")) +
  scale_size_manual(name = "", values = c(1, 1, 1, 1)) +
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  theme_minimal()+
  theme(legend.position="bottom", text = element_text(size=20), plot.title = element_text(size = 20, hjust = 0.5))

#CPUE
p3 <- ggplot(data=patch_population_compare %>% filter(Patch=="CPUE",
                                                      FishingEffort == 0.1 | FishingEffort ==0.3 | FishingEffort == 0.7),
             aes(x = as.numeric(Patch2Percent), y = as.numeric(Biomass), color = Scenario)) +
  geom_line(aes(size = as.factor(Scenario))) +
  facet_wrap(~FishingEffort) +
  scale_color_viridis_d(option="turbo", name="", 
                        breaks = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K"),
                        labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K")) +
  scale_size_manual(name = "", values = c(1, 1, 1, 1)) +
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  theme_minimal()+
  theme(legend.position="bottom", text = element_text(size=20), plot.title = element_text(size = 20, hjust = 0.5))

ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3,
          nrow=3, common.legend = TRUE, legend="bottom")


#Supplemental Plots with full range of fishing efforts
#Biomass
ggplot(data=patch_population_compare %>% filter(Patch=="Total Biomass"), aes(x = as.numeric(Patch2Percent), 
                                                                             y = as.numeric(Biomass), color = Scenario)) +
  geom_line(aes(size = as.factor(Scenario))) +
  facet_wrap(~FishingEffort) +
  scale_color_viridis_d(option="turbo", name="", 
                        breaks = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K"),
                        labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K")) +
  scale_size_manual(name = "", values = c(1, 1, 1, 1)) +
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  theme_minimal()+
  theme(legend.position="bottom", text = element_text(size=20), plot.title = element_text(size = 20, hjust = 0.5)) +
  ggtitle("Fishing Effort")


#Harvest       

ggplot(data=patch_population_compare %>% filter(Patch=="Total Harvest"), aes(x = as.numeric(Patch2Percent), 
                                                                             y = as.numeric(Biomass), color = Scenario)) +
  geom_line(aes(size = as.factor(Scenario))) +
  facet_wrap(~FishingEffort) +
  scale_color_viridis_d(option="turbo", name="", 
                        breaks = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K"),
                        labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K")) +
  scale_size_manual(name = "", values = c(1, 1, 1, 1)) +
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  theme_minimal()+
  theme(legend.position="bottom", text = element_text(size=20), plot.title = element_text(size = 20, hjust = 0.5)) +
  ggtitle("Fishing Effort")

#CPUE       

ggplot(data=patch_population_compare %>% filter(Patch=="CPUE"), aes(x = as.numeric(Patch2Percent), 
                                                                    y = as.numeric(Biomass), color = Scenario)) +
  geom_line(aes(size = as.factor(Scenario))) +
  facet_wrap(~FishingEffort) +
  scale_color_viridis_d(option="turbo", name="", 
                        breaks = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K"),
                        labels = c("No Temperature", "Temp-dependent r", "Temp-dependent K", "Temp-dependent r + K")) +
  scale_size_manual(name = "", values = c(1, 1, 1, 1)) +
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  theme_minimal()+
  theme(legend.position="bottom", text = element_text(size=20), plot.title = element_text(size = 20, hjust = 0.5)) +
  ggtitle("Fishing Effort")




