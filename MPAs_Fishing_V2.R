library(viridis)
library(ggpubr)
library(grid)
library(purrr)

source('pop_dynamics_fns.R')


timesteps <- 50
##Area - TEST DIFFERENT AREA COMBINATIONS
patch_area_sequences <- list(seq(0, 1, by = 0.1))
patch_area_grid <- do.call(expand.grid, patch_area_sequences)
patch_area_grid$Var2 <- 1 - patch_area_grid$Var1
patch_area_list <- split(patch_area_grid, 1:nrow(patch_area_grid))
number_patches <- ncol(patch_area_grid)

#fishing
fishing_effort_sequences <- list(seq(0, 1, by = 0.1), 0)
fishing_effort_grid <- do.call(expand.grid, fishing_effort_sequences)
fishing_effort_list <- split(fishing_effort_grid, 1:nrow(fishing_effort_grid))

catchability <- 1

#fish
intrinsic_growth_rate <- 0.3 
carrying_capacity <- c(342.4, 342.4) 

###Economic inputs
source('seafood_prices.R')
price <- prices_sum$PricePerG_2021 #price per g in 2021, calculated in seafoodprices.R script 
value_added_ratio <- 0.9  #from Grafeld et al. 2017

##dispersal
S <- 0.7
S_recruit <- 0.1
source('dispersal_fn.R')

##Grid of all different possible combinations
parameter_grid <- expand.grid(patch_area = patch_area_list,
                              fishing_effort = fishing_effort_list)


#outputs
outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_revenue <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))



#Model code
for (iter in 1:nrow(parameter_grid)) {
  
  source('initialize_state_vars.R')
  
  for(t in 2:timesteps){
    for(i in 1:number_patches){
      tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]),
                        number_patches, S)
      tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]),
                        number_patches, S_recruit)
      
      if (as.numeric(parameter_grid[['patch_area']][[iter]])[i] == 0) { 
        population[1, ][i] <- 0
      }
      if (as.numeric(parameter_grid[['patch_area']][[iter]])[i] > 0) {
        population[t, i] <- calculate_population_growth(population[t-1, i], 
                                                        intrinsic_growth_rate, carrying_capacity[i])
      } else {
        population[t, i] <- 0
      }
      recruits[t, i] <- population[t, i] - population[t-1, i]
      adults[t, i] <- population[t, i] - recruits[t, i]
    }
    
  ######DISPERSAL######
   # if (population[t,1] <= carrying_capacity[1] & population[t,2] <= carrying_capacity[2]) { #
    recruits[t,] <- recruits[t,] %*% tmp2
    #}

    
    
    for(i in 1:number_patches){
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      harvest[t, i] <- calculate_fisheries_harvest(adults[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8))
      escapement[t, i] <- calculate_escaped_stock_biomass(adults[t, i], fraction_harvested[t, i])
    }
    
  ######DISPERSAL######
   # if (population[t,1] <= carrying_capacity[1] & population[t,2] <= carrying_capacity[2]) {
    adults[t,] <- escapement[t,] %*% tmp1
    #}
    
    revenue <- (price * harvest) * value_added_ratio
    
    for (i in 1:number_patches) {
      population[t, i] <- adults[t, i] + recruits[t, i]
     # if (population[t, i] > carrying_capacity[i]) {
    #    population[t, i] = carrying_capacity[i]
    #  }
    }
  }
  
  outcome_biomass[iter, ] <- population[t, ]
  outcome_harvest[iter, ] <- harvest[t, ]
  outcome_revenue[iter, ] <- revenue[t, ]
  
}


# cbind outcome with parameter grid
colnames(outcome_biomass) <- c("patch_1_population", "patch_2_population")
colnames(outcome_harvest) <- c("patch_1_harvest", "patch_2_harvest")
colnames(outcome_revenue) <- c("patch_1_revenue", "patch_2_revenue")

outcome <- cbind(parameter_grid, outcome_biomass, outcome_harvest, outcome_revenue)
outcome$fishing_effort_patch1 <- map_dbl(outcome$fishing_effort, 1)
outcome$area_patch1 <- map_dbl(outcome$patch_area, 1)
outcome$area_patch2 <- map_dbl(outcome$patch_area, 2)

outcome$Patch_1_Area_m <- outcome$area_patch1 * 5.04e+8
outcome$Patch_2_Area_m <- outcome$area_patch2 * 5.04e+8

#population x area to get total g in the patch. convert to kg.
outcome$total_population_p1 <- outcome$patch_1_population * outcome$Patch_1_Area_m / 1000
outcome$total_population_p2 <- outcome$patch_2_population * outcome$Patch_2_Area_m / 1000






##Plots
#biomass
outcome$total_population <- outcome$total_population_p1 + outcome$total_population_p2
p1<-ggplot(data=outcome, aes(x=area_patch2*100, y=total_population, col=as.factor(fishing_effort_patch1)))+
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Fish Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

#harvest
p2 <- ggplot(data=outcome, aes(x=area_patch2*100, y=patch_1_harvest/1000, col=as.factor(fishing_effort_patch1)))+
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

##CPUE
outcome$CPUE <- (outcome$patch_1_harvest/1000) / outcome$fishing_effort_patch1
outcome$CPUE[is.na(outcome$CPUE)] <- 0

p3 <- ggplot(data=outcome, aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

#revenue
p4 <- ggplot(data=outcome, aes(x=area_patch2*100, y=patch_1_revenue, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Revenue (USD)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p3 + rremove("xlab"), p2 + rremove("xlab"), p4 + rremove("xlab"),
                  nrow=2, ncol=2, common.legend = TRUE)
annotate_figure(plot, bottom = textGrob("% Coastal Waters Closed to Fishing", gp = gpar(fontsize = 20)))
