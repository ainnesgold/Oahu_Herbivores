#for some reason it throws an error if you press source to run it all, but if you run all the parameter stuff first,
#and then run the loop, it works fine!

#sensitivity analysis on temp dependent r and k models
library(purrr)
library(tidyverse)
library(ggpubr)
library(grid)

calculate_population_growth <- function(population, intrinsic_growth_rate, carrying_capacity) {
  population <-population * exp(intrinsic_growth_rate * (1 - population / carrying_capacity))
  return(population)
}

#varying the constants in temp dependent r equation
a = seq(0.1, 0.5, by = 0.1)
b = 0
c = seq(-0.01, 0, by = 0.0025)

calculate_r_temp <- function(SST_dev, a, b, c) {
  r_temp <- a + b*SST_dev + c*SST_dev^2
  return(r_temp)
}

calculate_population_growth_temp_r <- function(population, r_temp, carrying_capacity) {
  population <-population * exp(r_temp * (1 - population / carrying_capacity))
  return(population)
}
#varying the constants in temp dependent K equation
m <- seq(-14.7, -6.7, by = 2) #-10.76933 
d <- seq(535, 735, by = 50) #635.5131

calculate_CC_temp <- function(m, mean_temp_series, d) {
  CC_temp <- (m * mean_temp_series) + d
  return(CC_temp)
}

calculate_population_growth_temp_CC <- function(population, intrinsic_growth_rate, CC_temp) {
  population <-population * exp(intrinsic_growth_rate * (1 - population / CC_temp))
  return(population)
}

calculate_population_growth_temp_r_CC <- function(population, r_temp, CC_temp) {
  population <-population * exp(r_temp * (1 - population / CC_temp))
  return(population)
}

calculate_fraction_harvested <- function(fishing_effort, catchability) {
  fraction_harvested <- 1 - exp(-fishing_effort * catchability)
  return(fraction_harvested)
}

calculate_fisheries_harvest <- function(population, fraction_harvested, patch_area_m2) {
  harvest <- population * fraction_harvested * patch_area_m2
  return(harvest)
}

calculate_escaped_stock_biomass <- function(population, fraction_harvested) {
  escaped_stock_biomass <- population * (1 - fraction_harvested)
  return(escaped_stock_biomass)
}

#parameters
#CLIMATE 
source('SST_Projections.R')
SST_dev <- projections[['anomaly']]
mean_temps_series <- mean_temps_crop[['extracted_mean_temps']]

#varying the constants in temp dependent r equation
a = seq(0.1, 0.5, by = 0.1)
b = 0
c = seq(-0.01, 0, by = 0.0025)

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
carrying_capacity <- c(342, 342)

source('seafood_prices.R')
price <- prices_sum$PricePerG_2021 #price per g in 2021, calculated in seafoodprices.R script 
value_added_ratio <- 0.9  #from Grafeld et al. 2017

S <- 0.7
S_recruit <- 0.1
source('dispersal_fn.R')

##Grid of all different possible combinations
parameter_grid <- expand.grid(a = a, c = c, m = m, d = d,
                              patch_area = patch_area_list,
                              fishing_effort = fishing_effort_list)


outcome_biomass <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))


for (iter in 1:nrow(parameter_grid)) {
  
  tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  tmp2 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S_recruit)
  
  source('initialize_state_vars.R')
  
  ## Run Population Model ---------------------------------------------------------------------------------------------------------
  for(t in 2:timesteps){ # start at 2 because we set the initial starting value which is t = 1
    for(i in 1:number_patches){
    
      r_temp[t, i] <- calculate_r_temp(SST_dev[t], parameter_grid[['a']][[iter]], b, parameter_grid[['c']][[iter]])
      CC_temp[t, i] <- calculate_CC_temp(parameter_grid[['m']][[iter]], mean_temps_series[t], parameter_grid[['d']][[iter]])
      population[t, i] <- calculate_population_growth_temp_r_CC(population[t-1, i], r_temp[t, i], 
                                                                CC_temp[t, i])
      
      recruits[t, i] <- population[t, i] - population[t-1, i]
    }
    recruits_dispersal[t,] <- recruits[t,] %*% tmp2
    
    for(i in 1:number_patches){
      population[t, i] <- population[t, i] - recruits[t, i] + recruits_dispersal[t, i]
      fraction_harvested[t, i] <- calculate_fraction_harvested(as.numeric(parameter_grid[['fishing_effort']][[iter]])[i],
                                                               catchability)
      harvest[t, i] <- calculate_fisheries_harvest(population[t, i], fraction_harvested[t, i], 
                                                   (as.numeric(parameter_grid[['patch_area']][[iter]])[i]*5.04e+8)) 
      escapement[t, i] <- calculate_escaped_stock_biomass(population[t, i], fraction_harvested[t, i])
    }  
    population <- escapement[,] %*% tmp1
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

outcome$p1_area_m <- outcome$area_patch1 * 5.04e+8
outcome$p2_area_m <- outcome$area_patch2 * 5.04e+8

outcome$total_pop_p1 <- outcome$patch_1_population * outcome$p1_area_m
outcome$total_pop_p2 <- outcome$patch_2_population * outcome$p2_area_m
outcome$total_pop <- outcome$total_pop_p1 + outcome$total_pop_p2
outcome$CPUE <- (outcome$patch_1_harvest/1000) / outcome$fishing_effort_patch1
outcome$CPUE[is.na(outcome$CPUE)] <- 0




##plots
#names of constants are different in the code than they are in the titles/writing because I changed them when writing for clarity.
#it is confusing. sorry. the title of the figure matches how the constants are named in the paper, but the coded variables don't.

#first constant: r (coded as "a")
p1 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005, round(m,1)==-10.7,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=total_pop, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~a) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p2 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005, round(m,1)==-10.7,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=patch_1_harvest, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~a) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p3 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005, round(m,1)==-10.7,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~a) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3, nrow=3, common.legend = TRUE, legend = "right")
annotate_figure(plot, bottom = text_grob("a = -0.005, m = -10.7, b = 635.", hjust = 1, x=1, size = 10),
                top = textGrob("r", gp = gpar(fontsize = 20)))


#second constant: a (coded as "c")
p1 <- ggplot(data=outcome%>%filter(round(a,1)==0.3, round(m,1)==-10.7,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=total_pop, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~c) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p2 <- ggplot(data=outcome%>%filter(round(a,1)==0.3, round(m,1)==-10.7,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=patch_1_harvest, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~c) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p3 <- ggplot(data=outcome%>%filter(round(a,1)==0.3, round(m,1)==-10.7,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~c) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3, nrow=3, common.legend = TRUE, legend = "right")
annotate_figure(plot, bottom = text_grob("r = 0.3, m = -10.7, b = 635.", hjust = 1, x=1, size = 10),
                top = textGrob("a", gp = gpar(fontsize = 20)))


#third constant: m (coded as "m")
p1 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005, round(a, 1)==0.3,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=total_pop, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~m) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p2 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005,  round(a, 1)==0.3,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=patch_1_harvest, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~m) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p3 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005,  round(a, 1)==0.3,
                                   round(d, 1)==635), 
             aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~m) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3, nrow=3, common.legend = TRUE, legend = "right")
annotate_figure(plot, bottom = text_grob("r = 0.3, a = -0.005, b = 635.", hjust = 1, x=1, size = 10),
                top = textGrob("m", gp = gpar(fontsize = 20)))



#fourth constant: b (coded as "d")
p1 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005, round(a, 1)==0.3,
                                   round(m, 1)==-10.7), 
             aes(x=area_patch2*100, y=total_pop, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~d) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p2 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005,  round(a, 1)==0.3,
                                   round(m, 1)==-10.7), 
             aes(x=area_patch2*100, y=patch_1_harvest, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~d) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

p3 <- ggplot(data=outcome%>%filter(round(c,3)==-0.005,  round(a, 1)==0.3,
                                   round(m, 1)==-10.7), 
             aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~d) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="bottom", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3, nrow=3, common.legend = TRUE, legend = "right")
annotate_figure(plot, bottom = text_grob("r = 0.3, a = -0.005, m = -10.7", hjust = 1, x=1, size = 10),
                top = textGrob("b", gp = gpar(fontsize = 20)))







