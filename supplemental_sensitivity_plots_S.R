library(grid)
library(ggpubr)

#run model
source('supplemental_sensitivity_s.R')


outcome$p1_area_m <- outcome$area_patch1 * 5.04e+8
outcome$p2_area_m <- outcome$area_patch2 * 5.04e+8

outcome$total_pop_p1 <- outcome$patch_1_population * outcome$p1_area_m
outcome$total_pop_p2 <- outcome$patch_2_population * outcome$p2_area_m
outcome$total_pop <- outcome$total_pop_p1 + outcome$total_pop_p2
outcome$CPUE <- (outcome$patch_1_harvest/1000) / outcome$fishing_effort_patch1
outcome$CPUE[is.na(outcome$CPUE)] <- 0

##S
p1<-ggplot(data=outcome%>%filter(round(S_recruit,1)==0.1,
                                 round(intrinsic_growth_rate, 1) == 0.3,
                                 carrying_capacity==350), 
           aes(x=area_patch2*100, y=total_pop, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~S) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

p2<-ggplot(data=outcome%>%filter(round(S_recruit,1)==0.1,
                                 round(intrinsic_growth_rate, 1) == 0.3,
                                 carrying_capacity==350), 
           aes(x=area_patch2*100, y=patch_1_harvest, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~S) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

p3<-ggplot(data=outcome%>%filter(round(S_recruit,1)==0.1,
                                 round(intrinsic_growth_rate, 1) == 0.3,
                                 carrying_capacity==350), 
           aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~S) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3, nrow=3, common.legend = TRUE, legend = "right")
annotate_figure(plot, bottom = text_grob("S (recruit) = 0.1, r = 0.3, K = 350.", hjust = 1, x=1, size = 10),
                top = textGrob("Site Fidelity (adult)", gp = gpar(fontsize = 20)))


##S recruit

p1<-ggplot(data=outcome%>%filter(round(S,1)==0.7,
                                 round(intrinsic_growth_rate, 1) == 0.3,
                                 carrying_capacity==350), 
           aes(x=area_patch2*100, y=total_pop, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~S_recruit) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Biomass (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

p2<-ggplot(data=outcome%>%filter(round(S,1)==0.7,
                                 round(intrinsic_growth_rate, 1) == 0.3,
                                 carrying_capacity==350), 
           aes(x=area_patch2*100, y=patch_1_harvest, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~S_recruit) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="Total Harvest (kg)") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

p3<-ggplot(data=outcome%>%filter(round(S,1)==0.7,
                                 round(intrinsic_growth_rate, 1) == 0.3,
                                 carrying_capacity==350), 
           aes(x=area_patch2*100, y=CPUE, col=as.factor(fishing_effort_patch1))) +
  geom_line(aes(size = as.factor(fishing_effort_patch1))) +
  facet_wrap(~S_recruit) +
  theme_minimal() + 
  labs(x="% Coastal Waters Closed to Fishing", y="CPUE") +
  scale_color_viridis_d(name="Fishing Effort") +
  scale_size_manual(name = "Fishing Effort", values = c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5)) +
  theme(legend.position="top", text = element_text(size=20))

plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3,
                  nrow=3, common.legend = TRUE, legend="right")
annotate_figure(plot, bottom = text_grob("S (adult) = 0.7, r = 0.3, K = 350.", hjust = 1, x=1, size = 10),
                top = textGrob("Site Fidelity (recruit)", gp = gpar(fontsize = 20)))



