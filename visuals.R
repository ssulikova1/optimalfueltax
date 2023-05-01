
rm(list = ls())
library(esquisse)
library(tidyverse)
# uk

colname <- c("Externality", "Values", "Country")

externality <- as.data.frame(c("Fuel pol.", "Dist. pol.", "Congestion",
                             "Accidents in.", "Accidents act.", 
                             "Inactivity"))
#uk
values_uk <- as.data.frame(c(150.1/28, 4.3, 6, 1.9, 1.9, 291.3)) # converting per gallon to per mile. 
#us
values_us <- as.data.frame(c(160.1/24, 5.4, 11.9, 7.6, 6.3, 825))

baseline <- cbind(externality, values_uk)
baseline$country <- "UK"
colnames(baseline) <- colname


baseline2 <- cbind(externality, values_us)
baseline2$country <- "USA"
colnames(baseline2) <- colname

baseline_own <- rbind(baseline, baseline2)
baseline_own$group <- "Own calculations"

## adding parry and small 


inflation_2000 = 1.75 # dollar in 2000 vs dollar in 2022

#uk
values_uk <- as.data.frame(c(6/30, 2, 7, 2.4, NA, NA))
values_uk <- values_uk * inflation_2000
#us
values_us <- as.data.frame(c(6/20, 2, 3.5, 3, NA, NA))
values_us <- values_us* inflation_2000

baseline <- cbind(externality, values_uk)
baseline$country <- "UK"
colnames(baseline) <- colname


baseline2 <- cbind(externality, values_us)
baseline2$country <- "USA"
colnames(baseline2) <- colname

baseline_p <- rbind(baseline, baseline2)
baseline_p$group <- "Parry & Small (2005)"

baselines <- rbind(baseline_own, baseline_p)
baselines_temp <- baselines
baselines_temp$Values[baselines_temp$Externality=="Inactivity"] <- 0

esquisser(baselines)




library(dplyr)
library(ggplot2)

## overall plot

baselines %>%
#  filter(!(group %in% "Parry & Small (2005)")) %>%
  ggplot() +
  aes(x = Externality, weight = Values, fill = Country) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(UK = "#6DABF8",
               USA = "#1018E2")) +
 # scale_fill_hue(direction = 1) +
  labs(y = "Cost, US cents/mile") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(vars(group),  ncol = 1L, 
             nrow = 2L, scales = "free_y")


## excluding physical inactivity


baselines_temp %>%
  #  filter(!(group %in% "Parry & Small (2005)")) %>%
  ggplot() +
  aes(x = Externality, weight = Values, fill = Country) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(UK = "#6DABF8",
               USA = "#1018E2")) +
  # scale_fill_hue(direction = 1) +
  labs(y = "Cost, US cents/mile") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(vars(group),  ncol = 1L, 
             nrow = 2L)


baselines_temp %>%
    filter((group %in% "Parry & Small (2005)")) %>%
  ggplot() +
  aes(x = Externality, weight = Values, fill = Country) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(UK = "#6DABF8",
               USA = "#1018E2")) +
  # scale_fill_hue(direction = 1) +
  labs(y = "Cost, US cents/mile") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 


baselines_temp %>%
  filter(!(group %in% "Parry & Small (2005)")) %>%
  ggplot() +
  aes(x = Externality, weight = Values, fill = Country) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(UK = "#6DABF8",
               USA = "#1018E2")) +
  # scale_fill_hue(direction = 1) +
  labs(y = "Cost, US cents/mile") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 

baselines %>%
  filter(!(group %in% "Parry & Small (2005)")) %>%
  ggplot() +
  aes(x = Externality, weight = Values, fill = Country) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(UK = "#6DABF8",
               USA = "#1018E2")) +
  # scale_fill_hue(direction = 1) +
  labs(y = "Cost, US cents/mile") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 
