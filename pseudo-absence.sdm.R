setwd("~/Library/CloudStorage/OneDrive-NTNU/NARM/MSNARM/Footprint/Data/R.master/leirelv")

library(tidyverse)
library(MASS)
library(mgcv)
library(sf)
library(performance)
library(readxl)
library(ggplot2)

#baseline data moose 
moose_baseline <- st_read("moose4.shp")
#baseline data roedeer
roedeer_baseline <- st_read("roedeer4.shp")
#baseline data redfox 
redfox_baseline <- st_read("redfox4.shp")


#---------------- script for fixing pseudo absence models -----------------------------

#---- MOOSE:

moose_baseline$presence <- as.integer(moose_baseline$obs_count > 0)
moose_baseline_presence <- subset(moose_baseline, presence == 1)
#View(moose_baseline_presence)

moose_baseline_absence <- subset(moose_baseline, presence == 0)
#View(moose_baseline_absence)

#randomly sample absences
set.seed(123)
abs_moose_sample <- moose_baseline_absence[sample(nrow(moose_baseline_absence), size = 10*nrow(moose_baseline_presence), replace = FALSE),]

#View(abs_moose_sample)

#merge absence and presence data
moose_sdm <- rbind(moose_baseline_presence, abs_moose_sample)

#weighting pseudo absences
np_moose <- sum(moose_sdm$presence == 1)
na_moose <- sum(moose_sdm$presence == 0)

moose_sdm$w <- ifelse(
  moose_sdm$presence == 1,
  1,
  np_moose/na_moose
)

#fit weighted GLM
model_moose_sdm <- glm(
  presence ~ building_c + road_area  + areatype20 + areatype30 + areatype50 + areatype60 + areatype81,
  data = moose_sdm,
  family = binomial, weights = moose_sdm$w
)

summary(model_moose_sdm)

#predict on the full grid
moose_baseline$p_true <- predict(model_moose_sdm, newdata = moose_baseline, type = "response")

summary(moose_baseline$p_true)

#View(moose_baseline)

#------ ROEDEER:
roedeer_baseline$presence <- as.integer(roedeer_baseline$obs_count > 0)
roedeer_baseline_presence <- subset(roedeer_baseline, presence == 1)
#View(roedeer_baseline_presence)

roedeer_baseline_absence <- subset(roedeer_baseline, presence == 0)
#View(roedeer_baseline_absence)

#randomly sample absences
set.seed(123)
abs_roedeer_sample <- roedeer_baseline_absence[sample(nrow(roedeer_baseline_absence), size = 10*nrow(roedeer_baseline_presence), replace = FALSE),]

#View(abs_roedeer_sample)

#merge absence and presence data
roedeer_sdm <- rbind(roedeer_baseline_presence, abs_roedeer_sample)

#weighting pseudo absences
np_roedeer <- sum(roedeer_sdm$presence == 1)
na_roedeer <- sum(roedeer_sdm$presence == 0)

roedeer_sdm$w <- ifelse(
  roedeer_sdm$presence == 1,
  1,
  np_roedeer/na_roedeer
)

#fit weighted GLM
model_roedeer_sdm <- glm(
  presence ~ building_c + road_area  + areatype20 + areatype30 + areatype50 + areatype60 + areatype81,
  data = roedeer_sdm,
  family = binomial, weights = roedeer_sdm$w
)

summary(model_roedeer_sdm)

#predict on the full grid
roedeer_baseline$p_true <- predict(model_roedeer_sdm, newdata = roedeer_baseline, type = "response")

summary(roedeer_baseline$p_true)

#View(roedeer_baseline)

#------ REDFOX:
redfox_baseline$presence <- as.integer(redfox_baseline$obs_count > 0)
redfox_baseline_presence <- subset(redfox_baseline, presence == 1)
#View(redfox_baseline_presence)

redfox_baseline_absence <- subset(redfox_baseline, presence == 0)
#View(redfox_baseline_absence)

#randomly sample absences
set.seed(123)
abs_redfox_sample <- redfox_baseline_absence[sample(nrow(redfox_baseline_absence), size = 10*nrow(redfox_baseline_presence), replace = FALSE),]

#View(abs_redfox_sample)

#merge absence and presence data
redfox_sdm <- rbind(redfox_baseline_presence, abs_redfox_sample)

#weighting pseudo absences
np_redfox <- sum(redfox_sdm$presence == 1)
na_redfox <- sum(redfox_sdm$presence == 0)

redfox_sdm$w <- ifelse(
  redfox_sdm$presence == 1,
  1,
  np_redfox/na_redfox
)

#fit weighted GLM
model_redfox_sdm <- glm(
  presence ~ building_c + road_area  + areatype20 + areatype30 + areatype50 + areatype60 + areatype81,
  data = redfox_sdm,
  family = binomial, weights = redfox_sdm$w
)

summary(model_redfox_sdm)

#predict on the full grid
redfox_baseline$p_true <- predict(model_redfox_sdm, newdata = redfox_baseline, type = "response")

summary(redfox_baseline$p_true)

# --------------- species baseline SDM models ---------------

#MOOSE
model_moose_sdm <- glm(
  presence ~ building_c + road_area  + areatype20 + areatype30 + areatype50 + areatype60 + areatype81,
  data = moose_sdm,
  family = binomial, weights = moose_sdm$w
)
summary(model_moose_sdm)

#ROEDEER
model_roedeer_sdm <- glm(
  presence ~ building_c + road_area  + areatype20 + areatype30 + areatype50 + areatype60 + areatype81,
  data = roedeer_sdm,
  family = binomial, weights = roedeer_sdm$w
)
summary(model_roedeer_sdm)

#REDFOX
model_redfox_sdm <- glm(
  presence ~ building_c + road_area + areatype20 + areatype30  + areatype50 + areatype60 + areatype81,
  data = redfox_sdm,
  family = binomial, weights = redfox_sdm$w
)
summary(model_redfox_sdm)

sum(moose_sdm$presence)
sum(roedeer_sdm$presence)
sum(redfox_sdm$presence)

# --------------- species predictions-------------------------
#MOOSE
moose_baseline$p_true <- predict(model_moose_sdm, newdata = moose_baseline, type = "response")
summary(moose_baseline$p_true)

#ROEDEER
roedeer_baseline$p_true <- predict(model_roedeer_sdm, newdata = roedeer_baseline, type = "response")
summary(roedeer_baseline$p_true)

#REDFOX
redfox_baseline$p_true <- predict(model_redfox_sdm, newdata = redfox_baseline, type = "response")
summary(redfox_baseline$p_true)

