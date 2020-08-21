#!/usr/bin/env Rscript
#
#     Set up cases for running model on high and low jellyfish biomass
#
#     Jessica Luo, 2018
#
#
#------------------------------------------------------------------------------
library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("gridExtra")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())
parallel <- TRUE
source("functions.R")


# --------------

# biomass / numeric density
rd <- read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv", stringsAsFactors=FALSE)

# save the originals
rd$nd_orig <- rd$numeric_density
rd$biomass_orig <- rd$Biomass.mgCm3

biomass = biomass_by_biome(rd)
biomass[which(biomass$biome=="TOTAL"),]

# --------------
# Biomass variability
# mean +/- standard deviation / variance
# ctenos: -0.008 +/- 0.604/0.365
# salps: -0.017 +/- 0.586/0.344
# medusae: 0.010 +/- 0.425/0.181
# if you calculate the % using variance then the values are 36.5%, 34.4% and 18.1% for ctenos, salps and medusae

# CASE 0: baseline
# CASE 1: High variance
CASE = '1-high_biomass'
dir.create(str_c('data/',CASE), showWarnings = FALSE)
dir.create(str_c('plots/',CASE), showWarnings = FALSE)
dir.create(str_c('data/',CASE,"/baseline"), showWarnings = FALSE)
dir.create(str_c('plots/',CASE,"/baseline"), showWarnings = FALSE)

rd[rd$rank_phylum=="Ctenophora",c("numeric_density", "Biomass.mgCm3")] = (1 + 0.604) * rd[rd$rank_phylum=="Ctenophora",c("numeric_density", "Biomass.mgCm3")]
rd[rd$rank_phylum=="Chordata",c("numeric_density", "Biomass.mgCm3")] = (1 + 0.586) * rd[rd$rank_phylum=="Chordata",c("numeric_density", "Biomass.mgCm3")]
rd[rd$rank_phylum=="Cnidaria",c("numeric_density", "Biomass.mgCm3")] = (1 + 0.425) * rd[rd$rank_phylum=="Cnidaria",c("numeric_density", "Biomass.mgCm3")]

write.csv(rd, str_c("data/",CASE,"/jelly_biomass_1_deg_grid.csv"), row.names=FALSE)

# CASE 2: Low variance
CASE = '2-low_biomass'
dir.create(str_c('data/',CASE), showWarnings = FALSE)
dir.create(str_c('plots/',CASE), showWarnings = FALSE)
dir.create(str_c('data/',CASE,"/baseline"), showWarnings = FALSE)
dir.create(str_c('plots/',CASE,"/baseline"), showWarnings = FALSE)

rd[rd$rank_phylum=="Ctenophora",c("numeric_density", "Biomass.mgCm3")] = (1 - 0.604) * rd[rd$rank_phylum=="Ctenophora",c("numeric_density", "Biomass.mgCm3")]
rd[rd$rank_phylum=="Chordata",c("numeric_density", "Biomass.mgCm3")] = (1 - 0.586) * rd[rd$rank_phylum=="Chordata",c("numeric_density", "Biomass.mgCm3")]
rd[rd$rank_phylum=="Cnidaria",c("numeric_density", "Biomass.mgCm3")] = (1 - 0.425) * rd[rd$rank_phylum=="Cnidaria",c("numeric_density", "Biomass.mgCm3")]

write.csv(rd, str_c("data/",CASE,"/jelly_biomass_1_deg_grid.csv"), row.names=FALSE)

