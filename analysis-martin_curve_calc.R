#!/usr/bin/env Rscript

library("plyr")
library("reshape2")
library("ggplot2")
library("gridExtra")
library("stringr")
library("viridis")
library("RColorBrewer")

source("functions.R")
Spectral <- brewer.pal(8, "Spectral")

# load biomes
load('data/biomes_chl_mixedlayer.Rdata')
N_BIOME=4

seafloor <- read.csv("data/woa_seafloor_depth.csv")


# Import Export to depth data for all three biomes
exportres = data.frame()
uores = data.frame()

for (C_ID in 1:3){
  CASES = c("0-baseline","1-high_biomass","2-low_biomass")
  CASE = CASES[C_ID]
  
  load(file=str_c("data/",CASE,"/baseline/allexport.Rdata")) # units are all g C m^-3 y^-1
  load(file=str_c("data/",CASE,"/baseline/upperocean_ens_mean.Rdata"))
  
  allexport$case <- CASE
  ens_mean$case <- CASE
  
  exportres <- rbind(exportres, allexport)
  uores <- rbind(uores, ens_mean)
}

exportres$case = factor(exportres$case, levels=c("2-low_biomass", "0-baseline", "1-high_biomass"))
exportres <- join(exportres, biomes_df, by=c("lat","lon"))

sumres=exportres[exportres$taxon=="sum",]
exprt=sumres[,c("lat",'lon','export','case','biome')]
seq_d=sumres[,c("lat",'lon','seq_d','case','biome')]
seaflr=sumres[,c("lat",'lon','seafloor','case','biome')]

exprt <- join(exprt, seafloor, by=c("lat","lon"))
exprt$depth = ifelse(exprt$depth > 100, 100, exprt$depth)
seq_d <- join(seq_d, seafloor, by=c("lat","lon"))
seq_d$depth = ifelse(seq_d$depth > 1000, 1000, seq_d$depth)
seaflr = join(seaflr, seafloor, by=c("lat","lon"))
names(seaflr) = names(seq_d) = names(exprt)

db=rbind(exprt, seq_d, seaflr)
db$log10depth=log10(db$depth)
db$log10export=log10(db$export)
db = db[is.finite(db$log10export) & is.finite(db$log10depth), ]

m = lm(-log10depth~log10export, data=db)
print('Results from all sinking POC:')
print(summary(m))
# b = 0.18 +/- 0.003


###
exprt=sumres[,c("lat",'lon','export_cf','case','biome')]
seq_d=sumres[,c("lat",'lon','seq_d_cf','case','biome')]
seaflr=sumres[,c("lat",'lon','seafloor_cf','case','biome')]

exprt <- join(exprt, seafloor, by=c("lat","lon"))
exprt$depth = ifelse(exprt$depth > 100, 100, exprt$depth)
seq_d <- join(seq_d, seafloor, by=c("lat","lon"))
seq_d$depth = ifelse(seq_d$depth > 1000, 1000, seq_d$depth)
seaflr = join(seaflr, seafloor, by=c("lat","lon"))
names(seaflr) = names(seq_d) = names(exprt)

db=rbind(exprt, seq_d, seaflr)
db$log10depth=log10(db$depth)
db$log10export=log10(db$export)
db = db[is.finite(db$log10export) & is.finite(db$log10depth), ]

m = lm(-log10depth~log10export, data=db)
print('Results from carcass flux only:')
print(summary(m))
# b = 0.15 +/- 0.003


###
exprt=sumres[,c("lat",'lon','export_eg','case','biome')]
seq_d=sumres[,c("lat",'lon','seq_d_eg','case','biome')]
seaflr=sumres[,c("lat",'lon','seafloor_eg','case','biome')]

exprt <- join(exprt, seafloor, by=c("lat","lon"))
exprt$depth = ifelse(exprt$depth > 100, 100, exprt$depth)
seq_d <- join(seq_d, seafloor, by=c("lat","lon"))
seq_d$depth = ifelse(seq_d$depth > 1000, 1000, seq_d$depth)
seaflr = join(seaflr, seafloor, by=c("lat","lon"))
names(seaflr) = names(seq_d) = names(exprt)

db=rbind(exprt, seq_d, seaflr)
db$log10depth=log10(db$depth)
db$log10export=log10(db$export)
db = db[is.finite(db$log10export) & is.finite(db$log10depth), ]

m = lm(-log10depth~log10export, data=db)
print('Results from fecal matter only:')
print(summary(m))

# b = 0.12 +/- 0.002
