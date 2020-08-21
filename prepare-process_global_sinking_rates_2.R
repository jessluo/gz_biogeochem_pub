#
#         Calculate the export flux for GZ
#             - use results from ML model for biomass and egestion
#
#         
#--------------------------------------------------------------------
library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())
parallel <- FALSE

source("functions.R")


## { Load files ---------------------------------------------------------
# import seafloor depth
seafloor <- read.csv("data/woa_seafloor_depth.csv")
seafloor$latlon <- paste0(seafloor$lat, ", ", seafloor$lon)

# import export depth
export.depth <- read.csv("data/woa_export_depth.csv")

# import sequestration depth
seq.depth <- read.csv("data/woa_sequestration_depth.csv")

rd = read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv", as.is=TRUE)
rdc = dcast(rd, lat+lon~rank_phylum, value.var="Biomass.mgCm3")

# mimic export flux dataset
cf <- cbind(rdc, Cnidaria_eg=rdc$Cnidaria, Ctenophora_eg=rdc$Ctenophora, Chordata_eg=rdc$Chordata)
cf_matrix <- as.matrix(cf[,3:8])

# create a template df for the biological data
cfcpy <- cf_matrix
cfcpy[which(!is.na(cfcpy))] <- 9999

cfSAVE <- cf # save a copy
cf[,3:8] <- cfcpy # matrix of NAs and 9999's into cf

# }

## { Set up the starting depths - DO IT ONCE AND SAVE OUTPUT --------------------------------------------

load(file="data/global_particle_sinking.Rdata")

#
bio_depth <- NULL
# LOOKS like the depths are cnid=20, cten=20, chor=60
bio_depth$taxon <- c("Cnidaria", "Ctenophora", "Chordata")
bio_depth$woa_depths <- c(20, 20, 50)
bio_depth <- as.data.frame(bio_depth, stringsAsFactors=FALSE)

# select out the locations where the seafloor depths are 60 m and shallower
shallow <- seafloor[seafloor$depth <= 60 & seafloor$depth > 30,]
very_shallow <- seafloor[seafloor$depth <= 30,]

shallow <- join(shallow[,c("lon","lat","latlon")], dm[,c("depth", "temp", "C", "time", "m_ratio", "latlon")], by="latlon")
shallow <- shallow[!is.na(shallow$depth),]
# tested this and this way of joining is exactly the same as joining by c("lat", "lon")
very_shallow <- join(very_shallow[,c("lon","lat","latlon")], dm[,c("depth", "temp", "C", "time", "m_ratio", "latlon")], by="latlon")

dm <- dm[which(dm$latlon %ni% c(shallow$latlon, very_shallow$latlon)),]


# match up the global particle sinking values with the depths at which the jelly falls will start
# for the non-shallow locations
for (i in 1:3) {
  cf$depth <- bio_depth[i, "woa_depths"]
  idy <- which(names(cf) %in% c("lat", "lon", "latlon", bio_depth$taxon[i], str_c(bio_depth$taxon[i], "_eg"), "depth"))
  subset <- cf[,idy]
  dm <- join(dm, subset, by=c("lat","lon", "depth"), type="left")
}

# for the shallow locations
for (i in 1:3){
  cf$depth <- bio_depth[i, "woa_depths"] / 2
  idy <- which(names(cf) %in% c("lat", "lon", "latlon", bio_depth$taxon[i], str_c(bio_depth$taxon[i], "_eg"), "depth"))
  subset <- cf[,idy]
  shallow <- join(shallow, subset, by=c("lat","lon", "depth"), type="left")
}

cf$depth <- 0
very_shallow <- join(very_shallow, cf, by=c("lat","lon", "depth"), type="left")

# combine
dm <- rbind(dm, shallow, very_shallow)

# get rid of the locations where there is no jelly biomass at all
dm <- join(cf[,c("lat", "lon")], dm, by=c("lat", "lon"), type="left", match="all")

dm <- dm[,which(names(dm) != "latlon")]


# save as some kind of template
head(dm)
# lat   lon depth    temp   C time m_ratio Cnidaria Cnidaria_cf Ctenophora Ctenophora_cf Chordata Chordata_cf
# 1 -62.5 -57.5     0 -0.9266 100   NA  1.0000       NA          NA         NA            NA       NA          NA
# 2 -62.5 -57.5     5 -0.9925 100 0.05  0.9939       NA          NA         NA            NA       NA          NA
# 3 -62.5 -57.5    10 -0.9988 100 0.05  0.9879       NA          NA         NA            NA       NA          NA
# 4 -62.5 -57.5    15 -0.9723 100 0.05  0.9819       NA          NA         NA            NA       NA          NA
# 5 -62.5 -57.5    20 -0.9989 100 0.05  0.9759     9999        9999         NA            NA       NA          NA
# 6 -62.5 -57.5    25 -1.0150 100 0.05  0.9700       NA          NA         NA            NA       NA          NA
save(dm, file="data/dm_particle_sinking_matched_template.Rdata")

# }