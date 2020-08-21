#!/usr/bin/env Rscript
#
#         Preparing observational GZ data for modeling
# 
#
# --------------------------------------------------------------------------------------------


library("plyr")
library("reshape2")
library("stringr")

LATLON <- c("lat", "lon")
source("functions.R")

NEW_DOWNLOAD=FALSE
dir.create("data")
dir.create("data/0-baseline")

## { Import and process raw jellyfish biomass data --------------------------------------------

# read in ungridded biomass data (raw data)
rd <- read.csv("raw_data/rawdata.csv", stringsAsFactors=F)
rd$numeric_density <- as.numeric(rd$numeric_density)

# exclude Appendiculariansdownload.file(filename, loc_temp)
rd <- rd[which(rd$rank_class != "Appendicularia"),]

# read in original JEDI database
if(NEW_DOWNLOAD){
	url = "https://erddap.bco-dmo.org/erddap/tabledap/bcodmo_dataset_526852.csv"
	loc = "raw_data/JeDI_BCODMO_output.csv"
	download.file(url, loc)
	
	# if this fails due to SSL verification, just download it manually and save it to the designated spot
}

d <- read.csv("raw_data/JeDI_BCODMO_output.csv", stringsAsFactors=FALSE)
# units are in row 1
units = d[1,]
d <- d[-1,]

d$density <- as.numeric(d$density)
d$actual_count <- as.numeric(d$actual_count)
d$depth_integrated_density <- as.numeric(d$depth_integrated_density)
d$latitude <- as.numeric(d$latitude)
d$longitude <- as.numeric(d$longitude)

# plot data types from JeDI
d$lat=floor(d$latitude) + 0.5
d$lon=floor(d$longitude) + 0.5
dd = ddply(d[d$rank_phylum %in% c("Cnidaria", "Ctenophora", "Chordata") & d$taxon != "larvacean",], 
           ~project_title+owner_dataset+location_name+rank_phylum+data_type+lat+lon, function(x){
             res=ddply(x, ~year, function(xi){return(data.frame(numeric_density=mean(xi$density)))})
             n=nrow(res)
             return(data.frame(length=n, numeric_density=mean(x$density)))
           }, .progress='text')

# load coastline
coastline.world <- read.csv("data/gshhg_world_c.csv")
# download coastline data from https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/ and convert to csv

p = ggplot(mapping=aes(x=lon, y=lat)) + 
  geom_point(aes(color=data_type, size=length), data=dd, alpha=0.5) + 
  scale_x_continuous("", breaks=c(-90, 0, 90), expand=c(0,0)) +
  scale_y_continuous("", breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + 
  scale_size_continuous("# years in series", range=c(0.15,4), breaks=c(1,10,25,60)) +
  labs(color="Data Type") + 
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                     panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  facet_grid(rank_phylum~.)

ggsave("plots/JEDI_data_types.pdf", p, width=7.5, height=9)


## { Remove and add data ----------------------------------------------------------------------
# 1) Clean up unknown taxa
# assign all those records that do not have a rank_phylum name into Cnidaria
# only 377 entires, reasonable since they are marked as jellies/medusa
rd[which(rd$rank_phylum %in% c("None", "")), "rank_phylum"] <- "Cnidaria" 

# remove anomalously high data from the central north pacific
anom <- d[d$latitude > 12 & d$latitude < 24 & d$longitude > -164 & d$longitude < -147,]
anom <- anom[which(anom$depth_integrated_density > 2000),]
# --> these are all values from the COPEPOD-BCF-POFI database, from 1959. Seems highly suspect. definitely outliers

# # join and pull out the index numbers that match the really high numbers
match_by <- c("latitude", "longitude", "taxon", "rank_phylum", "rank_class", "rank_order", "numeric_density")
anom <- join(cbind(rd, id=seq_len(nrow(rd))), anom, by=match_by, type = "inner")
rd <- rd[which(as.numeric(row.names(rd)) %ni% anom$id),]

# check
NPTG <- rd[rd$Long.Prov=="NPTG",]
NPTG <- arrange(NPTG, desc(numeric_density))
head(NPTG)
# --> there are still three records that are very high of siphonophores, which do look suspect because it seems to match the kinds of records that we excluded. Will go ahead and exclude them
rd <- rd[which(!(rd$Long.Prov=="NPTG" & rd$Biomass.mgCm3 > 2000)),]

# round the data to the nearest 1-degree grid cell
rd$lat <- floor(rd$latitude) + 0.5
rd$lon <- floor(rd$longitude) + 0.5
rd$lon <- ifelse(rd$lon == 180.5, 179.5, rd$lon)

# pull together
# one row per phyla per lat/lon
rd <- ddply(rd, rank_phylum~lat+lon, function(x){
  latitude <- mean(x$latitude)
  longitude <- mean(x$longitude)
  numeric_density <- gm_mean(x$numeric_density, na.rm=TRUE)
  Biomass.mgCm3 <- gm_mean(x$Biomass.mgCm3, na.rm=TRUE)
  n <- nrow(x)
  return(data.frame(latitude, longitude, numeric_density, Biomass.mgCm3, n))
}, .progress="text")


# 2) ADD GOM DATA FROM KELLY ROBINSON
gom <- read.csv("raw_data/GoMex_25yr mean LargeJellyfish Biomass_AL-LA_1985-2009.csv")
names(gom) <- c("longitude", "latitude", "WW", "n")

# gom trawls are bottom trawls. grab bottom depth from bathymetry data
# download from: http://gcoos.tamu.edu/products/topography/SRTM30PLUS.html
gom_bathy <- read.table("raw_data/gom_srtm30_plus.txt", header=TRUE)

gom$z <- apply(gom, 1, function(x){
  idx <- which_closest_value(x[1], gom_bathy$long)
  idy <- which_closest_value(x[2], gom_bathy$lat[idx])
  z <- gom_bathy[idx,][idy,"z"]
  return(z)
})

gom <- cbind(gom, data.frame(rank_phylum="Cnidaria"))

# WW in kg, convert to g
gom$WW <- gom$WW * 1000

# GOM jellies mix of Aurelia and Chrysaora
# - from Kylie Pitt's C/N/ESD spreadsheet
# Aurelia, 0.202 g WW = 0.000201 g C, so 0.000995 g C / g WW
# Chrysaora, 880 g WW = 2.47 g C, so 0.0028 g C / g WW
# assume 50/50 mix??? 
conversion <- mean(c(0.0009950495, 0.002806818))
# -> 0.0019 g C / g WW
gom$Biomass.mgCm3 <- gom$WW * conversion * 1000 * abs(gom$z)

# grid into 1 degree grids
gom$lon <- floor(gom$longitude) + 0.5
gom$lat <- floor(gom$latitude) + 0.5

# one row per phyla per lat/lon
gom <- ddply(gom, rank_phylum~lat+lon, function(x){
  latitude <- mean(x$latitude)
  longitude <- mean(x$longitude)
  Biomass.mgCm3 <- gm_mean(x$Biomass.mgCm3)
  n <- sum(x$n)
  return(data.frame(latitude, longitude, Biomass.mgCm3, n))
})

# on average, one Aurelia aurita is ~0.2 mg C
# average of juv and adult chrysaora is 1008.5 mg
# so the two averaged together is 504.35 mg C
gom$numeric_density <- gom$Biomass.mgCm3 / 504.35

# 3) ADD NCC DATA FROM RIC BRODEUR
# CHRYSAORA FUSCESCENS - units are kg WW / 1000 m^3
ncc1 <- read.csv("raw_data/BPAcatch_4kelly_species_6.csv")
ncc2 <- read.csv("raw_data/BPAcatch_4kelly_species_9.csv")
col <- c("YEAR", "MONTH", "LONGITUDE", "LATITUDE", "seanett")
ncc <- rbind(ncc1[,col], ncc2[,col])

# convert to carbon biomass (g)
ncc$Biomass <- ncc$seanett * 0.00280 # Shenker 1985, C. fuscescens carbon as 0.280% of wet weight.
ncc$Biomass <- ncc$Biomass * 1000 # convert to mg C
# brodeur et al. 2014: effective mouth area of 123 m^2, and each trawl was towed over the upper 20 m at ~6 km/hr for 30 min. ==> vol per tow = 369000 m^3

# Shenker 1985 studied Chrysaora off the NCC in June and in Aug/Sep. June Chrysaora were immature typically, and had less carbon content. Aug/Sept ones were mature and larger
# June: 28 g WW individual with 0.202% C ==> 57 mg C on average
# Aug/Sept: 700 g WW individual with 0.280% C ==> 1960 mg C on average
ncc$numeric_density <- NA
ncc[ncc$MONTH %in% 6:7,]$numeric_density <- ncc[ncc$MONTH %in% 6:7,]$Biomass / 57
ncc[ncc$MONTH %in% 8:10,]$numeric_density <- ncc[ncc$MONTH %in% 8:10,]$Biomass / 1960

names(ncc) <- tolower(names(ncc))
ncc <- rename(ncc, replace=c("biomass" = "Biomass.mgCm3"))
ncc$rank_phylum <- "Cnidaria"
ncc <- ncc[,which(names(ncc) != "seanett")]
#ncc <- cbind(ncc, data.frame(local_time=NA, taxon="jellies", rank_phylum="Cnidaria", rank_class="Scyphozoa", rank_order="Semaeostomeae", rank_family="Pelagiidae", rank_genus="Chrysaora", rank_species="fuscescens", GeoStd.Biomass=NA, Long.Prov="CCAL"))

ncc$lat = floor(ncc$latitude) + 0.5
ncc$lon = floor(ncc$longitude) + 0.5
# one row per phyla per lat/lon
ncc <- ddply(ncc, rank_phylum~lat+lon, function(x){
  latitude <- mean(x$latitude)
  longitude <- mean(x$longitude)
  Biomass.mgCm3 <- mean(x$Biomass.mgCm3)
  numeric_density <- mean(x$numeric_density)
  n <- sum(x$n)
  return(data.frame(latitude, longitude, Biomass.mgCm3, numeric_density, n))
})

# 4) ADD SALP DATA FROM BATS AND PALMER LTER
# add BATS salps data
# year 1994-2011
# salps only
bats <- data.frame(latitude=31.6667, longitude=-64.01667, lat=31.5, lon=-64.5, rank_phylum="Chordata", numeric_density=2.635, Biomass.mgCm3=0.0423, n=1)

# 5) ADD PALMER LTER SALP DATA
# download data from: https://oceaninformatics.ucsd.edu/datazoo/catalogs/pallter/datasets/212

season_mapping = data.frame(Month=1:12, Season=c(rep("DJF",2), rep("MAM",3), rep("JJA",3), rep("SON",3),"DJF"))

palmerLTER <- read.csv("raw_data/PalmerLTER_zooplankton_1993_2008.csv", stringsAsFactors=F) # "historical data"
names(palmerLTER) <- c("studyName", "CruiseTow", "CruiseName", "DateGMT", "TimeStartGMT", "TimeLocalCLST", "LatitudeStart", "LongitudeStart", "TimeEndGMT", 
                     "LatitudeEnd", "LongitudeEnd", "TowDurationmin", "TowType", "DepthTarget", "DepthMaximumm", "VolumeFilteredM3m", "SiphonNum", "SiphonVol", 
                     "CnidariaNum", "CnidariaVol", "CtenophNum","CtenophVol", "SalpEmbNum", "SalpEmbVol","SalpAggNum", "SalpAggVol")
palmerLTER <- palmerLTER[,c("DateGMT", "LatitudeStart", "LongitudeStart", "LatitudeEnd", "LongitudeEnd", "SalpEmbNum", "SalpEmbVol", "SalpAggNum", "SalpAggVol")]
palmerLTER <- palmerLTER[palmerLTER$DateGMT != "",]
palmerLTER$DateGMT <- as.Date(palmerLTER$DateGMT, format = "%m/%d/%y")
palmerLTER$Year <- format(palmerLTER$DateGMT, "%Y")
palmerLTER$Month <- as.numeric(format(palmerLTER$DateGMT, "%m"))
palmerLTER <- join(palmerLTER, season_mapping, by="Month")

palmerLTER[which(palmerLTER$SalpEmbNum == -1 & palmerLTER$SalpEmbVol == -1), c("SalpEmbNum","SalpEmbVol")] <- NA
palmerLTER[which(palmerLTER$SalpAggNum == -1 & palmerLTER$SalpAggVol == -1), c("SalpAggNum","SalpAggVol")] <- NA

palmerLTER <- palmerLTER[!apply(palmerLTER[,4:7],1, function(x){all(is.na(x))}),]


# fill in rows where Vol or Num is negative with mean value
non_negs = which(palmerLTER$SalpEmbNum > 0 & palmerLTER$SalpEmbVol > 0)
mean_Emb_vol = mean(palmerLTER$SalpEmbVol[non_negs] / palmerLTER$SalpEmbNum[non_negs])
palmerLTER[which(palmerLTER$SalpEmbNum == -1),"SalpEmbNum"] <- palmerLTER[which(palmerLTER$SalpEmbNum == -1),"SalpEmbVol"] / mean_Emb_vol
palmerLTER[which(palmerLTER$SalpEmbVol == -1),"SalpEmbVol"] <- palmerLTER[which(palmerLTER$SalpEmbVol == -1),"SalpEmbNum"] * mean_Emb_vol

non_negs = which(palmerLTER$SalpAggNum > 0 & palmerLTER$SalpAggVol > 0)
mean_Agg_vol = mean(palmerLTER$SalpAggVol[non_negs] / palmerLTER$SalpAggNum[non_negs])
palmerLTER[which(palmerLTER$SalpAggNum == -1),"SalpAggNum"] <- palmerLTER[which(palmerLTER$SalpAggNum == -1),"SalpAggVol"] / mean_Agg_vol
palmerLTER[which(palmerLTER$SalpAggVol == -1),"SalpAggVol"] <- palmerLTER[which(palmerLTER$SalpAggVol == -1),"SalpAggNum"] * mean_Agg_vol


palmerLTER$Total.Salp.Vol <- rowSums(palmerLTER[,c("SalpEmbVol", "SalpAggVol")], na.rm=T) # ml / 1000 m^3
palmerLTER$Total.Salp.Num <- rowSums(palmerLTER[,c("SalpEmbNum", "SalpAggNum")], na.rm=T) # number / 1000 m^3

# mean individual volume, excluding negatives (must be times in which the value wasn't recorded)
non_negs = which(palmerLTER$Total.Salp.Num >0 & palmerLTER$Total.Salp.Vol > 0)
mean_salp_vol = mean(palmerLTER$Total.Salp.Vol / palmerLTER$Total.Salp.Num, na.rm=TRUE) # units: ml / indiv.

# determine center lat/lon coordinate of each tow
palmerLTER$Latitude = rowMeans(palmerLTER[,c("LatitudeStart", "LatitudeEnd")])
palmerLTER$Longitude = rowMeans(palmerLTER[,c("LongitudeStart", "LongitudeEnd")])

palmerLTER$lat <- floor(-palmerLTER$Latitude) + 0.5
palmerLTER$lon <- floor(-palmerLTER$Longitude) + 0.5

# take the daily sums over the grid cell
palmerAgg <- ddply(palmerLTER, ~DateGMT+Year+Month+Season+lat+lon, function(x){
  Sum.Daily.Vol <- sum(x$Total.Salp.Vol)
  Sum.Daily.Num <- sum(x$Total.Salp.Num)
  n = nrow(x)
  x <- data.frame(Sum.Daily.Vol, Sum.Daily.Num, n)
  return(x)
})

# average over every month
palmerAgg <- ddply(palmerAgg, ~Year+Month+Season+lat+lon, function(x){
  Mean.Total.Vol <- mean(x$Sum.Daily.Vol)
  Mean.Total.Num <- mean(x$Sum.Daily.Num)
  x <- data.frame(Mean.Total.Vol, Mean.Total.Num, n=sum(x$n))
  return(x)
})

# average over the years
palmerAgg <- ddply(palmerAgg, ~Month+Season+lat+lon, function(x){
  Mean.Total.Vol <- mean(x$Mean.Total.Vol)
  Mean.Total.Num <- mean(x$Mean.Total.Num)
  x <- data.frame(Mean.Total.Vol, Mean.Total.Num, n=sum(x$n))
  return(x)
})

# average over seasons
palmerAgg <- ddply(palmerAgg, ~Season+lat+lon, function(x){
  Mean.Total.Vol <- mean(x$Mean.Total.Vol)
  Mean.Total.Num <- mean(x$Mean.Total.Num)
  x <- data.frame(Mean.Total.Vol, Mean.Total.Num, n=sum(x$n))
  return(x)
})

# average per lat/lon
palmerAgg <- ddply(palmerAgg, ~lat+lon, function(x){
  Mean.Total.Vol <- mean(x$Mean.Total.Vol)
  Mean.Total.Num <- mean(x$Mean.Total.Num)
  x <- data.frame(Mean.Total.Vol, Mean.Total.Num, n=sum(x$n))
  return(x)
})

palmerAgg$Biomass <- palmerAgg$Mean.Total.Vol*1.02686 # from ml/1000 m^3 to mg WW / m^3
palmerAgg$Biomass <- palmerAgg$Biomass * 0.04 # DW as a % of WW. used values from Ikeda and Bruce (1986), Huntley et al. (1989), Reinke (1987)
palmerAgg$Biomass <- palmerAgg$Biomass * 0.11 # C as a % of WW - mean of 10.86
palmerAgg$numeric_density <- palmerAgg$Mean.Total.Num / 1000 # original units in #/1000 m^3

palmer <- data.frame(lat=palmerAgg$lat, lon=palmerAgg$lon, latitude=palmerAgg$lat, longitude=palmerAgg$lon,  rank_phylum="Chordata", numeric_density=palmerAgg$numeric_density, Biomass.mgCm3=palmerAgg$Biomass, n=palmerAgg$n)

# Add in southern ocean krillbase data
# download data from: https://data.bas.ac.uk/full-record.php?id=GB/NERC/BAS/PDC/00915
krillbase <- read.csv("raw_data/krillbase.csv", as.is=TRUE)

# exclude data from PalmerLTER
krillbase <- krillbase[which (! str_detect(krillbase$Station,"pal")),]
krillbase <- krillbase[,! str_detect(names(krillbase),"krill")]

# convert from column integrated invidiauls to individuals / m^3
krillbase <- krillbase[which(! is.na(krillbase$No..ofsalpsunder.1m2)),]

krillbase$samplingdepth = krillbase$Bottomsamplingdepth..m. - krillbase$Topsamplingdepth..m.
krillbase[which(krillbase$samplingdepth == 0),"samplingdepth"] = 200 # mean
krillbase$salps.cnt.m3 = krillbase$No..ofsalpsunder.1m2 / krillbase$samplingdepth

krillbase$lat = floor(krillbase$Latitude) + 0.5
krillbase$lon = floor(krillbase$Longitude) + 0.5
krillbase$Date = as.Date(krillbase$Date, format = "%d-%b-%Y")
krillbase$Month = as.numeric(format(krillbase$Date, "%m"))
krillbase$Year = as.numeric(format(krillbase$Date, "%Y"))
krillbase = join(krillbase, season_mapping)

# sequentially aggregating the data
krillbaseAgg = ddply(krillbase, ~Year+Month+Season+lat+lon, function(x){
  res=data.frame(salps.cnt.m3=mean(x$salps.cnt.m3), n=nrow(x))
  return(res)
})

krillbaseAgg = ddply(krillbaseAgg, ~Month+Season+lat+lon, function(x){
  res=data.frame(salps.cnt.m3=mean(x$salps.cnt.m3), n=sum(x$n))
  return(res)
})

krillbaseAgg = ddply(krillbaseAgg, ~Season+lat+lon, function(x){
  res=data.frame(salps.cnt.m3=mean(x$salps.cnt.m3), n=sum(x$n))
  return(res)
})

# ggplot(krillbaseAgg) + geom_point(aes(x=lon, y=lat, color=Season, size=salps.cnt.m2), alpha=0.5) + theme_bw()

krillbaseAgg = ddply(krillbaseAgg, ~lat+lon, function(x){
  res=data.frame(salps.cnt.m3=mean(x$salps.cnt.m3), n=sum(x$n))
  return(res)
})

# convert counts to biomass, use data from palmerLTER to convert
palmer_mean_indiv_biomass = mean(palmerAgg$Biomass / palmerAgg$numeric_density, na.rm=TRUE)
krillbaseAgg$Biomass <- krillbaseAgg$salps.cnt.m3 * palmer_mean_indiv_biomass

southernocean <- data.frame(latitude=krillbaseAgg$lat, longitude=krillbaseAgg$lon, lat=krillbaseAgg$lat, lon=krillbaseAgg$lon, rank_phylum="Chordata", numeric_density=krillbaseAgg$salps.cnt.m3, Biomass.mgCm3=krillbaseAgg$Biomass, n=krillbaseAgg$n)
southernocean <- rbind(palmer, southernocean)

krillbaseDups = which(duplicated(southernocean[,c("latitude","longitude")]))
palmerDups = which(duplicated(southernocean[,c("latitude","longitude")], fromLast=TRUE))
# do we average the two sets of data?

# put together
rd <- rbind(rd, bats, southernocean, gom, ncc)

rd <- ddply(rd, rank_phylum~lat+lon, function(x){
  latitude <- mean(x$latitude)
  longitude <- mean(x$longitude)
  numeric_density <- mean(x$numeric_density, na.rm=TRUE)
  Biomass.mgCm3 <- mean(x$Biomass.mgCm3, na.rm=TRUE)
  n <- sum(x$n)
  return(data.frame(latitude, longitude, numeric_density, Biomass.mgCm3, n))
}, .progress="text")



## { Process raw jelly data --------------------------------------------
# deal with the cases where no numeric density is recorded
# first calculate an average biomass per individual per row
rd$stdBiomass <- rd$Biomass.mgCm3 / rd$numeric_density

# change the Inf values to NA
rd[which(!is.finite(rd$stdBiomass)), "stdBiomass"] <- NA

# calculate an average size for the three groups
avg <- ddply(rd, ~rank_phylum, function(x){
  x <- x[complete.cases(x),]
  return(mean(x$stdBiomass, na.rm=TRUE))
})

for (i in 1:3){ rd[which(is.na(rd$stdBiomass) & rd$rank_phylum==avg$rank_phylum[i]),"stdBiomass"] <- avg$V1[i] }

# calculate total Biomass values from average stdBiomass
temp <- which(is.na(rd$Biomass.mgCm3))
rd[temp,"Biomass.mgCm3"] <- rd[temp,"stdBiomass"] * rd[temp,"numeric_density"] 

# recalculate a numeric density for the rows where you don't have this data based on the average biomass size for each group (Chordates, ctenophores, cnidarians)
temp <- which(is.na(rd$numeric_density))
rd[temp,"numeric_density"] <- rd[temp,"Biomass.mgCm3"] /rd[temp,"stdBiomass"]

# }

write.csv(rd, "data/0-baseline/jelly_biomass_1_deg_grid.csv", row.names=F)

# } 
