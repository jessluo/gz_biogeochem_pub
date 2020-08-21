#!/usr/bin/env Rscript
#
#
#   Extract single point comparisons
#
#
# -------------------------------------------

# load libraries ----------------------------
library("plyr")
library("reshape2")
library("stringr")

# load data
args <- commandArgs(trailingOnly=TRUE)
C_ID <- as.numeric(args[1])
# C_ID = 1
CASES = c("0-baseline","1-high_biomass","2-low_biomass")
CASE = CASES[C_ID]

load("data/surfaceArea_volume_1-deg_filled_in.Rdata")
av$depth <- av$vol / av$area.m2

# extract lat's and lon's for each taxa
load(str_c("data/",CASE,"/mc_res/MC_ExportError_exprt.Rdata"))
cn_lat = exprt_all$Cnidarians$lat
cn_lon = exprt_all$Cnidarians$lon
ct_lat = exprt_all$Ctenophores$lat
ct_lon = exprt_all$Ctenophores$lon
ch_lat = exprt_all$Chordates$lat
ch_lon = exprt_all$Chordates$lon

# load biomes
load('data/biomes_chl_mixedlayer.Rdata')

biomes_cn = join(data.frame(lat=cn_lat, lon=cn_lon), biomes_df)
biomes_ct = join(data.frame(lat=ct_lat, lon=ct_lon), biomes_df)
biomes_ch = join(data.frame(lat=ch_lat, lon=ch_lon), biomes_df)

# load ensemble data
load(file=str_c("data/",CASE,"/mc_res/cn_cf_seafl_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/cn_eg_seafl_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/cn_cf_exprt_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/cn_eg_exprt_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/cn_cf_seq_d_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/cn_eg_seq_d_all.Rdata"))

load(file=str_c("data/",CASE,"/mc_res/ct_cf_seafl_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ct_eg_seafl_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ct_cf_exprt_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ct_eg_exprt_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ct_cf_seq_d_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ct_eg_seq_d_all.Rdata"))

load(file=str_c("data/",CASE,"/mc_res/ch_cf_seafl_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ch_eg_seafl_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ch_cf_exprt_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ch_eg_exprt_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ch_cf_seq_d_all.Rdata"))
load(file=str_c("data/",CASE,"/mc_res/ch_eg_seq_d_all.Rdata"))


## ------------- Extract Single Point -------------- ##

Mode <- function(df){
  return(names(which.max(table(df))))
}

singlept = as.data.frame(array(NA, dim=c(500,1)))

# cnidarian carcasses @ seafloor
# chesapeake bay
idx = which(cn_lat >= 36 & cn_lat < 40 & cn_lon >= -78 & cn_lon < -75)
chesapeake_bay = c(apply(cn_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 36 & av$lat < 40 & av$lon >= -78 & av$lon < -75),"depth"])
chesapeake_bay = chesapeake_bay * mean_depth
singlept$chesapeake_bay = NA; singlept$chesapeake_bay[1:length(chesapeake_bay)]=chesapeake_bay
singleptinfo <- data.frame(site="chesapeake_bay", taxa="cnid", type="carcass", biome=Mode(biomes_cn[idx,'biome']))

# norway fjord
idx = which(cn_lat >= 58 & cn_lat < 70 & cn_lon >= 5 & cn_lon < 20) # exclude baltic sea - no points north of Tromso; target western Norway fjords
norway_fjord = c(apply(cn_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE)) # take the mean over the grid points, retaining all the values from the diff. sinking speeds and ensembles
mean_depth = 200 # mean depth from dataset is 90 m, fjords are very very deep
norway_fjord = norway_fjord * mean_depth
singlept$norway_fjord = NA; singlept$norway_fjord[1:length(norway_fjord)]=norway_fjord
singleptinfo <- rbind(singleptinfo, data.frame(site="norway_fjord", taxa="cnid", type="carcass", biome=Mode(biomes_cn[idx,'biome'])))

# sea of japan
idx = which(cn_lat >= 33 & cn_lat < 50 & cn_lon >= 128 & cn_lon < 142)
sea_of_japan = c(apply(cn_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 33 & av$lat < 50 & av$lon >= 128 & av$lon < 142),"depth"])
sea_of_japan = sea_of_japan * mean_depth
singlept$sea_of_japan = NA; singlept$sea_of_japan[1:length(sea_of_japan)]=sea_of_japan
singleptinfo <- rbind(singleptinfo, data.frame(site="sea_of_japan", taxa="cnid", type="carcass", biome=Mode(biomes_cn[idx,'biome'])))

# gulf of oman
idx = which(cn_lat >= 22 & cn_lat < 26 & cn_lon >= 56 & cn_lon < 68)
gulf_of_oman = c(apply(cn_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 20 & av$lat < 26 & av$lon >= 56 & av$lon < 68),"depth"])
gulf_of_oman = gulf_of_oman * mean_depth
singlept$gulf_of_oman = NA; singlept$gulf_of_oman[1:length(gulf_of_oman)]=gulf_of_oman
singleptinfo <- rbind(singleptinfo, data.frame(site="gulf_of_oman", taxa="cnid", type="carcass", biome=Mode(biomes_cn[idx,'biome'])))

# # hudson bay
# idx = which(cn_lat >= 63 & cn_lat < 64 & cn_lon >= -81 & cn_lon < -79) # only one point with cnidarian obs
# hudson_bay = c(apply(cn_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
# mean_depth = mean(av[which(av$lat >= 63 & av$lat < 64 & av$lon >= -81 & av$lon < -79),"depth"])
# hudson_bay = hudson_bay * mean_depth
# singlept$hudson_bay = NA; singlept$hudson_bay[1:length(hudson_bay)]=hudson_bay
# singleptinfo <- rbind(singleptinfo, data.frame(site="hudson_bay", taxa="cnid", type="carcass"))

# tunicate carcasses & fecal pellets
# west mediterranean
idx = which(ch_lat >= 34 & ch_lat < 42 & ch_lon >= -5 & ch_lon < 15)
nw_med = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 34 & av$lat < 42 & av$lon >= -5 & av$lon < 15),"depth"])
nw_med = nw_med * mean_depth
singlept$nw_med = NA; singlept$nw_med[1:length(nw_med)]=nw_med
singleptinfo <- rbind(singleptinfo, data.frame(site="nw_med", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

# ivory coast
idx = which(ch_lat >= 0 & ch_lat < 6 & ch_lon >= -7.5 & ch_lon < 3)
ivory_coast = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 0 & av$lat < 6 & av$lon >= -7.5 & av$lon < 3),"depth"])
ivory_coast = ivory_coast * mean_depth
singlept$ivory_coast = NA; singlept$ivory_coast[1:length(ivory_coast)]=ivory_coast
singleptinfo <- rbind(singleptinfo, data.frame(site="ivory_coast", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

# tasman sea
idx = which(ch_lat >= -50 & ch_lat < -35 & ch_lon >= 148 & ch_lon < 180)
tasman_sea = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= -40 & av$lat < -37 & av$lon >= 148 & av$lon < 151),"depth"])
tasman_sea = tasman_sea * mean_depth
singlept$tasman_sea = NA; singlept$tasman_sea[1:length(tasman_sea)]=tasman_sea
singleptinfo <- rbind(singleptinfo, data.frame(site="tasman_sea", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

# North Pacific - Smith et al.
idx = which(ch_lat >= 30 & ch_lat < 40 & ch_lon >= -130 & ch_lon < -120)
stn_M_carcass = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 30 & av$lat < 40 & av$lon >= -130 & av$lon < -120),"depth"])
stn_M_carcass = stn_M_carcass * mean_depth
singlept$stn_M_carcass = NA; singlept$stn_M_carcass[1:length(stn_M_carcass)]=stn_M_carcass
singleptinfo <- rbind(singleptinfo, data.frame(site="stn_M_carcass", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 30 & ch_lat < 40 & ch_lon >= -130 & ch_lon < -120)
stn_M_fecalpellets = c(apply(ch_eg_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 30 & av$lat < 40 & av$lon >= -130 & av$lon < -120),"depth"])
stn_M_fecalpellets = stn_M_fecalpellets * mean_depth
singlept$stn_M_fecalpellets = NA; singlept$stn_M_fecalpellets[1:length(stn_M_fecalpellets)]=stn_M_fecalpellets
singleptinfo <- rbind(singleptinfo, data.frame(site="stn_M_fecalpellets", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

# west antarctic peninsula
idx = which(ch_lat >= -68 & ch_lat < -62 & ch_lon >= -70 & ch_lon < -56)
west_antarctica = c(apply(ch_eg_exprt_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= -68 & av$lat < 60 & av$lon >= -70 & av$lon < -50),"depth"])
west_antarctica = west_antarctica * mean_depth
singlept$west_antarctica = NA; singlept$west_antarctica[1:length(west_antarctica)]=west_antarctica
singleptinfo <- rbind(singleptinfo, data.frame(site="west_antarctica", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

# southern ocean - lazarev sea
idx = which(ch_lat >= -71 & ch_lat < -65 & ch_lon >= -5 & ch_lon < 5)
lazarev_sea = c(apply(ch_eg_exprt_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= -71 & av$lat < -65 & av$lon >= -5 & av$lon < 5),"depth"])
lazarev_sea = lazarev_sea * mean_depth
singlept$lazarev_sea = NA; singlept$lazarev_sea[1:length(lazarev_sea)]=lazarev_sea
singleptinfo <- rbind(singleptinfo, data.frame(site="lazarev_sea", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

# southern ocean - south atlantic
idx = which(ch_lat >= -55 & ch_lat < -50 & ch_lon >= -15 & ch_lon < -10)
southern_ocean_atl = c(apply(ch_eg_exprt_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= -55 & av$lat < -50 & av$lon >= -15 & av$lon < -10),"depth"])
southern_ocean_atl = southern_ocean_atl * mean_depth
singlept$southern_ocean_atl = NA; singlept$southern_ocean_atl[1:length(southern_ocean_atl)]=southern_ocean_atl
singleptinfo <- rbind(singleptinfo, data.frame(site="southern_ocean_atl", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

# NE pacific - Matsueda
idx = which(ch_lat >= 35 & ch_lat < 40 & ch_lon >= -130 & ch_lon < -125) # single point
ne_pacific = c(ch_eg_seafl_all[idx,,])
mean_depth = mean(av[which(av$lat >= 35 & av$lat < 40 & av$lon >= -130 & av$lon < -125),"depth"])
ne_pacific = ne_pacific * mean_depth
singlept$ne_pacific = NA; singlept$ne_pacific[1:length(ne_pacific)]=ne_pacific
singleptinfo <- rbind(singleptinfo, data.frame(site="ne_pacific", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

# North Atlantic
idx = which(ch_lat >= 37 & ch_lat < 40 & ch_lon >= -74 & ch_lon < -69)
atl_wiebe_cf = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 37 & av$lat < 40 & av$lon >= -74 & av$lon < -69),"depth"])
atl_wiebe_cf = atl_wiebe_cf * mean_depth
singlept$atl_wiebe_cf = NA; singlept$atl_wiebe_cf[1:length(atl_wiebe_cf)]=atl_wiebe_cf
singleptinfo <- rbind(singleptinfo, data.frame(site="atl_wiebe_cf", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 37 & ch_lat < 40 & ch_lon >= -74 & ch_lon < -69)
atl_wiebe_eg = c(apply(ch_eg_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 37 & av$lat < 40 & av$lon >= -74 & av$lon < -69),"depth"])
atl_wiebe_eg = atl_wiebe_eg * mean_depth
singlept$atl_wiebe_eg = NA; singlept$atl_wiebe_eg[1:length(atl_wiebe_eg)]=atl_wiebe_eg
singleptinfo <- rbind(singleptinfo, data.frame(site="atl_wiebe_eg", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 37 & ch_lat < 42 & ch_lon >= -75 & ch_lon < -68)
atl_madin_eg = c(apply(ch_eg_exprt_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 37 & av$lat < 42 & av$lon >= -75 & av$lon < -68),"depth"])
atl_madin_eg = atl_madin_eg * mean_depth
singlept$atl_madin_eg = NA; singlept$atl_madin_eg[1:length(atl_madin_eg)]=atl_madin_eg
singleptinfo <- rbind(singleptinfo, data.frame(site="atl_madin_eg", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 26 & ch_lat < 37 & ch_lon >= -78 & ch_lon < -73)
atl_caron_eg = c(apply(ch_eg_exprt_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 26 & av$lat < 37 & av$lon >= -78 & av$lon < -73),"depth"])
atl_caron_eg = atl_caron_eg * mean_depth
singlept$atl_caron_eg = NA; singlept$atl_caron_eg[1:length(atl_caron_eg)]=atl_caron_eg
singleptinfo <- rbind(singleptinfo, data.frame(site="atl_caron_eg", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

# BATS & sargasso sea - stone and steinberg 2016
idx = which(ch_lat >= 30 & ch_lat < 34 & ch_lon >= -66 & ch_lon < -62)
bats_cf = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 30 & av$lat < 34 & av$lon >= -65 & av$lon < -62),"depth"])
bats_cf = bats_cf * mean_depth
singlept$bats_cf = NA; singlept$bats_cf[1:length(bats_cf)]=bats_cf
singleptinfo <- rbind(singleptinfo, data.frame(site="bats_cf", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 30 & ch_lat < 34 & ch_lon >= -66 & ch_lon < -62)
bats_eg = c(apply(ch_eg_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 30 & av$lat < 34 & av$lon >= -65 & av$lon < -62),"depth"])
bats_eg = bats_eg * mean_depth
singlept$bats_eg = NA; singlept$bats_eg[1:length(bats_eg)]=bats_eg
singleptinfo <- rbind(singleptinfo, data.frame(site="bats_eg", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 22 & ch_lat < 40 & ch_lon >= -77 & ch_lon < -45)
sargasso_sea_cf = c(apply(ch_cf_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 22 & av$lat < 40 & av$lon >= -77 & av$lon < -45),"depth"])
sargasso_sea_cf = sargasso_sea_cf * mean_depth
singlept$sargasso_sea_cf = NA; singlept$sargasso_sea_cf[1:length(sargasso_sea_cf)]=sargasso_sea_cf
singleptinfo <- rbind(singleptinfo, data.frame(site="sargasso_sea_cf", taxa="chor", type="carcass", biome=Mode(biomes_ch[idx,'biome'])))

idx = which(ch_lat >= 22 & ch_lat < 40 & ch_lon >= -77 & ch_lon < -45)
sargasso_sea_eg = c(apply(ch_eg_seafl_all[idx,,], 2:3, mean, na.rm=TRUE))
mean_depth = mean(av[which(av$lat >= 22 & av$lat < 40 & av$lon >= -77 & av$lon < -45),"depth"])
sargasso_sea_eg = sargasso_sea_eg * mean_depth
singlept$sargasso_sea_eg = NA; singlept$sargasso_sea_eg[1:length(sargasso_sea_eg)]=sargasso_sea_eg
singleptinfo <- rbind(singleptinfo, data.frame(site="sargasso_sea_eg", taxa="chor", type="fecalpellets", biome=Mode(biomes_ch[idx,'biome'])))

singlept <- singlept[,which(names(singlept) != "V1")]

write.csv(singlept, str_c("data/",CASE,"/site_comparisons_allpoints.csv"), row.names=FALSE)
write.csv(singleptinfo, "data/site_comparisons_info.csv", row.names=FALSE)
