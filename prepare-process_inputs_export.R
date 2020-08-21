#!/usr/bin/env Rscript
#
#      Prepare COBALT model export files
#
#-------------------------------------------------------

library("plyr")
library("reshape2")
library("ncdf4")
library("matlab")
library("oce")
library("ggplot2")
library("stringr")

## ------------------ FUNCTIONS ------------------------
source(file = "functions.R")

## -----------------------------------------------------

nc <- nc_open("raw_data/cobalt/export_data_luo_ann.nc")

print(nc)


# lat-lon grid
lon <- ncvar_get(nc, "XT_OCEAN")
lat <- ncvar_get(nc, "YT_OCEAN")

# grab variables; this represents POC flux at 100m, 1000m, and at the seafloor bottom
det_100m <- ncvar_get(nc, "FCDET_100")
det_1000m <- ncvar_get(nc, "FCDET_1000")
det_btm <- ncvar_get(nc, "FCDET_BTM")

# units?? should be mg C m^-2 d^-1

# process data and summarize into seasonal and yearly flux values
det_100m_df <- process.nc.biomass.data(det_100m, times = 'ann')
det_1000m_df <- process.nc.biomass.data(det_1000m, times = 'ann')
det_btm_df <- process.nc.biomass.data(det_btm, times = 'ann')

# collapse into 1 degree
det_100m_df <- collapse.into.1deg(det_100m_df)
det_1000m_df <- collapse.into.1deg(det_1000m_df)
det_btm_df <- collapse.into.1deg(det_btm_df)

# visualizations
plotting("ann", det_100m_df, size=1)

plotting("ann", det_1000m_df, size=1)

plotting("ann", det_btm_df, size=1)

# join all together and save
cobalt_det <- data.frame(lon=det_100m_df$lon, lat=det_100m_df$lat, det_100=det_100m_df$ann, det_1000=det_1000m_df$ann, det_btm=det_btm_df$ann)

save(cobalt_det, file="data/gfdl_POCflux_100m_1000m_seafloor.Rdata")
# units of mg C m^-2 d^-1
