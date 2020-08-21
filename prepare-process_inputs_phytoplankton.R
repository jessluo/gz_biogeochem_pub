#!/usr/bin/env Rscript
#
#           Global input data processing - for phytoplankton
#
#------------------------------------------------------------------------------------------

# load libraries
library("plyr")
library("reshape2")
library("ncdf4")

# load functions
source("functions.R")


##{ Process COBALT NPP data --------------------------------------------
# load raw netCDF data
nc <- nc_open("raw_data/cobalt/primary_production_luo_ann.nc")

# units:
# mg C m-2 day-1

print(nc)

lon <- ncvar_get(nc, "XT_OCEAN")
lat <- ncvar_get(nc, "YT_OCEAN") # note uneven lats

# get data
prim_prod <- ncvar_get (nc, "PRIMARY_PROD") # num [1:360, 1:200]

# process data
pp <- process.nc.biomass.data(prim_prod, times = "ann")

# collapse into 1 degree bins
pp <- collapse.into.1deg(pp)
pp <- pp[,c("lon", "lat", "ann")]

# convert from mg C m-2 day-1 to mg C m-3 day-1 (assume 200 m depth for all; consistent with Charlie's conversion originally)
pp$ann = pp$ann / 200

save(pp, file="data/gfdl_primary_prod_all_1deg.Rdata")
# }

### {Process separated NPP data -----------------------------------------------
# load raw netCDF data
nc <- nc_open("raw_data/cobalt/primary_production_separated_luo_ann.nc")

# units:
# mg C m-2 day-1
print(nc)

lon <- ncvar_get(nc, "xt_ocean")
lat <- ncvar_get(nc, "yt_ocean") # note uneven lats

# get data
smp_prod <- ncvar_get (nc, "sm_phyto_prod") # num [1:360, 1:200]
lgp_prod <- ncvar_get (nc, "lg_phyto_prod") # num [1:360, 1:200]
diaz_prod <- ncvar_get (nc, "diaz_prod") # num [1:360, 1:200]

# pull together
# smp + 0.5 lgp + 0.5 diaz
prim_prod_salps = smp_prod + 0.5 * lgp_prod + 0.5 * diaz_prod

# process data
pp_salps <- process.nc.biomass.data(prim_prod_salps, times = "ann")
pp_salps <- collapse.into.1deg(pp_salps)

# convert from mg C m-2 day-1 to mg C m-3 day-1 (assume 200 m depth for all; consistent with Charlie's conversion originally)
pp_salps$ann = pp_salps$ann / 200

# save
save(pp_salps, file="data/gfdl_primary_prod_for_salps_1deg.Rdata")

# }


##{ Process COBALT Phytoplankton Biomass data --------------------------------------------
# load raw netCDF data
nc <- nc_open("raw_data/cobalt/phytoplankton_luo_ann.nc")
print(nc)

# file has data on three phytoplankton types: small phytoplankton, large phytoplankton, and diazotrophs
# averaged 20 depth dimensions, over time

# units:
# mg C m-2

# all phyto instead
allphyto <- ncvar_get(nc, "all_phyto")

lon <- ncvar_get(nc, "xt_ocean")
lat <- ncvar_get(nc, "yt_ocean") # note uneven lats

# process data, & collapse into 1 degree bins
allphyto_df <- collapse.into.1deg(process.nc.biomass.data(allphyto, times = 'ann'))

# convert from mg C m-2 to mg C m-3 (assume 200 m depth for all; consistent with Charlie's conversion)
allphyto_df$ann = allphyto_df$ann / 200

save(allphyto_df, file="data/gfdl_all_phytoplankton_1deg.Rdata")

### do this again but calculate a salp specific available phytoplankton biomass
# salp_phyto = small phyto + (1/2) large phyto + (1/2) diazotrophs
smp <- ncvar_get(nc, "sm_phyto")
lgp <- ncvar_get(nc, "lg_phyto")
diaz <- ncvar_get(nc, "diaz")

salp_phyto = smp + 0.5 * lgp + 0.5 * diaz


# process data, & collapse into 1 degree bins
salp_phyto_df = process.nc.biomass.data(salp_phyto, times = 'ann')
salp_phyto_df = collapse.into.1deg(salp_phyto_df)

# convert from mg C m-2 to mg C m-3 (assume 200 m depth for all; consistent with Charlie's conversion)
salp_phyto_df$ann = salp_phyto_df$ann / 200

save(salp_phyto_df, file="data/gfdl_phytoplankton_for_salps_1deg.Rdata")


# }




