#!/usr/bin/env Rscript
#
#             Process raw zooplankton files
#         - output from GFDL COBALT model (Stock et al. 2014)
#
#
#---------------------------------------------------------

library("plyr")
library("reshape2")
library("ncdf4")
library("stringr")

## ------------------ FUNCTIONS ------------------------
source(file = "functions.R")

## ------------------ LOADING MODEL DATA ---------------------------
# load raw netCDF data
nc <- nc_open("raw_data/cobalt/zooplankton_data_luo_ann.nc")
print(nc)
# units:
# biomass (mg C m-2 over the top 200m) 
# production (mg C m-2 day-1 over the top 200m)


# SMALLZOO_BIOMASS
small_biomass <- ncvar_get(nc, "SMALLZOO_BIOMASS")
# MEDZOO_BIOMASS
med_biomass <- ncvar_get(nc, "MEDZOO_BIOMASS")
# LARGEZOO_BIOMASS
large_biomass <- ncvar_get(nc, "LARGEZOO_BIOMASS")
# SMALLZOO_PROD
small_prod <- ncvar_get(nc, "SMALLZOO_PROD")
# MEDZOO_PROD
med_prod <- ncvar_get(nc, "MEDZOO_PROD")
# LARGEZOO_PROD
large_prod <- ncvar_get(nc, "LARGEZOO_PROD")

lon <- nc$dim$XT_OCEAN$vals
lat <- nc$dim$YT_OCEAN$vals

## ----------- PROCESSING MODEL DATA -------------------------
#### Biomass
small_biomass_df <- collapse.into.1deg(process.nc.biomass.data(small_biomass, times = 'ann'))
med_biomass_df <- collapse.into.1deg(process.nc.biomass.data(med_biomass, times = 'ann'))
large_biomass_df <- collapse.into.1deg(process.nc.biomass.data(large_biomass, times = 'ann'))

small_biomass_df$ann = small_biomass_df$ann / 200 # convert to mg C m-3
med_biomass_df$ann = med_biomass_df$ann / 200
large_biomass_df$ann = large_biomass_df$ann / 200

# Write to csv
write.csv(small_biomass_df, "data/gfdl_small_zoo_biomass_1deg.csv", row.names=F)
write.csv(med_biomass_df, "data/gfdl_med_zoo_biomass_1deg.csv", row.names=F)
write.csv(large_biomass_df, "data/gfdl_large_zoo_biomass_1deg.csv", row.names=F)
# biomass values are daily mean biomass integrated over the top 200 m, mg C m^-3

rzoomass <- cbind(small_biomass_df[,1:2], sm_zoomass=small_biomass_df$ann, md_zoomass=med_biomass_df$ann, lg_zoomass=large_biomass_df$ann)
save(rzoomass, file="data/gfdl_all_zoo_biomass_1deg.Rdata")

#### Production
small_prod_df <- collapse.into.1deg(process.nc.biomass.data(small_prod, times = 'ann'))
med_prod_df <- collapse.into.1deg(process.nc.biomass.data(med_prod, times = 'ann'))
large_prod_df <- collapse.into.1deg(process.nc.biomass.data(large_prod, times = 'ann'))
# values of production are in mg C m^-2 d^-1 integrated over the top 200 m

small_prod_df$ann = small_prod_df$ann / 200 # now convert to mg C m-3 d-1
med_prod_df$ann = med_prod_df$ann / 200
large_prod_df$ann = large_prod_df$ann / 200

mean(med_prod_df$ann + large_prod_df$ann, na.rm=TRUE) * 365 * 200 * 3.4E14 / 1000
# --> 4.063264e+15

# Write to csv
write.csv(small_prod_df, "data/gfdl_small_zoo_prod_1deg.csv", row.names=F)
write.csv(med_prod_df, "data/gfdl_med_zoo_prod_1deg.csv", row.names=F)
write.csv(large_prod_df, "data/gfdl_large_zoo_prod_1deg.csv", row.names=F)

rzooprod <- cbind(small_prod_df[,1:2], sm_zooprod=small_prod_df$ann, md_zooprod=med_prod_df$ann, lg_zooprod=large_prod_df$ann)
rzooprod$total_sum <- rzooprod$sm_zooprod + rzooprod$md_zooprod + rzooprod$lg_zooprod

save(rzooprod, file="data/gfdl_all_zoo_production_1deg.Rdata")

## process total production over the ocean
load("data/surfaceArea_volume_1-deg_filled_in.Rdata")
all_prod <- rzooprod
names(all_prod) <- c("lat", "lon", "sm_ann", "md_ann", "lg_ann", "all_ann")
all_prod <- all_prod[complete.cases(all_prod),]
all_prod <- join(all_prod, av, by=c("lat", "lon"), type="left")
all_prod$sm_ann_sum <- all_prod$sm_ann * all_prod$area.m2
all_prod$md_ann_sum <- all_prod$md_ann * all_prod$area.m2
all_prod$lg_ann_sum <- all_prod$lg_ann * all_prod$area.m2
sum(all_prod[, c("sm_ann_sum", "md_ann_sum", "lg_ann_sum")], na.rm=T) * 365
sum(all_prod[, c("sm_ann_sum", "md_ann_sum", "lg_ann_sum")], na.rm=T) * 365 / 1000 # in grams
# -->  3.477625e+15
mean(rowSums(all_prod[,c("md_ann", "lg_ann")])) * 365 / 1000 * 3.4E14
# --> 4.07402e+15


all_biomass <- cbind(small_biomass_df[,c("lat", "lon", "ann")], med_biomass_df$ann, large_biomass_df$ann)
names(all_biomass) <- c("lat", "lon", "sm_ann", "md_ann", "lg_ann")
all_biomass <- all_biomass[complete.cases(all_biomass),]
all_biomass <- join(all_biomass, av, by=c("lat", "lon"), type="left")
all_biomass$sm_ann_sum <- all_biomass$sm_ann * all_biomass$area.m2
all_biomass$md_ann_sum <- all_biomass$md_ann * all_biomass$area.m2
all_biomass$lg_ann_sum <- all_biomass$lg_ann * all_biomass$area.m2
all_biomass$sum <- rowSums(all_biomass[,c("sm_ann", "md_ann", "lg_ann")])
mean(all_biomass$sum) / 1000 * 3.4E14
# -->6 .038434e+14

# }
