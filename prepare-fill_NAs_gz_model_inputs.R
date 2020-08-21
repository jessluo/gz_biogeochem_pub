# 
#          Standardize grids and fill in NAs 
#            for upper ocean model input files
#
#-----------------------------------------------------

# reading packages
library(plyr)
library(reshape2)
library(stringr)
library(oce)


rd <- read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv", stringsAsFactors=FALSE)

LATLON <- c("lat", "lon")


##{ Read in data -------------------------------------------------------
load("data/sst.ltm.1-deg.annual.Rdata")

# COBALT Primary Productivity, units of mg C m^-3 d^-1, average over top 200 m
load("data/gfdl_primary_prod_for_salps_1deg.Rdata")
pp_salps <- rename(pp_salps, replace=c("ann"="pp_ch"))

load("data/gfdl_primary_prod_all_1deg.Rdata")
pp <- rename(pp, replace=c("ann"="pp"))

# read in phytoplankton biomass, units are mg C m-3, average over top 200 m
load("data/gfdl_phytoplankton_for_salps_1deg.Rdata")
salp_phyto_df <- rename(salp_phyto_df, replace=c("ann"="phytoC_ch"))

load("data/gfdl_all_phytoplankton_1deg.Rdata")
allphyto_df <- rename(allphyto_df, replace=c("ann"="phytoC"))

# read in COBALT zooplankton production, mg C m^-3 d^-1, average over top 200 m
load("data/gfdl_all_zoo_production_1deg.Rdata")

# read in COBALT zooplankton biomass, mg C m-3, average over top 200 m
load("data/gfdl_all_zoo_biomass_1deg.Rdata")

# # read in the surface area & volume of the 1 deg grid
# av <- read.csv("data/surfaceArea_volume_1-deg.csv")

# 
rzooprod <- rzooprod[complete.cases(rzooprod),]

rsst <- join(rzooprod[,LATLON], rsst, by=LATLON, type="left")
rzoomass <- join(rzooprod[,LATLON], rzoomass, by=LATLON, type="left")
pp_salps <- join(rzooprod[,LATLON], pp_salps, by=LATLON, type="left")
phyto_salps <- join(rzooprod[,LATLON], salp_phyto_df, by=LATLON, type="left")
pp <- join(rzooprod[,LATLON], pp, by=LATLON, type="left")
phyto <- join(rzooprod[,LATLON], allphyto_df, by=LATLON, type="left")

# }

## { Process raw jelly data --------------------------------------------
# deal with the cases where no numeric density is recorded
# first calculate an average biomass per individual per row
rd$stdSize <- rd$Biomass.mgCm3 / rd$numeric_density

# calculate an average size for the three groups
avg <- ddply(rd, ~rank_phylum, function(x){
  return(mean(x$stdSize)) #   return(gm_mean(x$stdSize))
})

# }

##{ Standardize and fill in data for other forcings ----------------------------------------------------------------------------- 

## FILL IN MISSING VALUES ##
fill_in_missing_grids <- function(data, jellyrd=rd, var="ann"){
  data=data[,c("lat","lon",var)]
  tmpdata = join(x=jellyrd[,c("lat","lon")], y=data, by=c("lat","lon"), type="left")
  
  nodata = tmpdata[which(!complete.cases(tmpdata)),]
  nodata = nodata[!duplicated(nodata),]
  
  if(nrow(nodata) == 0){
    return(data)
  }
  
  varlength = length(var)
  if(varlength==1){
    data$sum=data[,var]
  } else if(varlength > 1){
    data$sum=rowSums(data[,var])
  }
  
  tmp <- data[complete.cases(data),]
  
  # use loop to fill in missing values
  # if there are values within a nearly 3-degree radius, then fill in the average of the values within a 3-degree grid cell
  # otherwise look for the closest grid cell with model output and fill it in
  
  for (i in 1:nrow(nodata)){
    loopdata <- data[which(data$lat <= (nodata$lat[i] + 1.5) & data$lat >= (nodata$lat[i] - 1.5) & 
                             data$lon <= (nodata$lon[i] + 1.5) & data$lon >= (nodata$lon[i] - 1.5)),]
    
    if(all(is.na(loopdata$sum))){
      
      print(paste0('searching nearby, row:', i))
      
      loopdata <- tmp[which(geodDist(tmp$lon, tmp$lat, nodata$lon[i], nodata$lat[i]) == 
                              min(geodDist(tmp$lon, tmp$lat, nodata$lon[i], nodata$lat[i]))),]
      
    }
    
    if(varlength==1){
      nodata[i,var] <- mean(loopdata[,var], na.rm=TRUE)
    } else if (varlength > 1){
      nodata[i,var] <- colMeans(loopdata[,var], na.rm=TRUE)
    }
    # assign back to data
    idx <- which(tmpdata$lat == nodata$lat[i] & tmpdata$lon == nodata$lon[i])
    tmpdata[idx,] <- nodata[i,]  
  }
  
  return(tmpdata)
}

rsst <- fill_in_missing_grids(rsst, rd, var="ann")
rzoomass <- fill_in_missing_grids(rzoomass, rd, var=c("sm_zoomass","md_zoomass","lg_zoomass"))
rzooprod <- fill_in_missing_grids(rzooprod, rd, var=c("sm_zooprod","md_zooprod","lg_zooprod"))
phyto_salps <- fill_in_missing_grids(phyto_salps, rd, var=c("phytoC_ch"))
pp_salps <- fill_in_missing_grids(pp_salps, rd, var=c("pp_ch"))
phyto <- fill_in_missing_grids(phyto, rd, var=c("phytoC"))
pp <- fill_in_missing_grids(pp, rd, var=c("pp"))

# save files
dir.create("data/gz_model_inputs/", showWarnings = FALSE)

save(rsst, file = "data/gz_model_inputs/rsst.Rdata")
save(rzoomass, file = "data/gz_model_inputs/rzoomass.Rdata")
save(rzooprod, file = "data/gz_model_inputs/rzooprod.Rdata")
save(phyto_salps, file = "data/gz_model_inputs/phyto_salps.Rdata")
save(pp_salps, file = "data/gz_model_inputs/pp_salps.Rdata")
save(phyto, file = "data/gz_model_inputs/phyto.Rdata")
save(pp, file = "data/gz_model_inputs/pp.Rdata")

# }