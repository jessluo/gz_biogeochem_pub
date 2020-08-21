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

PLOTTING <- FALSE

source("functions.R")


## { Load files ---------------------------------------------------------
# import seafloor depth
seafloor <- read.csv("data/woa_seafloor_depth.csv")
seafloor$latlon <- paste0(seafloor$lat, ", ", seafloor$lon)

# import export depth
export.depth <- read.csv("data/woa_export_depth.csv")

# import sequestration depth
seq.depth <- read.csv("data/woa_sequestration_depth.csv")

if(any(cf[,c("lat", "lon")]!=eg[,c("lat", "lon")])){print("WARNING: carcasses and fecal pellet data do not match")}

cf <- cbind(cf, Cnidaria_eg=eg$Cnidaria, Ctenophora_eg=eg$Ctenophora, Chordata_eg=eg$Chordata)
cf_matrix <- as.matrix(cf[,3:8])


# }

## { Import particle sinking df and match it up with raw data ------------------------------

load("data/dm_particle_sinking_matched_template.Rdata")

# check to make sure that all the same
if (! (all(which(dm$Cnidaria==9999)==which(dm$Cnidaria_eg==9999)) &
    all(which(dm$Ctenophora==9999)==which(dm$Ctenophora_eg==9999)) &
    all(which(dm$Chordata==9999)==which(dm$Chordata_eg==9999))))
  {print ("WARNING: sinking particles do not match")}

idx <- which(dm$Cnidaria == 9999)
dm[idx,c("lat", "lon", "Cnidaria", "Cnidaria_eg")] <- temp <- join(dm[idx,c("lat", "lon")], cf[,c("lat", "lon", "Cnidaria", "Cnidaria_eg")], by=c("lat", "lon"))

idx <- which(dm$Ctenophora == 9999)
dm[idx,c("lat", "lon", "Ctenophora", "Ctenophora_eg")] <- temp <- join(dm[idx,c("lat", "lon")], cf[,c("lat", "lon", "Ctenophora", "Ctenophora_eg")], by=c("lat", "lon"))

idx <- which(dm$Chordata == 9999)
dm[idx,c("lat", "lon", "Chordata", "Chordata_eg")] <- temp <- join(dm[idx,c("lat", "lon")], cf[,c("lat", "lon", "Chordata", "Chordata_eg")], by=c("lat", "lon"))

# }

## { Calculate jelly decomposition rates -------------------------------

# calculate the amount of jelly biomass that arrives at each depth 
dm <- ddply(dm, ~lat+lon+C, function(x){
  
  cnTMP <- which(!is.na(x$Cnidaria))

  if(length(cnTMP) != 0){
    idx <- which(x$depth >= x$depth[cnTMP])
    temp <- x[idx,]
    temp$m_ratio_mod <- temp$m_ratio/temp$m_ratio[1]
    biomass <- temp$Cnidaria[1]
    egestion <- temp$Cnidaria_eg[1]
    temp$Cnidaria <- temp$m_ratio_mod * biomass
    temp$Cnidaria_eg <- temp$m_ratio_mod * egestion
    x[idx,"Cnidaria"] <- temp$Cnidaria
    x[idx,"Cnidaria_eg"] <- temp$Cnidaria_eg
  }
  
  ctTMP <- which(!is.na(x$Ctenophora))
  
  if(length(ctTMP) != 0){
    idx <- which(x$depth >= x$depth[ctTMP])
    temp <- x[idx,]
    temp$m_ratio_mod <- temp$m_ratio/temp$m_ratio[1]
    biomass <- temp$Ctenophora[1]
    egestion <- temp$Ctenophora_eg[1]
    temp$Ctenophora <- temp$m_ratio_mod * biomass
    temp$Ctenophora_eg <- temp$m_ratio_mod * egestion
    x[idx,"Ctenophora"] <- temp$Ctenophora
    x[idx,"Ctenophora_eg"] <- temp$Ctenophora_eg
  }
  
  chTMP <- which(!is.na(x$Chordata))
  
  if(length(chTMP)){
    idx <- which(x$depth >= x$depth[chTMP])
    temp <- x[idx,]
    temp$m_ratio_mod <- temp$m_ratio/temp$m_ratio[1]
    biomass <- temp$Chordata[1]
    egestion <- temp$Chordata_eg[1]
    temp$Chordata <- temp$m_ratio_mod * biomass
    temp$Chordata_eg <- temp$m_ratio_mod * egestion
    x[idx,"Chordata"] <- temp$Chordata
    x[idx,"Chordata_eg"] <- temp$Chordata_eg
  }
  return(x)
})#, .progress="text")

# calculate the amount of biomass ending up past the export depth
d.export <- join(dm, export.depth, by = c("lat", "lon", "depth"), type = "right", match="all")


# TODO: MAKE THIS MORE EFFICIENT
for (i in 1:nrow(cf)){
  lat <- cf$lat[i]
  lon <- cf$lon[i]
  # if there is a value in cf for a taxa but NOT a value in d.export for a taxa, then import the first value below 100m
  for (j in c("Cnidaria", "Ctenophora", "Chordata")){
    if(!is.na(cf[i,j]) & is.na(d.export[d.export$lat==lat & d.export$lon==lon,j][1])){
      j_eg <- paste0(j, "_eg")
      # subset
      temp <- dm[dm$lat==lat & dm$lon==lon, c("lat", "lon", "depth", "C", j, j_eg)]
      temp <- temp[complete.cases(temp),]
      
      # subset d.export and join
      temp2 <- d.export[d.export$lat==lat & d.export$lon==lon,c("lat", "lon", "C", j, j_eg)]
      temp2 <- join(temp2[c("lat", "lon", "C")], temp, by=c("lat", "lon", "C"), match="first")
      temp2 <- temp2[,which(names(temp2) != "depth")]
      
      # reassign back into d.export
      d.export[d.export$lat==lat & d.export$lon==lon,c("lat", "lon", "C", j, j_eg)] <- temp2
    }    
  }
}


# calculate the amount of biomass ending up past the sequestration depth
seq.depth <- join(dm, seq.depth, by = c("lat", "lon", "depth"), type = "right", match="all")

# calculate the amount of biomass ending up in the seafloor
seafloor <- ddply(dm, ~lat+lon+C, function(x){
  final.row <- x[nrow(x),]
  return(final.row)
}, .progress="text")


# write to disk
# save(dm, file=str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_sinking.Rdata")) # take off comments later
rm("dm")
save(seafloor, file=str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_reaching_seafloor.Rdata"))
save(d.export, file=str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_100m.Rdata"))
save(seq.depth, file=str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_seq_depth.Rdata"))


# don't run this for the MC error propagation
if(!ENSEMBLE_RUN){
  
  # calculate a biomass mean and standard deviation
  exportflux_mean_sd <- function(x){
    require(plyr)
    mean <- ddply(x, ~lat+lon, function(y){
      Cnidaria_cf <- mean(y[y$C %in% c(1000, 1100, 1200),"Cnidaria"])
      Cnidaria_eg <- mean(y[y$C == 100,"Cnidaria_eg"])
      Ctenophora_cf <- mean(y[y$C %in% c(800, 900, 1000),"Ctenophora"])
      Ctenophora_eg <- mean(y[y$C == 100,"Ctenophora_eg"])
      Chordata_cf <- mean(y[y$C %in% c(800, 900, 1000, 1100, 1200),"Chordata"])
      Chordata_eg <- mean(y[y$C %in% c(100, 800, 1200),"Chordata_eg"])
      Cnidaria <- sum(c(Cnidaria_cf, Cnidaria_eg))
      Ctenophora <- sum(c(Ctenophora_cf, Ctenophora_eg))
      Chordata <- sum(c(Chordata_cf, Chordata_eg))
      return(data.frame(depth=unique(y$depth), Cnidaria, Ctenophora, Chordata,
                        Cnidaria_cf, Ctenophora_cf, Chordata_cf,
                        Cnidaria_eg, Ctenophora_eg, Chordata_eg))
    }, .parallel=parallel, .progress="text")
    
    sd <- ddply(x, ~lat+lon, function(y){
      Cnidaria_cf <- sd(y[y$C %in% c(1000, 1100, 1200),"Cnidaria"])
      Cnidaria_eg <- NA
      Ctenophora_cf <- sd(y[y$C %in% c(800, 900, 1000),"Ctenophora"])
      Ctenophora_eg <- NA
      Chordata_cf <- sd(y[y$C %in% c(800, 900, 1000, 1100, 1200),"Chordata"])
      Chordata_eg <- sd(y[y$C %in% c(800, 900, 1000, 1100),"Chordata_eg"])
      Cnidaria <- Cnidaria_cf
      Ctenophora <- Ctenophora_cf
      Chordata <- sd(c(y[y$C %in% c(800, 900, 1000, 1100, 1200),"Chordata"], y[y$C %in% c(100, 800, 1200),"Chordata_eg"]))
      return(data.frame(depth=unique(y$depth), Cnidaria, Ctenophora, Chordata,
                        Cnidaria_cf, Ctenophora_cf, Chordata_cf,
                        Cnidaria_eg, Ctenophora_eg, Chordata_eg))
    }, .parallel=parallel, .progress="text")
    
    require(reshape2)
    mean <- melt(mean, id.vars = c("lat", "lon", "depth"), value.name = "mean", variable.name = "taxon")
    sd <- melt(sd, id.vars = c("lat", "lon", "depth"), value.name = "sd", variable.name = "taxon")
    mean_and_sd <- join(mean, sd, by=c("lat", "lon", "depth", "taxon"))
    # mean_and_sd <- cbind(mean, sd=sd[,5])
    
    mean_and_sd$taxon <- as.character(mean_and_sd$taxon)
    
    mean_and_sd <- mean_and_sd[!is.na(mean_and_sd$mean),]
    
    mean_and_sd$mean[which(mean_and_sd$mean==0)] <- NA
    mean_and_sd$sd[which(mean_and_sd$sd==0)] <- NA
    
    return(mean_and_sd)
  }
  
  d.export_mean_sd <- exportflux_mean_sd(d.export)
  seq.depth_mean_sd <- exportflux_mean_sd(seq.depth)
  seafloor_mean_sd <- exportflux_mean_sd(seafloor)
  
  # fill in the NAs with zeros (because those are areas of negative flux)
  d.export_mean_sd[which(is.na(d.export_mean_sd$mean)), "mean"] <- 0
  seq.depth_mean_sd[which(is.na(seq.depth_mean_sd$mean)), "mean"] <- 0
  seafloor_mean_sd[which(is.na(seafloor_mean_sd$mean)), "mean"] <- 0
  
  write.csv(d.export_mean_sd, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_100m_mean_sd.csv"), row.names=FALSE)
  write.csv(seq.depth_mean_sd, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_seq_depth_mean_sd.csv"), row.names=FALSE)
  write.csv(seafloor_mean_sd, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_seafloor_mean_sd.csv"), row.names=FALSE)
  
  # Regridding and plotting the means
  d.export_mean <- dcast(d.export_mean_sd[,1:5], lat+lon+depth~taxon, value.var = "mean")
  d.export_sd <- dcast(d.export_mean_sd[,c(1:4, 6)], lat+lon+depth~taxon, value.var="sd")
  seq.depth_mean <- dcast(seq.depth_mean_sd[,1:5], lat+lon+depth~taxon, value.var = "mean")
  seq.depth_sd <- dcast(seq.depth_mean_sd[,c(1:4, 6)], lat+lon+depth~taxon, value.var="sd")
  seafloor_mean <- dcast(seafloor_mean_sd[,1:5], lat+lon+depth~taxon, value.var="mean")
  seafloor_sd <- dcast(seafloor_mean_sd[,c(1:4, 6)], lat+lon+depth~taxon, value.var = "sd")
  
  
  if(PLOTTING==TRUE) {
    
    biovars <- c("Chordata", "Cnidaria", "Ctenophora") # these are the sum of the biomass carbon flux and the egestion flux
    # the taxa_cf and taxa_eg were combined in the earlier ddply loop
    
    d.export_mean <- join(regrid(d.export_mean, "depth"), regrid_geomean(d.export_mean, biovars), by=c("lat", "lon"))
    d.export_sd <- join(regrid(d.export_sd, "depth"), regrid_geomean(d.export_sd, biovars), by=c("lat", "lon"))
    seafloor_mean <- join(regrid(seafloor_mean, "depth"), regrid_geomean(seafloor_mean, biovars), by=c("lat", "lon"))
    seafloor_sd <- join(regrid(seafloor_sd, "depth"), regrid_geomean(seafloor_sd, biovars), by=c("lat", "lon"))
    
    # log and plot
    d.export_mean$logCnidaria <- log10(d.export_mean$Cnidaria + 1)
    d.export_mean$logChordata <- log10(d.export_mean$Chordata + 1)
    d.export_mean$logCtenophora <- log10(d.export_mean$Ctenophora + 1)
    
    p <- plotting(data=d.export_mean, var = "logCnidaria", plot.title="Biomass of Cnidaria past 100m (log-scale)", na.value=NA)
    ggsave(filename = str_c("plots/",CASE,"/",DIRNAME,"/ExportFlux_thermo_egcf_Cnidaria.pdf"), plot = p, width=16, height=9, units = "in")
    p <- plotting(data=d.export_mean, var = "logChordata", plot.title="Biomass of Chordata past 100m (log-scale)", na.value=NA)
    ggsave(filename = str_c("plots/",CASE,"/",DIRNAME,"/ExportFlux_thermo_egcf_Chordata.pdf"), plot = p, width=16, height=9, units = "in")
    p <- plotting(data=d.export_mean, var = "logCtenophora", plot.title="Biomass of Ctenophora at past 100m (log-scale)", na.value=NA)
    ggsave(filename = str_c("plots/",CASE,"/",DIRNAME,"/ExportFlux_thermo_egcf_Ctenophora.pdf"), plot = p, width=16, height=9, units = "in")
    
    
    seafloor_mean$logCnidaria <- log10(seafloor_mean$Cnidaria + 1)
    seafloor_mean$logChordata <- log10(seafloor_mean$Chordata + 1)
    seafloor_mean$logCtenophora <- log10(seafloor_mean$Ctenophora + 1)
    
    p <- plotting(data=seafloor_mean, var = "logCnidaria", plot.title="Biomass of Cnidaria at seafloor (log-scale)", na.value=NA)
    ggsave(filename = str_c("plots/",CASE,"/",DIRNAME,"/ExportFlux_seafloor_egcf_Cnidaria.pdf"), plot = p, width=16, height=9, units = "in")
    p <- plotting(data=seafloor_mean, var = "logChordata", plot.title="Biomass of Chordata at seafloor (log-scale)", na.value=NA)
    ggsave(filename = str_c("plots/",CASE,"/",DIRNAME,"/ExportFlux_seafloor_egcf_Chordata.pdf"), plot = p, width=16, height=9, units = "in")
    p <- plotting(data=seafloor_mean, var = "logCtenophora", plot.title="Biomass of Ctenophora at seafloor (log-scale)", na.value=NA)
    ggsave(filename = str_c("plots/",CASE,"/",DIRNAME,"/ExportFlux_seafloor_egcf_Ctenophora.pdf"), plot = p, width=16, height=9, units = "in")
    
  }
  
  
  # save regridded files
  write.csv(d.export_mean, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_100m_mean_rg.csv"), row.names=FALSE)
  write.csv(d.export_sd, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_100m_sd_rg.csv"), row.names=FALSE)
  write.csv(seq.depth_mean, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_seq_depth_mean_rg.csv"), row.names=FALSE)
  write.csv(seq.depth_sd, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_past_seq_depth_sd_rg.csv"), row.names=FALSE)
  write.csv(seafloor_mean, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_seafloor_mean_rg.csv"), row.names=FALSE)
  write.csv(seafloor_sd, str_c("data/",CASE,"/",DIRNAME,"/jelly_flux_seafloor_sd_rg.csv"), row.names=FALSE)
  
  
  
}



