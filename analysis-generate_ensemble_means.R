#!/usr/bin/env Rscript
#
#      Generate Ensemble Means for export fluxes
#
#-----------------------------------------------

# load libraries ----------------------------
library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")

source("functions.R")

for (C_ID in 1:3){
  #C_ID = 1
  
  CASES = c("0-baseline","1-high_biomass","2-low_biomass")
  CASE = CASES[C_ID]
  print(str_c('calculating ensemble means for case: ', CASE))
  
  #### UPPER OCEAN ####
  load(str_c("data/",CASE,"/mc_res/MC_ExportError_mc_uores.Rdata")) #variable name is mc_uores
  load(str_c("data/",CASE,"/mc_res/MC_ExportError_mc_uores_mat.Rdata")) #taxon levels: sum = 1, cnid = 2, cten = 3, chor = 4
  
  # double check the lat-lon's
  latlons <- apply(mc_uores_mat[,1:2,1],1:2,as.numeric)
  if(! all(latlons[,1]==mc_uores[[1]]$lat & latlons[,2]==mc_uores[[1]]$lon)){
    print("Warning: check dimensions")}

  # derive the ensemble mean
  rates <- apply(mc_uores_mat[,3:9,],1:3,as.numeric)
  ens_mean <- apply(rates,1:2, mean)
  
  ens_mean <- as.data.frame(ens_mean)
  ens_mean <- cbind(mc_uores[[1]][,c("lat","lon")], ens_mean, mc_uores[[1]][,c("taxon", "biome")])
  names(ens_mean) <- names(mc_uores[[1]])
  
  # recalculate sum (because took geometric mean)
  ens_mean <- ens_mean[ens_mean$taxon != "sum",]
  ens_sum <- ddply(ens_mean, ~lat+lon+biome, function(x){
    ingestion=sum(x$ingestion)
    respiration=sum(x$respiration)
    DOC=sum(x$DOC)
    reproduction=sum(x$reproduction)
    predation=sum(x$predation)
    flux=sum(x$flux)
    egestion=sum(x$egestion)
    return(data.frame(ingestion, respiration, DOC, reproduction, predation, flux, egestion))
  })
  ens_sum$taxon <- "sum"
  ens_mean <- rbind(ens_mean, ens_sum)
  
  save(ens_mean, file=str_c("data/",CASE,"/baseline/upperocean_ens_mean.Rdata"))
  
  #### EXPORT ####
  load(str_c("data/",CASE,"/mc_res/MC_ExportError_exprt.Rdata"))
  cn_lat = exprt_all$Cnidarians$lat
  cn_lon = exprt_all$Cnidarians$lon
  ct_lat = exprt_all$Ctenophores$lat
  ct_lon = exprt_all$Ctenophores$lon
  ch_lat = exprt_all$Chordates$lat
  ch_lon = exprt_all$Chordates$lon
  
  # extract variables
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
  
  
  
  ## # straight means
  mean_cn_cf_exprt = apply(cn_cf_exprt_all,1,mean)
  mean_cn_eg_exprt = apply(cn_eg_exprt_all,1,mean)
  mean_cn_exprt = mean_cn_cf_exprt + mean_cn_eg_exprt
  mean_ct_cf_exprt = apply(ct_cf_exprt_all,1,mean)
  mean_ct_eg_exprt = apply(ct_eg_exprt_all,1,mean)
  mean_ct_exprt = mean_ct_cf_exprt + mean_ct_eg_exprt
  mean_ch_cf_exprt = apply(ch_cf_exprt_all,1,mean)
  mean_ch_eg_exprt = apply(ch_eg_exprt_all,1,mean)
  mean_ch_exprt = mean_ch_cf_exprt + mean_ch_eg_exprt
  
  mean_cn_cf_seq_d = apply(cn_cf_seq_d_all,1,mean)
  mean_cn_eg_seq_d = apply(cn_eg_seq_d_all,1,mean)
  mean_cn_seq_d = mean_cn_cf_seq_d + mean_cn_eg_seq_d
  mean_ct_cf_seq_d = apply(ct_cf_seq_d_all,1,mean)
  mean_ct_eg_seq_d = apply(ct_eg_seq_d_all,1,mean)
  mean_ct_seq_d = mean_ct_cf_seq_d + mean_ct_eg_seq_d
  mean_ch_cf_seq_d = apply(ch_cf_seq_d_all,1,mean)
  mean_ch_eg_seq_d = apply(ch_eg_seq_d_all,1,mean)
  mean_ch_seq_d = mean_ch_cf_seq_d + mean_ch_eg_seq_d
  
  mean_cn_cf_seafl = apply(cn_cf_seafl_all,1,mean)
  mean_cn_eg_seafl = apply(cn_eg_seafl_all,1,mean)
  mean_cn_seafl = mean_cn_cf_seafl + mean_cn_eg_seafl
  mean_ct_cf_seafl = apply(ct_cf_seafl_all,1,mean)
  mean_ct_eg_seafl = apply(ct_eg_seafl_all,1,mean)
  mean_ct_seafl = mean_ct_cf_seafl + mean_ct_eg_seafl
  mean_ch_cf_seafl = apply(ch_cf_seafl_all,1,mean)
  mean_ch_eg_seafl = apply(ch_eg_seafl_all,1,mean)
  mean_ch_seafl = mean_ch_cf_seafl + mean_ch_eg_seafl
  
  
  export = rbind(data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria_cf", mean=mean_cn_cf_exprt),
                 data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria_eg", mean=mean_cn_eg_exprt),
                 data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria", mean=mean_cn_exprt),
                 data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora_cf", mean=mean_ct_cf_exprt),
                 data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora_eg", mean=mean_ct_eg_exprt),
                 data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora", mean=mean_ct_exprt),
                 data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata_cf", mean=mean_ch_cf_exprt),
                 data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata_eg", mean=mean_ch_eg_exprt),
                 data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata", mean=mean_ch_exprt))
  export$depth = "Export"
  
  seq_d = rbind(data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria_cf", mean=mean_cn_cf_seq_d),
                data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria_eg", mean=mean_cn_eg_seq_d),
                data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria", mean=mean_cn_seq_d),
                data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora_cf", mean=mean_ct_cf_seq_d),
                data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora_eg", mean=mean_ct_eg_seq_d),
                data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora", mean=mean_ct_seq_d),
                data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata_cf", mean=mean_ch_cf_seq_d),
                data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata_eg", mean=mean_ch_eg_seq_d),
                data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata", mean=mean_ch_seq_d))
  seq_d$depth = "Sequestration"
  
  seafloor = rbind(data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria_cf", mean=mean_cn_cf_seafl),
                   data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria_eg", mean=mean_cn_eg_seafl),
                   data.frame(lat=cn_lat, lon=cn_lon, taxon="Cnidaria", mean=mean_cn_seafl),
                   data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora_cf", mean=mean_ct_cf_seafl),
                   data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora_eg", mean=mean_ct_eg_seafl),
                   data.frame(lat=ct_lat, lon=ct_lon, taxon="Ctenophora", mean=mean_ct_seafl),
                   data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata_cf", mean=mean_ch_cf_seafl),
                   data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata_eg", mean=mean_ch_eg_seafl),
                   data.frame(lat=ch_lat, lon=ch_lon, taxon="Chordata", mean=mean_ch_seafl))
  seafloor$depth = "Seafloor"
  
  # save
  write.csv(export, paste0("data/",CASE,"/mc_res/jelly_flux_past_100m_ens_means.csv"), row.names=FALSE)
  write.csv(seq_d, paste0("data/",CASE,"/mc_res/jelly_flux_past_seq_depth_ens_means.csv"), row.names=FALSE)
  write.csv(seafloor, paste0("data/",CASE,"/mc_res/jelly_flux_seafloor_ens_means.csv"), row.names=FALSE)
  
  # in a different shape
  cnid_export = data.frame(lat=cn_lat, lon=cn_lon, taxon='cnid',
                           export=mean_cn_exprt, export_cf=mean_cn_cf_exprt, export_eg=mean_cn_eg_exprt, 
                           seq_d=mean_cn_seq_d, seq_d_cf=mean_cn_cf_seq_d, seq_d_eg=mean_cn_eg_seq_d,
                           seafloor=mean_cn_seafl, seafloor_cf=mean_cn_cf_seafl, seafloor_eg=mean_cn_eg_seafl)
  cten_export = data.frame(lat=ct_lat, lon=ct_lon, taxon='cten', 
                           export=mean_ct_exprt, export_cf=mean_ct_cf_exprt, export_eg=mean_ct_eg_exprt, 
                           seq_d=mean_ct_seq_d, seq_d_cf=mean_ct_cf_seq_d, seq_d_eg=mean_ct_eg_seq_d,
                           seafloor=mean_ct_seafl, seafloor_cf=mean_ct_cf_seafl, seafloor_eg=mean_ct_eg_seafl)
  chor_export = data.frame(lat=ch_lat, lon=ch_lon, taxon='chor',
                           export=mean_ch_exprt, export_cf=mean_ch_cf_exprt, export_eg=mean_ch_eg_exprt, 
                           seq_d=mean_ch_seq_d, seq_d_cf=mean_ch_cf_seq_d, seq_d_eg=mean_ch_eg_seq_d,
                           seafloor=mean_ch_seafl, seafloor_cf=mean_ch_cf_seafl, seafloor_eg=mean_ch_eg_seafl)
  
  # combine taxa
  allexport <- rbind(cnid_export, cten_export, chor_export)
  
  # Sum up all the results from all three groups per grid cell
  # fyi - this is slightly misleading because summing this way implies that the grid cells with no data are zeros
  # but we can do it for plotting
  sumres <- ddply(allexport, ~lat+lon, function(x){
    x <- x[, which(names(x) %ni% c("taxon", "lat", "lon"))] # get rid of the taxon column
    x <- colSums(x)
    return(x)
  }, .progress="text")
  
  allexport <- rbind(allexport, data.frame(sumres, taxon="sum"))
  
  save(allexport, file=str_c("data/",CASE,"/baseline/allexport.Rdata"))
  
  
}