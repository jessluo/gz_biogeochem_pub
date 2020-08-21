#!/usr/bin/env Rscript
#
#           Extract results from ensemble runs
#
#           Jessica Luo
#
#------------------------------------------------------------------------

### load libraries
source("functions.R")
library("plyr")
library("reshape2")

### set up case info
args <- commandArgs(trailingOnly=TRUE)
C_ID <- as.numeric(args[1])
#C_ID = 1
CASES = c("0-baseline","1-high_biomass","2-low_biomass")
CASE = CASES[C_ID]
dir.create(paste0("data/",CASE,"/mc_res/"), showWarnings = FALSE)

### set up constants and parameters
N_ENS <- 100

cnid_cf_sinking = c(1000, 1100, 1200)
cnid_eg_sinking = 100
cten_cf_sinking = c(800, 900, 1000)
cten_eg_sinking = 100
chor_cf_sinking = c(800, 900, 1000, 1100, 1200)
chor_eg_sinking = c(100, 800, 1200)

START <- 1
END <- N_ENS
print(paste("Extracting results from case:", CASE, "and running from", START, "to", END))

### load standard datasets (across all simulations)
rd <- read.csv(str_c("data/", CASE, "/jelly_biomass_1_deg_grid.csv"), stringsAsFactors=FALSE)
rd$latlon=str_c(rd$lat,',',rd$lon)

# load biomes
load('data/biomes_chl_mixedlayer.Rdata')

mc_uores <- list()
# load one and get dim
load(paste0("data/",CASE,"/MC_ExportError/1/results_MCsim_uo_res.Rdata"))
mc_uores_mat <- array(NA, dim=c(dim(uo_res), N_ENS))


### Ensemble results for Upper Ocean only
print("Compiling upper ocean results")
ptm <- proc.time()

# read data in
for(i in START:END){
  
  load(paste0("data/",CASE,"/MC_ExportError/",i,"/results_MCsim_uo_res.Rdata"))
  # save upper ocean results all together
  uo_res$taxon <- as.character(uo_res$taxon)
  mc_uores_mat[,,i] <- as.matrix(uo_res)
  
  # add biomes, save into list
  uo_res <- join(uo_res, biomes_df, by=c("lat","lon"))
  mc_uores[[i]] <- uo_res
  
}

print(proc.time()-ptm)


print("Compiling export model results")
ptm <- proc.time()

# initialize storage arrays - store them separately
# putting them all together just makes accessing them more difficult
cnid_cf = nrow(rd[rd$rank_phylum=="Cnidaria",])
cnid_eg = nrow(rd[rd$rank_phylum=="Cnidaria",])
cten_cf = nrow(rd[rd$rank_phylum=="Ctenophora",])
cten_eg = nrow(rd[rd$rank_phylum=="Ctenophora",])
chor_cf = nrow(rd[rd$rank_phylum=="Chordata",])
chor_eg = nrow(rd[rd$rank_phylum=="Chordata",])

cn_cf_seafl_all <- array(NA, dim=c(cnid_cf, length(cnid_cf_sinking), N_ENS))
cn_cf_exprt_all <- array(NA, dim=c(cnid_cf, length(cnid_cf_sinking), N_ENS))
cn_cf_seq_d_all <- array(NA, dim=c(cnid_cf, length(cnid_cf_sinking), N_ENS))
ct_cf_seafl_all <- array(NA, dim=c(cten_cf, length(cten_cf_sinking), N_ENS))
ct_cf_exprt_all <- array(NA, dim=c(cten_cf, length(cten_cf_sinking), N_ENS))
ct_cf_seq_d_all <- array(NA, dim=c(cten_cf, length(cten_cf_sinking), N_ENS))
ch_cf_seafl_all <- array(NA, dim=c(chor_cf, length(chor_cf_sinking), N_ENS))
ch_cf_exprt_all <- array(NA, dim=c(chor_cf, length(chor_cf_sinking), N_ENS))
ch_cf_seq_d_all <- array(NA, dim=c(chor_cf, length(chor_cf_sinking), N_ENS))

cn_eg_seafl_all <- array(NA, dim=c(cnid_eg, length(cnid_eg_sinking), N_ENS))
cn_eg_exprt_all <- array(NA, dim=c(cnid_eg, length(cnid_eg_sinking), N_ENS))
cn_eg_seq_d_all <- array(NA, dim=c(cnid_eg, length(cnid_eg_sinking), N_ENS))
ct_eg_seafl_all <- array(NA, dim=c(cten_eg, length(cten_eg_sinking), N_ENS))
ct_eg_exprt_all <- array(NA, dim=c(cten_eg, length(cten_eg_sinking), N_ENS))
ct_eg_seq_d_all <- array(NA, dim=c(cten_eg, length(cten_eg_sinking), N_ENS))
ch_eg_seafl_all <- array(NA, dim=c(chor_eg, length(chor_eg_sinking), N_ENS))
ch_eg_exprt_all <- array(NA, dim=c(chor_eg, length(chor_eg_sinking), N_ENS))
ch_eg_seq_d_all <- array(NA, dim=c(chor_eg, length(chor_eg_sinking), N_ENS))

# read data in
for(i in START:END){
  
  print(i)
  
  load(paste0("data/",CASE,"/MC_ExportError/",i,"/jelly_flux_reaching_seafloor.Rdata"))
  load(paste0("data/",CASE,"/MC_ExportError/",i,"/jelly_flux_past_100m.Rdata"))
  load(paste0("data/",CASE,"/MC_ExportError/",i,"/jelly_flux_past_seq_depth.Rdata"))

  d.export$latlon=str_c(d.export$lat,',',d.export$lon)
  seq.depth$latlon=str_c(seq.depth$lat,',',seq.depth$lon)
  seafloor$latlon=str_c(seafloor$lat,',',seafloor$lon)
  
  # this step should not be necessary... but JUST IN CASE, add the true zeros back in
  for (taxon in c("Cnidaria", "Ctenophora", "Chordata")){
    orig = rd[rd$rank_phylum==taxon,"latlon"]
    
    d.export[which(d.export$latlon %in% orig & is.na(d.export[,taxon])), taxon] <- 0
    d.export[which(d.export$latlon %in% orig & is.na(d.export[,str_c(taxon,'_eg')])), str_c(taxon,'_eg')] <- 0
    seq.depth[which(seq.depth$latlon %in% orig & is.na(seq.depth[,taxon])), taxon] <- 0
    seq.depth[which(seq.depth$latlon %in% orig & is.na(seq.depth[,str_c(taxon,'_eg')])), str_c(taxon,'_eg')] <- 0
    seafloor[which(seafloor$latlon %in% orig & is.na(seafloor[,taxon])), taxon] <- 0
    seafloor[which(seafloor$latlon %in% orig & is.na(seafloor[,str_c(taxon,'_eg')])), str_c(taxon,'_eg')] <- 0
    
  }


  # -- Cnidarians --
  cn_cf_seafl <- seafloor[which(seafloor$C %in% cnid_cf_sinking), c("lat", "lon", "C", "Cnidaria")]
  cn_cf_exprt <- d.export[which(d.export$C %in% cnid_cf_sinking), c("lat", "lon", "C", "Cnidaria")]
  cn_cf_seq_d <- seq.depth[which(seq.depth$C %in% cnid_cf_sinking), c("lat", "lon", "C", "Cnidaria")]

  cn_eg_seafl <- seafloor[which(seafloor$C %in% cnid_eg_sinking), c("lat", "lon", "C", "Cnidaria_eg")]
  cn_eg_exprt <- d.export[which(d.export$C %in% cnid_eg_sinking), c("lat", "lon", "C", "Cnidaria_eg")]
  cn_eg_seq_d <- seq.depth[which(seq.depth$C %in% cnid_eg_sinking), c("lat", "lon", "C", "Cnidaria_eg")]

  # remove all NAs
  cn_cf_seafl <- cn_cf_seafl[complete.cases(cn_cf_seafl),]
  cn_cf_exprt <- cn_cf_exprt[complete.cases(cn_cf_exprt),]
  cn_cf_seq_d <- cn_cf_seq_d[complete.cases(cn_cf_seq_d),]
  cn_eg_seafl <- cn_eg_seafl[complete.cases(cn_eg_seafl),]
  cn_eg_exprt <- cn_eg_exprt[complete.cases(cn_eg_exprt),]
  cn_eg_seq_d <- cn_eg_seq_d[complete.cases(cn_eg_seq_d),]
  
  # cast into wide format
  cn_cf_seafl <- dcast(cn_cf_seafl, formula = lat+lon ~ C, value.var="Cnidaria")
  cn_cf_exprt <- dcast(cn_cf_exprt, formula = lat+lon ~ C, value.var="Cnidaria")
  cn_cf_seq_d <- dcast(cn_cf_seq_d, formula = lat+lon ~ C, value.var="Cnidaria")

  cn_eg_seafl <- dcast(cn_eg_seafl, formula = lat+lon ~ C, value.var="Cnidaria_eg")
  cn_eg_exprt <- dcast(cn_eg_exprt, formula = lat+lon ~ C, value.var="Cnidaria_eg")
  cn_eg_seq_d <- dcast(cn_eg_seq_d, formula = lat+lon ~ C, value.var="Cnidaria_eg")

  # save into array
  cn_cf_seafl_all[,,i] <- as.matrix(cn_cf_seafl[,3:ncol(cn_cf_seafl)])
  cn_cf_exprt_all[,,i] <- as.matrix(cn_cf_exprt[,3:ncol(cn_cf_exprt)])
  cn_cf_seq_d_all[,,i] <- as.matrix(cn_cf_seq_d[,3:ncol(cn_cf_seq_d)])

  cn_eg_seafl_all[,,i] <- as.matrix(cn_eg_seafl[,3:ncol(cn_eg_seafl)])
  cn_eg_exprt_all[,,i] <- as.matrix(cn_eg_exprt[,3:ncol(cn_eg_exprt)])
  cn_eg_seq_d_all[,,i] <- as.matrix(cn_eg_seq_d[,3:ncol(cn_eg_seq_d)])


  # -- Ctenophores --
  ct_cf_seafl <- seafloor[which(seafloor$C %in% cten_cf_sinking), c("lat", "lon", "C", "Ctenophora")]
  ct_cf_exprt <- d.export[which(d.export$C %in% cten_cf_sinking), c("lat", "lon", "C", "Ctenophora")]
  ct_cf_seq_d <- seq.depth[which(seq.depth$C %in% cten_cf_sinking), c("lat", "lon", "C", "Ctenophora")]

  ct_eg_seafl <- seafloor[which(seafloor$C %in% cten_eg_sinking), c("lat", "lon", "C", "Ctenophora_eg")]
  ct_eg_exprt <- d.export[which(d.export$C %in% cten_eg_sinking), c("lat", "lon", "C", "Ctenophora_eg")]
  ct_eg_seq_d <- seq.depth[which(seq.depth$C %in% cten_eg_sinking), c("lat", "lon", "C", "Ctenophora_eg")]

  # remove all NAs
  ct_cf_seafl <- ct_cf_seafl[complete.cases(ct_cf_seafl),]
  ct_cf_exprt <- ct_cf_exprt[complete.cases(ct_cf_exprt),]
  ct_cf_seq_d <- ct_cf_seq_d[complete.cases(ct_cf_seq_d),]
  ct_eg_seafl <- ct_eg_seafl[complete.cases(ct_eg_seafl),]
  ct_eg_exprt <- ct_eg_exprt[complete.cases(ct_eg_exprt),]
  ct_eg_seq_d <- ct_eg_seq_d[complete.cases(ct_eg_seq_d),]

  # cast into wide format
  ct_cf_seafl <- dcast(ct_cf_seafl, formula = lat+lon ~ C, value.var="Ctenophora")
  ct_cf_exprt <- dcast(ct_cf_exprt, formula = lat+lon ~ C, value.var="Ctenophora")
  ct_cf_seq_d <- dcast(ct_cf_seq_d, formula = lat+lon ~ C, value.var="Ctenophora")

  ct_eg_seafl <- dcast(ct_eg_seafl, formula = lat+lon ~ C, value.var="Ctenophora_eg")
  ct_eg_exprt <- dcast(ct_eg_exprt, formula = lat+lon ~ C, value.var="Ctenophora_eg")
  ct_eg_seq_d <- dcast(ct_eg_seq_d, formula = lat+lon ~ C, value.var="Ctenophora_eg")

  # save into array
  ct_cf_seafl_all[,,i] <- as.matrix(ct_cf_seafl[,3:ncol(ct_cf_seafl)])
  ct_cf_exprt_all[,,i] <- as.matrix(ct_cf_exprt[,3:ncol(ct_cf_exprt)])
  ct_cf_seq_d_all[,,i] <- as.matrix(ct_cf_seq_d[,3:ncol(ct_cf_seq_d)])

  ct_eg_seafl_all[,,i] <- as.matrix(ct_eg_seafl[,3:ncol(ct_eg_seafl)])
  ct_eg_exprt_all[,,i] <- as.matrix(ct_eg_exprt[,3:ncol(ct_eg_exprt)])
  ct_eg_seq_d_all[,,i] <- as.matrix(ct_eg_seq_d[,3:ncol(ct_eg_seq_d)])


  # -- Tunicates --
  ch_cf_seafl <- seafloor[which(seafloor$C %in% chor_cf_sinking), c("lat", "lon", "C", "Chordata")]
  ch_cf_exprt <- d.export[which(d.export$C %in% chor_cf_sinking), c("lat", "lon", "C", "Chordata")]
  ch_cf_seq_d <- seq.depth[which(seq.depth$C %in% chor_cf_sinking), c("lat", "lon", "C", "Chordata")]

  ch_eg_seafl <- seafloor[which(seafloor$C %in% chor_eg_sinking), c("lat", "lon", "C", "Chordata_eg")]
  ch_eg_exprt <- d.export[which(d.export$C %in% chor_eg_sinking), c("lat", "lon", "C", "Chordata_eg")]
  ch_eg_seq_d <- seq.depth[which(seq.depth$C %in% chor_eg_sinking), c("lat", "lon", "C", "Chordata_eg")]

  # remove all NAs
  ch_cf_seafl <- ch_cf_seafl[complete.cases(ch_cf_seafl),]
  ch_cf_exprt <- ch_cf_exprt[complete.cases(ch_cf_exprt),]
  ch_cf_seq_d <- ch_cf_seq_d[complete.cases(ch_cf_seq_d),]
  ch_eg_seafl <- ch_eg_seafl[complete.cases(ch_eg_seafl),]
  ch_eg_exprt <- ch_eg_exprt[complete.cases(ch_eg_exprt),]
  ch_eg_seq_d <- ch_eg_seq_d[complete.cases(ch_eg_seq_d),]
  
  # cast into wide format
  ch_cf_seafl <- dcast(ch_cf_seafl, formula = lat+lon ~ C, value.var="Chordata")
  ch_cf_exprt <- dcast(ch_cf_exprt, formula = lat+lon ~ C, value.var="Chordata")
  ch_cf_seq_d <- dcast(ch_cf_seq_d, formula = lat+lon ~ C, value.var="Chordata")

  ch_eg_seafl <- dcast(ch_eg_seafl, formula = lat+lon ~ C, value.var="Chordata_eg")
  ch_eg_exprt <- dcast(ch_eg_exprt, formula = lat+lon ~ C, value.var="Chordata_eg")
  ch_eg_seq_d <- dcast(ch_eg_seq_d, formula = lat+lon ~ C, value.var="Chordata_eg")

  # save into array
  ch_cf_seafl_all[,,i] <- as.matrix(ch_cf_seafl[,3:ncol(ch_cf_seafl)])
  ch_cf_exprt_all[,,i] <- as.matrix(ch_cf_exprt[,3:ncol(ch_cf_exprt)])
  ch_cf_seq_d_all[,,i] <- as.matrix(ch_cf_seq_d[,3:ncol(ch_cf_seq_d)])

  ch_eg_seafl_all[,,i] <- as.matrix(ch_eg_seafl[,3:ncol(ch_eg_seafl)])
  ch_eg_exprt_all[,,i] <- as.matrix(ch_eg_exprt[,3:ncol(ch_eg_exprt)])
  ch_eg_seq_d_all[,,i] <- as.matrix(ch_eg_seq_d[,3:ncol(ch_eg_seq_d)])



}

# 
save(cn_cf_seafl_all, file=paste0("data/",CASE,"/mc_res/cn_cf_seafl_all.Rdata"))
save(cn_eg_seafl_all, file=paste0("data/",CASE,"/mc_res/cn_eg_seafl_all.Rdata"))
save(cn_cf_exprt_all, file=paste0("data/",CASE,"/mc_res/cn_cf_exprt_all.Rdata"))
save(cn_eg_exprt_all, file=paste0("data/",CASE,"/mc_res/cn_eg_exprt_all.Rdata"))
save(cn_cf_seq_d_all, file=paste0("data/",CASE,"/mc_res/cn_cf_seq_d_all.Rdata"))
save(cn_eg_seq_d_all, file=paste0("data/",CASE,"/mc_res/cn_eg_seq_d_all.Rdata"))

save(ct_cf_seafl_all, file=paste0("data/",CASE,"/mc_res/ct_cf_seafl_all.Rdata"))
save(ct_eg_seafl_all, file=paste0("data/",CASE,"/mc_res/ct_eg_seafl_all.Rdata"))
save(ct_cf_exprt_all, file=paste0("data/",CASE,"/mc_res/ct_cf_exprt_all.Rdata"))
save(ct_eg_exprt_all, file=paste0("data/",CASE,"/mc_res/ct_eg_exprt_all.Rdata"))
save(ct_cf_seq_d_all, file=paste0("data/",CASE,"/mc_res/ct_cf_seq_d_all.Rdata"))
save(ct_eg_seq_d_all, file=paste0("data/",CASE,"/mc_res/ct_eg_seq_d_all.Rdata"))

save(ch_cf_seafl_all, file=paste0("data/",CASE,"/mc_res/ch_cf_seafl_all.Rdata"))
save(ch_eg_seafl_all, file=paste0("data/",CASE,"/mc_res/ch_eg_seafl_all.Rdata"))
save(ch_cf_exprt_all, file=paste0("data/",CASE,"/mc_res/ch_cf_exprt_all.Rdata"))
save(ch_eg_exprt_all, file=paste0("data/",CASE,"/mc_res/ch_eg_exprt_all.Rdata"))
save(ch_cf_seq_d_all, file=paste0("data/",CASE,"/mc_res/ch_cf_seq_d_all.Rdata"))
save(ch_eg_seq_d_all, file=paste0("data/",CASE,"/mc_res/ch_eg_seq_d_all.Rdata"))

print(proc.time() - ptm)

#----------------------------------------------------
# set dim names for all files

dimnames(cn_cf_seq_d_all) <- dimnames(cn_cf_exprt_all) <- dimnames(cn_cf_seafl_all) <- list(latlon=NULL, sinking_rates=cnid_cf_sinking, runNo=NULL)
dimnames(cn_eg_seq_d_all) <- dimnames(cn_eg_exprt_all) <- dimnames(cn_eg_seafl_all) <- list(latlon=NULL, sinking_rates=cnid_eg_sinking, runNo=NULL)

dimnames(ct_cf_seq_d_all) <- dimnames(ct_cf_exprt_all) <- dimnames(ct_cf_seafl_all) <- list(latlon=NULL, sinking_rates=cten_cf_sinking, runNo=NULL)
dimnames(ct_eg_seq_d_all) <- dimnames(ct_eg_exprt_all) <- dimnames(ct_eg_seafl_all) <- list(latlon=NULL, sinking_rates=cten_eg_sinking, runNo=NULL)

dimnames(ch_cf_seq_d_all) <- dimnames(ch_cf_exprt_all) <- dimnames(ch_cf_seafl_all) <- list(latlon=NULL, sinking_rates=chor_cf_sinking, runNo=NULL)
dimnames(ch_eg_seq_d_all) <- dimnames(ch_eg_exprt_all) <- dimnames(ch_eg_seafl_all) <- list(latlon=NULL, sinking_rates=chor_eg_sinking, runNo=NULL)


# combine all of them into lists
seafl_all <- list(Cnidarians=list(biomass=cn_cf_seafl_all, egestion=cn_eg_seafl_all, lat=cn_cf_seafl$lat, lon=cn_cf_seafl$lon),
                  Ctenophores=list(biomass=ct_cf_seafl_all, egestion=ct_eg_seafl_all, lat=ct_cf_seafl$lat, lon=ct_cf_seafl$lon),
                  Chordates=list(biomass=ch_cf_seafl_all, egestion=ch_eg_seafl_all, lat=ch_cf_seafl$lat, lon=ch_cf_seafl$lon))

exprt_all <- list(Cnidarians=list(biomass=cn_cf_exprt_all, egestion=cn_eg_exprt_all, lat=cn_cf_exprt$lat, lon=cn_cf_exprt$lon),
                  Ctenophores=list(biomass=ct_cf_exprt_all, egestion=ct_eg_exprt_all, lat=ct_cf_exprt$lat, lon=ct_cf_exprt$lon),
                  Chordates=list(biomass=ch_cf_exprt_all, egestion=ch_eg_exprt_all, lat=ch_cf_exprt$lat, lon=ch_cf_exprt$lon))

seq_d_all <- list(Cnidarians=list(biomass=cn_cf_seq_d_all, egestion=cn_eg_seq_d_all, lat=cn_cf_seq_d$lat, lon=cn_cf_seq_d$lon),
                  Ctenophores=list(biomass=ct_cf_seq_d_all, egestion=ct_eg_seq_d_all, lat=ct_cf_seq_d$lat, lon=ct_cf_seq_d$lon),
                  Chordates=list(biomass=ch_cf_seq_d_all, egestion=ch_eg_seq_d_all, lat=ch_cf_seq_d$lat, lon=ch_cf_seq_d$lon))

# save files
save(seafl_all, file=paste0("data/",CASE,"/mc_res/MC_ExportError_seafl.Rdata"))
save(exprt_all, file=paste0("data/",CASE,"/mc_res/MC_ExportError_exprt.Rdata"))
save(seq_d_all, file=paste0("data/",CASE,"/mc_res/MC_ExportError_seq_d.Rdata"))

save(mc_uores, file=paste0("data/",CASE,"/mc_res/MC_ExportError_mc_uores.Rdata"))
save(mc_uores_mat, file=paste0("data/",CASE,"/mc_res/MC_ExportError_mc_uores_mat.Rdata"))

