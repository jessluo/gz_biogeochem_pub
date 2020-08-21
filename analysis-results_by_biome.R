#!/usr/bin/env Rscript
#
#     Uses a biome mask to calculate
#       global means
#
#
# -------------------------------------------

# load libraries ----------------------------
library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")
library("abind")

# set parameters ------------------------------
args <- commandArgs(trailingOnly=TRUE)
C_ID <- as.numeric(args[1])
# C_ID = 1
CASES = c("0-baseline","1-high_biomass","2-low_biomass")
CASE = CASES[C_ID]

N_ENS = 100

# constants
OCEANS_SURFACE_AREA <- 3.619e+14
OCEANS_SURFACE_VOL <- 6.708e+16

saveFiles=TRUE

# load functions ------------------------------
source("functions.R")

calculate_biome_mean_UO <- function(data, taxon, biomes=biomes_df, n_biome=N_BIOME, n_ens=N_ENS, saveFiles=TRUE){
  
  mc_biomeres_mat <- array(NA, dim=c(n_biome, 8, n_ens))
  
  #print('Extract and compute averages for each biome')
  for (i in 1:length(data)){
    temp <- data[[i]]
    temp <- temp[which(temp$taxon==taxon),] # swap out cnid, cten, chor if you want those specific fluxes
    temp <- temp[,which(names(temp) != "taxon")]
    temp$biome = as.numeric(temp$biome)
    # calculate biome means first
    temp_b <- ddply(temp, ~biome, function(x){
      return(apply(x[,which(names(temp) %ni% c("lat","lon","biome"))], MARGIN = 2, mean))
    })
    mc_biomeres_mat[,,i] <- as.matrix(temp_b)
  }
  
  # dimensions of mc_biomeres_mat are: [biome, rates, ensemble]
  
  # Then compute min, max, and mean for each biome --
  mins <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
  maxs <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
  uo_mean <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
  bcis_upper <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
  bcis_lower <- array(NA, dim=dim(mc_biomeres_mat)[1:2])
  
  for (i in 1:dim(mc_biomeres_mat)[1]){ # loop over every biome
    temp <- mc_biomeres_mat[i,,]
    temp <- t(temp) # transpose
    
    rates <- temp[,2:8]

    min <- apply(rates, 2, min, na.omit=TRUE)
    min <- c(i, min)
    mins[i,] <- min
    
    max <- apply(rates, 2, max, na.omit=TRUE)
    max <- c(i, max)
    maxs[i,] <- max
    
    avg <- apply(rates, 2, mean, na.rm=TRUE) 
    avg <- c(i, avg)
    uo_mean[i,] <- avg
        
    bci <- apply(rates, 2, return_bci)
    bci_lower <- c(i, bci[1,])
    bci_upper <- c(i, bci[2,])
    bcis_lower[i,] <- bci_lower
    bcis_upper[i,] <- bci_upper
    
  }
  # if(saveFiles){
  #     save(bcis_lower, file=str_c("data/",CASE,"/mc_res/Biomes_UO_",taxon,"_QUANT_lower.Rdata"))
  #     save(bcis_upper, file=str_c("data/",CASE,"/mc_res/Biomes_UO_",taxon,"_QUANT_upper.Rdata"))
  # }
  
  # Put results together --
  res <- cbind(melt(as.data.frame(uo_mean), id.vars = "V1", value.name = "mean"),
               melt(as.data.frame(bcis_lower), id.vars = "V1", value.name = "lower")[,"lower"],
               melt(as.data.frame(bcis_upper), id.vars = "V1", value.name = "upper")[,"upper"])
  names(res) <- c("biome","var","mean","lower","upper")
  res$var <- factor(res$var, labels=c("ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion"))
  res$biome <- cut(res$biome, breaks=n_biome, labels=c("COAST","LC","HCSS","HCPS"))
  
  res <- join(res, av_biomesum, by="biome")
  res$mean_tot <- res$mean * res$vol
  res$lower_tot <- res$lower * res$vol
  res$upper_tot <- res$upper * res$vol
  
  if(saveFiles) write.csv(res, str_c("data/",CASE,"/mc_res/Biomes_UO_",taxon,"_MeanResults.csv"), row.names=FALSE)
  
  return(res)
  
}

sum_by_biome <- function(data, lat, lon, biomes=biomes_df, avb=av_biomesum, av_ind=av, exclude_shallow=FALSE){
  d = data.frame(lat, lon, data)
  d = join(d, biomes, by=c("lat","lon"))
  
  if(exclude_shallow){
    d <- join(d, av_ind, by=c("lat","lon", "biome"))
    d <- d[which(d$depth.m > 50), which(names(d) %ni% c("area.km2", "area.m2", "vol", "multiplier", "depth.m"))]
  }
  
  dm = melt(d, id.vars = c("lat","lon","biome"))
  dm = ddply(dm, ~biome, function(x){ return(data.frame(mean=mean(x$value))) })
  dm = join(dm, avb, by="biome")
  dm$tot = dm$mean * dm$vol
  return(dm[,c("biome","mean","tot")])
}



conf_int_by_biome <- function(data, lat, lon, biomes=biomes_df, avb=av_biomesum, av_ind=av, exclude_shallow=FALSE){
  d = data.frame(lat, lon, data)
  d = join(d, biomes, by=c("lat","lon"))
  
  if(exclude_shallow){
    d <- join(d, av_ind, by=c("lat","lon", "biome"))
    d <- d[which(d$depth.m > 50), which(names(d) %ni% c("area.km2", "area.m2", "vol", "multiplier", "depth.m"))]
  }
  
  # first calculate a biome mean
  d = ddply(d, ~biome, function(x){
    return(colMeans(x[,str_c("X",c(1:100))]))
  })
  
  # then calculate the confidence intervals
  dm = melt(d, id.vars = c("biome"))
  ci = ddply(dm, ~biome, function(x){ 
    return_bci(x$value) 
    })
  return(ci)
}


# Import standard files------------------------
# load biomes
load('data/biomes_chl_mixedlayer.Rdata')
N_BIOME=4#length(unique(biomes_df$biome[!is.na(biomes_df$biome)]))

# import coastline
coastline.world <- read.csv("data/gshhg_world_c.csv")
load("data/coastline.grid_1deg.Rdata") # 1 deg grid

load("data/1-DEGREE_LATLONGRID_SEAS.Rdata")
seas <- join(seas, read.csv("raw_data/oceans.csv", as.is=TRUE))

# import area/volume
av <- read.csv("data/surfaceArea_volume_1-deg.csv")
av <- join(av,biomes_df, by=c("lat","lon"))
av$depth.m = av$vol / av$area.m2
av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(vol=sum(x$vol)))})
names(av_biomesum) <- c("biome", "vol")
av_biomedepth <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(depth.m=mean(x$depth.m)))})
#sum(av_biomesum$vol)


# Upper Ocean Ensemble results ----------
print(str_c('processing case: ',CASE))

# all units are g C m^-3 y^-1
load(str_c("data/",CASE,"/mc_res/MC_ExportError_mc_uores.Rdata")) #variable name is mc_uores

res <- data.frame()
for (taxon in c("cnid","cten","chor")){
  tres <- calculate_biome_mean_UO(mc_uores, taxon=taxon, saveFiles=saveFiles)
  tres$taxon <- taxon
  res <- rbind(res, tres)
  print(str_c(CASE, ", ", taxon)); print(ddply(tres, ~var, function(x){return(sum(x$mean_tot))}))
}

sum_res <- ddply(res, ~biome+var, function(x){
  d = colSums(x[,c("mean", "lower", "upper", "vol", "mean_tot", "lower_tot", "upper_tot")])
  d = as.data.frame(t(d))
  d$vol = x$vol[1]
  d$taxon = "sum"
  return(d)
})

print(str_c(CASE, ", sum")); print(ddply(sum_res, ~var, function(x){return(sum(x$mean_tot))}))
write.csv(res, str_c("data/",CASE,"/mc_res/Biomes_UO_sum_MeanResults.csv"), row.names=FALSE)


# Export Ensemble results ----------
load(str_c("data/",CASE,"/mc_res/MC_ExportError_exprt.Rdata"))
load(str_c("data/",CASE,"/mc_res/MC_ExportError_seq_d.Rdata"))
load(str_c("data/",CASE,"/mc_res/MC_ExportError_seafl.Rdata"))

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

bci_names = c("slow_lo", "slow_hi", "fast_lo", "fast_hi")

# A. EXPORT DEPTH 100 m --
res <- data.frame()

# 1. extract mean results
exprt_mean_cn <- apply(cn_cf_exprt_all,c(1,3),mean) + apply(cn_eg_exprt_all,c(1,3),mean)
exprt_mean_ct <- apply(ct_cf_exprt_all,c(1,3),mean) + apply(ct_eg_exprt_all,c(1,3),mean)
exprt_mean_ch <- apply(ch_cf_exprt_all,c(1,3),mean) + apply(ch_eg_exprt_all,c(1,3),mean)

cnid=sum_by_biome(exprt_mean_cn, exprt_all$Cnidarians$lat, exprt_all$Cnidarians$lon)
cten=sum_by_biome(exprt_mean_ct, exprt_all$Ctenophores$lat, exprt_all$Ctenophores$lon)
chor=sum_by_biome(exprt_mean_ch, exprt_all$Chordates$lat, exprt_all$Chordates$lon)

exprt_mean_res <- join(data.frame(biome=cnid$biome, cnid[,2:3]+cten[,2:3]+chor[,2:3]), av_biomesum, by="biome")
names(exprt_mean_res) <- c("biome", "Export", "Export_total", "vol")

save(exprt_mean_res, file=str_c("data/",CASE,"/mc_res/Biomes_export_depth_meanResults.Rdata"))

# 2. extract values based on the slowest vs fastest sinking speeds
exprt_slow_cn <- cn_cf_exprt_all[,1,] + cn_eg_exprt_all[,1,]
exprt_slow_ct <- ct_cf_exprt_all[,1,] + ct_eg_exprt_all[,1,]
exprt_slow_ch <- ch_cf_exprt_all[,1,] + ch_eg_exprt_all[,1,]
exprt_fast_cn <- cn_cf_exprt_all[,dim(cn_cf_exprt_all)[2],] + cn_eg_exprt_all[,dim(cn_eg_exprt_all)[2],]
exprt_fast_ct <- ct_cf_exprt_all[,dim(ct_cf_exprt_all)[2],] + ct_eg_exprt_all[,dim(ct_eg_exprt_all)[2],]
exprt_fast_ch <- ch_cf_exprt_all[,dim(ch_cf_exprt_all)[2],] + ch_eg_exprt_all[,dim(ch_eg_exprt_all)[2],]

cnid_slow=conf_int_by_biome(exprt_slow_cn, exprt_all$Cnidarians$lat, exprt_all$Cnidarians$lon)
cnid_fast=conf_int_by_biome(exprt_fast_cn, exprt_all$Cnidarians$lat, exprt_all$Cnidarians$lon)
cten_slow=conf_int_by_biome(exprt_slow_ct, exprt_all$Ctenophores$lat, exprt_all$Ctenophores$lon)
cten_fast=conf_int_by_biome(exprt_fast_ct, exprt_all$Ctenophores$lat, exprt_all$Ctenophores$lon)
chor_slow=conf_int_by_biome(exprt_slow_ch, exprt_all$Chordates$lat, exprt_all$Chordates$lon)
chor_fast=conf_int_by_biome(exprt_fast_ch, exprt_all$Chordates$lat, exprt_all$Chordates$lon)

slow = cnid_slow[,2:3]+cten_slow[,2:3]+chor_slow[,2:3]
fast = cnid_fast[,2:3]+cten_fast[,2:3]+chor_fast[,2:3]

bci_exprt = data.frame(biome=cnid_slow$biome, slow, fast)
names(bci_exprt) = c("biome", bci_names)
bci_exprt <- join(bci_exprt, av_biomesum, by="biome")
bci_exprt <- cbind(bci_exprt, bci_exprt[,bci_names] * bci_exprt$vol)
names(bci_exprt)[7:10]<-str_c(bci_names,"_total")

# save(bci_exprt, file=str_c("data/",CASE,"/mc_res/Biomes_export_depth_slow_v_fast_quantile.Rdata"))
# load(str_c("data/",CASE,"/mc_res/Biomes_export_depth_slow_v_fast_quantile.Rdata"))

# 3. Construct a single table with all the sinking flux means and ranges
cnid_range = melt(rbind(data.frame(cnid_slow, speed="slow"), data.frame(cnid_fast, speed="fast")), id.vars = c("biome","speed"))
cten_range = melt(rbind(data.frame(cten_slow, speed="slow"), data.frame(cten_fast, speed="fast")), id.vars = c("biome","speed"))
chor_range = melt(rbind(data.frame(chor_slow, speed="slow"), data.frame(chor_fast, speed="fast")), id.vars = c("biome","speed"))

cnid_all = rbind(cnid_range, data.frame(melt(cnid[,c("biome","mean")], id.vars = "biome"), speed="all"))
cten_all = rbind(cten_range, data.frame(melt(cten[,c("biome","mean")], id.vars = "biome"), speed="all"))
chor_all = rbind(chor_range, data.frame(melt(chor[,c("biome","mean")], id.vars = "biome"), speed="all"))

exprt_all = rbind(data.frame(cnid_all, taxon="cnid"), data.frame(cten_all, taxon="cten"), data.frame(chor_all, taxon="chor"))
names(exprt_all) = c("biome","speed","value","flux","taxon")
exprt_all$value = factor(exprt_all$value, labels=c("low","high","mean"))

# 4. Extract global sums for slow vs. fast sinkers confidence intervals
res <- rbind(res, data.frame(depth="Export", speed="slow", t(unname(colSums(bci_exprt[c("slow_lo_total","slow_hi_total")])))))
res <- rbind(res, data.frame(depth="Export", speed="fast", t(unname(colSums(bci_exprt[c("fast_lo_total","fast_hi_total")])))))



# B. SEQUESTRATION DEPTH --
# 1. extract mean results
seq_d_mean_cn <- apply(cn_cf_seq_d_all,c(1,3),mean) + apply(cn_eg_seq_d_all,c(1,3),mean)
seq_d_mean_ct <- apply(ct_cf_seq_d_all,c(1,3),mean) + apply(ct_eg_seq_d_all,c(1,3),mean)
seq_d_mean_ch <- apply(ch_cf_seq_d_all,c(1,3),mean) + apply(ch_eg_seq_d_all,c(1,3),mean)

cnid=sum_by_biome(seq_d_mean_cn, seq_d_all$Cnidarians$lat, seq_d_all$Cnidarians$lon)
cten=sum_by_biome(seq_d_mean_ct, seq_d_all$Ctenophores$lat, seq_d_all$Ctenophores$lon)
chor=sum_by_biome(seq_d_mean_ch, seq_d_all$Chordates$lat, seq_d_all$Chordates$lon)

seq_d_mean_res <- join(data.frame(biome=cnid$biome, cnid[,2:3]+cten[,2:3]+chor[,2:3]), av_biomesum, by="biome")
names(seq_d_mean_res) <- c("biome", "Sequestration", "Sequestration_total", "vol")

save(seq_d_mean_res, file=str_c("data/",CASE,"/mc_res/Biomes_seq_depth_meanResults.Rdata"))

# 2. extract values based on the slowest vs fastest sinking speeds
seq_d_slow_cn <- cn_cf_seq_d_all[,1,] + cn_eg_seq_d_all[,1,]
seq_d_slow_ct <- ct_cf_seq_d_all[,1,] + ct_eg_seq_d_all[,1,]
seq_d_slow_ch <- ch_cf_seq_d_all[,1,] + ch_eg_seq_d_all[,1,]
seq_d_fast_cn <- cn_cf_seq_d_all[,dim(cn_cf_seq_d_all)[2],] + cn_eg_seq_d_all[,dim(cn_eg_seq_d_all)[2],]
seq_d_fast_ct <- ct_cf_seq_d_all[,dim(ct_cf_seq_d_all)[2],] + ct_eg_seq_d_all[,dim(ct_eg_seq_d_all)[2],]
seq_d_fast_ch <- ch_cf_seq_d_all[,dim(ch_cf_seq_d_all)[2],] + ch_eg_seq_d_all[,dim(ch_eg_seq_d_all)[2],]

cnid_slow=conf_int_by_biome(seq_d_slow_cn, seq_d_all$Cnidarians$lat, seq_d_all$Cnidarians$lon)
cnid_fast=conf_int_by_biome(seq_d_fast_cn, seq_d_all$Cnidarians$lat, seq_d_all$Cnidarians$lon)
cten_slow=conf_int_by_biome(seq_d_slow_ct, seq_d_all$Ctenophores$lat, seq_d_all$Ctenophores$lon)
cten_fast=conf_int_by_biome(seq_d_fast_ct, seq_d_all$Ctenophores$lat, seq_d_all$Ctenophores$lon)
chor_slow=conf_int_by_biome(seq_d_slow_ch, seq_d_all$Chordates$lat, seq_d_all$Chordates$lon)
chor_fast=conf_int_by_biome(seq_d_fast_ch, seq_d_all$Chordates$lat, seq_d_all$Chordates$lon)

slow = cnid_slow[,2:3]+cten_slow[,2:3]+chor_slow[,2:3]
fast = cnid_fast[,2:3]+cten_fast[,2:3]+chor_fast[,2:3]

bci_seq_d = data.frame(biome=cnid_slow$biome, slow, fast)
names(bci_seq_d) = c("biome", bci_names)
bci_seq_d <- join(bci_seq_d, av_biomesum, by="biome")
bci_seq_d <- cbind(bci_seq_d, bci_seq_d[,bci_names] * bci_seq_d$vol)
names(bci_seq_d)[7:10]<-str_c(bci_names,"_total")

# save(bci_seq_d, file=str_c("data/",CASE,"/mc_res/Biomes_seq_depth_slow_v_fast_quantile.Rdata"))
# load(str_c("data/",CASE,"/mc_res/Biomes_seq_depth_slow_v_fast_quantile.Rdata"))

# 3. Construct a single table with all the sinking flux means and ranges
cnid_range = melt(rbind(data.frame(cnid_slow, speed="slow"), data.frame(cnid_fast, speed="fast")), id.vars = c("biome","speed"))
cten_range = melt(rbind(data.frame(cten_slow, speed="slow"), data.frame(cten_fast, speed="fast")), id.vars = c("biome","speed"))
chor_range = melt(rbind(data.frame(chor_slow, speed="slow"), data.frame(chor_fast, speed="fast")), id.vars = c("biome","speed"))

cnid_all = rbind(cnid_range, data.frame(melt(cnid[,c("biome","mean")], id.vars = "biome"), speed="all"))
cten_all = rbind(cten_range, data.frame(melt(cten[,c("biome","mean")], id.vars = "biome"), speed="all"))
chor_all = rbind(chor_range, data.frame(melt(chor[,c("biome","mean")], id.vars = "biome"), speed="all"))

seq_d_all = rbind(data.frame(cnid_all, taxon="cnid"), data.frame(cten_all, taxon="cten"), data.frame(chor_all, taxon="chor"))
names(seq_d_all) = c("biome","speed","value","flux","taxon")
seq_d_all$value = factor(seq_d_all$value, labels=c("low","high","mean"))

# 4. Extract global sums for slow vs. fast sinkers confidence intervals
res <- rbind(res, data.frame(depth="Sequestration", speed="slow", t(unname(colSums(bci_seq_d[c("slow_lo_total","slow_hi_total")])))))
res <- rbind(res, data.frame(depth="Sequestration", speed="fast", t(unname(colSums(bci_seq_d[c("fast_lo_total","fast_hi_total")])))))



# C. SEAFLOOR --
# 1. extract mean results
seafl_mean_cn <- apply(cn_cf_seafl_all,c(1,3),mean) + apply(cn_eg_seafl_all,c(1,3),mean)
seafl_mean_ct <- apply(ct_cf_seafl_all,c(1,3),mean) + apply(ct_eg_seafl_all,c(1,3),mean)
seafl_mean_ch <- apply(ch_cf_seafl_all,c(1,3),mean) + apply(ch_eg_seafl_all,c(1,3),mean)

cnid=sum_by_biome(seafl_mean_cn, seafl_all$Cnidarians$lat, seafl_all$Cnidarians$lon)
cten=sum_by_biome(seafl_mean_ct, seafl_all$Ctenophores$lat, seafl_all$Ctenophores$lon)
chor=sum_by_biome(seafl_mean_ch, seafl_all$Chordates$lat, seafl_all$Chordates$lon)

seafl_mean_res <- join(data.frame(biome=cnid$biome, cnid[,2:3]+cten[,2:3]+chor[,2:3]), av_biomesum, by="biome")
names(seafl_mean_res) <- c("biome", "Seafloor", "Seafloor_total", "vol")

save(seafl_mean_res, file=str_c("data/",CASE,"/mc_res/Biomes_seafloor_meanResults.Rdata"))

# 2. extract values based on the slowest vs fastest sinking speeds
seafl_slow_cn <- cn_cf_seafl_all[,1,] + cn_eg_seafl_all[,1,]
seafl_slow_ct <- ct_cf_seafl_all[,1,] + ct_eg_seafl_all[,1,]
seafl_slow_ch <- ch_cf_seafl_all[,1,] + ch_eg_seafl_all[,1,]
seafl_fast_cn <- cn_cf_seafl_all[,dim(cn_cf_seafl_all)[2],] + cn_eg_seafl_all[,dim(cn_eg_seafl_all)[2],]
seafl_fast_ct <- ct_cf_seafl_all[,dim(ct_cf_seafl_all)[2],] + ct_eg_seafl_all[,dim(ct_eg_seafl_all)[2],]
seafl_fast_ch <- ch_cf_seafl_all[,dim(ch_cf_seafl_all)[2],] + ch_eg_seafl_all[,dim(ch_eg_seafl_all)[2],]

exclude_shallow=FALSE
cnid_slow=conf_int_by_biome(seafl_slow_cn, seafl_all$Cnidarians$lat, seafl_all$Cnidarians$lon, exclude_shallow = exclude_shallow)
cnid_fast=conf_int_by_biome(seafl_fast_cn, seafl_all$Cnidarians$lat, seafl_all$Cnidarians$lon, exclude_shallow = exclude_shallow)
cten_slow=conf_int_by_biome(seafl_slow_ct, seafl_all$Ctenophores$lat, seafl_all$Ctenophores$lon, exclude_shallow = exclude_shallow)
cten_fast=conf_int_by_biome(seafl_fast_ct, seafl_all$Ctenophores$lat, seafl_all$Ctenophores$lon, exclude_shallow = exclude_shallow)
chor_slow=conf_int_by_biome(seafl_slow_ch, seafl_all$Chordates$lat, seafl_all$Chordates$lon, exclude_shallow = exclude_shallow)
chor_fast=conf_int_by_biome(seafl_fast_ch, seafl_all$Chordates$lat, seafl_all$Chordates$lon, exclude_shallow = exclude_shallow)

slow = cnid_slow[,2:3]+cten_slow[,2:3]+chor_slow[,2:3]
fast = cnid_fast[,2:3]+cten_fast[,2:3]+chor_fast[,2:3]

bci_seafl = data.frame(biome=cnid_slow$biome, slow, fast)
names(bci_seafl) = c("biome", bci_names)
bci_seafl <- join(bci_seafl, av_biomesum, by="biome")
bci_seafl <- cbind(bci_seafl, bci_seafl[,bci_names] * bci_seafl$vol)
names(bci_seafl)[7:10]<-str_c(bci_names,"_total")


# exclude depths < 50 m for evaluating contribution of carcasses to seafloor flux
cnid_cf_slow=conf_int_by_biome(cn_cf_seafl_all[,1,], 
                               seafl_all$Cnidarians$lat, seafl_all$Cnidarians$lon, exclude_shallow = TRUE)
cnid_cf_fast=conf_int_by_biome(cn_cf_seafl_all[,dim(cn_cf_seafl_all)[2],], 
                               seafl_all$Cnidarians$lat, seafl_all$Cnidarians$lon, exclude_shallow = TRUE)
cten_cf_slow=conf_int_by_biome(ct_cf_seafl_all[,1,], 
                               seafl_all$Ctenophores$lat, seafl_all$Ctenophores$lon, exclude_shallow = TRUE)
cten_cf_fast=conf_int_by_biome(ct_cf_seafl_all[,dim(cn_cf_seafl_all)[2],], 
                               seafl_all$Ctenophores$lat, seafl_all$Ctenophores$lon, exclude_shallow = TRUE)
chor_cf_slow=conf_int_by_biome(ch_cf_seafl_all[,1,], 
                               seafl_all$Chordates$lat, seafl_all$Chordates$lon, exclude_shallow = TRUE)
chor_cf_fast=conf_int_by_biome(ch_cf_seafl_all[,dim(cn_cf_seafl_all)[2],], 
                               seafl_all$Chordates$lat, seafl_all$Chordates$lon, exclude_shallow = TRUE)

slow_cf = cnid_cf_slow[,2:3]+cten_cf_slow[,2:3]+chor_cf_slow[,2:3]
fast_cf = cnid_cf_fast[,2:3]+cten_cf_fast[,2:3]+chor_cf_fast[,2:3]

bci_seafl_cf = data.frame(biome=cnid_cf_slow$biome, slow_cf, fast_cf)
names(bci_seafl_cf) = c("biome", bci_names)
bci_seafl_cf <- join(bci_seafl_cf, av_biomesum, by="biome")
bci_seafl_cf <- cbind(bci_seafl_cf, bci_seafl_cf[,bci_names] * bci_seafl_cf$vol)
names(bci_seafl_cf)[7:10]<-str_c(bci_names,"_total")

colSums(bci_seafl_cf[2:10])
write.csv(bci_seafl_cf, str_c("data/",CASE,"/mc_res/seafloor_CIs_carcasses_exclude_shallow.csv"), row.names=FALSE)

# save(bci_seafl, file=str_c("data/",CASE,"/mc_res/Biomes_seafloor_slow_v_fast_quantile.Rdata"))
# load(str_c("data/",CASE,"/mc_res/Biomes_seafloor_slow_v_fast_quantile.Rdata"))

# 3. Construct a single table with all the sinking flux means and ranges
cnid_range = melt(rbind(data.frame(cnid_slow, speed="slow"), data.frame(cnid_fast, speed="fast")), id.vars = c("biome","speed"))
cten_range = melt(rbind(data.frame(cten_slow, speed="slow"), data.frame(cten_fast, speed="fast")), id.vars = c("biome","speed"))
chor_range = melt(rbind(data.frame(chor_slow, speed="slow"), data.frame(chor_fast, speed="fast")), id.vars = c("biome","speed"))

cnid_all = rbind(cnid_range, data.frame(melt(cnid[,c("biome","mean")], id.vars = "biome"), speed="all"))
cten_all = rbind(cten_range, data.frame(melt(cten[,c("biome","mean")], id.vars = "biome"), speed="all"))
chor_all = rbind(chor_range, data.frame(melt(chor[,c("biome","mean")], id.vars = "biome"), speed="all"))

seafl_all = rbind(data.frame(cnid_all, taxon="cnid"), data.frame(cten_all, taxon="cten"), data.frame(chor_all, taxon="chor"))
names(seafl_all) = c("biome","speed","value","flux","taxon")
seafl_all$value = factor(seafl_all$value, labels=c("low","high","mean"))

# 4. Extract global sums for slow vs. fast sinkers confidence intervals
res <- rbind(res, data.frame(depth="Seafloor", speed="slow", t(unname(colSums(bci_seafl[c("slow_lo_total","slow_hi_total")])))))
res <- rbind(res, data.frame(depth="Seafloor", speed="fast", t(unname(colSums(bci_seafl[c("fast_lo_total","fast_hi_total")])))))

names(res) <- c("depth","speed","lower","upper")

write.csv(res, str_c("data/",CASE,"/mc_res/Biomes_Depth_SummedFluxes_CIs.csv"), row.names=FALSE)

## save files
combined <- rbind(data.frame(depth="Export",exprt_all),
                 data.frame(depth="Sequestration",seq_d_all),
                 data.frame(depth="Seafloor",seafl_all))
combined <- join(combined, av_biomesum, by="biome")
combined$flux_tot <- combined$flux * combined$vol
combined <- combined[,c("depth","biome","taxon","speed","value","flux","flux_tot")]

write.csv(combined, str_c("data/",CASE,"/mc_res/Biomes_Depth_Combined_Fluxes_CIs.csv"), row.names=FALSE)

