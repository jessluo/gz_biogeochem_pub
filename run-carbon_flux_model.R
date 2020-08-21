#!/usr/bin/env Rscript
#
#           Calculating flux of C through GZ
#           Model formulation and construction
#
#           (c) Jessica Luo
#
#
#          To complete full runs with Monte Carlo simulations:
#          run from command line:
#              ./run-carbon_flux_model [C-ID NUMBER]
#------------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
C_ID <- as.numeric(args[1])
# C_ID = 1
CASES = c("0-baseline","1-high_biomass","2-low_biomass")
CASE = CASES[C_ID]

print(paste0("running case: ", CASE))

# reading packages
library(ggplot2)
library(plyr)
library(reshape2)
library(oce)
library(gridExtra)
library(stringr)

source("functions.R")
source("functions_c-flux-model.R")


DIRNAME <- "baseline"
rd <- read.csv(str_c("data/", CASE, "/jelly_biomass_1_deg_grid.csv"), stringsAsFactors=FALSE)

RUNMCSIMS <- TRUE

USESALPPHYTO <- TRUE

LATLON <- c("lat", "lon")
taxa <-  c("Cnidaria", "Ctenophora", "Chordata")


##{ Read in & set up data -------------------------------------------------------
# processed in prepare-fill_NAs_gz_model_inputs.R
load(file = "data/gz_model_inputs/rsst.Rdata")
load(file = "data/gz_model_inputs/rzoomass.Rdata")
load(file = "data/gz_model_inputs/rzooprod.Rdata")
load(file = "data/gz_model_inputs/phyto.Rdata")
load(file = "data/gz_model_inputs/pp.Rdata")
load(file = "data/gz_model_inputs/phyto_salps.Rdata")
load(file = "data/gz_model_inputs/pp_salps.Rdata")

# are we going to use the salp-specific phytoplankton fields?
if(USESALPPHYTO){
  phyto <- phyto_salps
  names(phyto)[3] <- 'phytoC'
  pp <- pp_salps
  names(pp)[3] <- 'pp'
}

# CNIDARIANS
# rd$stdSize is now a "Standardized biomass" term. Biomass per unit individual
whichcnid <- which(rd$rank_phylum=="Cnidaria")
cnid <- rd[whichcnid, c(LATLON, "numeric_density", "Biomass.mgCm3", "n", "stdBiomass")]
row.names(cnid) <- 1:nrow(cnid)
cnid <- plyr::rename(cnid, replace=c("n" = "nrows"))

# CTENOPHORES 
whichcten <- which(rd$rank_phylum=="Ctenophora")
cten <- rd[whichcten, c(LATLON, "numeric_density", "Biomass.mgCm3", "n", "stdBiomass")]
row.names(cten) <- 1:nrow(cten)
cten <- plyr::rename(cten, replace=c("n" = "nrows"))

# CHORDATES
whichchor <- which(rd$rank_phylum=="Chordata")
chor <- rd[whichchor, c(LATLON, "numeric_density", "Biomass.mgCm3", "n", "stdBiomass")]
row.names(chor) <- 1:nrow(chor)
chor <- plyr::rename(chor, replace=c("n" = "nrows"))

# }

## Simulation models for baseline ----------------------------------------------------------

if(RUNMCSIMS == TRUE){
  # read in simulations
  sim <- read.csv("data/MonteCarloSimulations.csv")

  saveidx1 = NULL
  ## Cnidarians
  dcnid_allres <- array(NA, dim = c(length(whichcnid), 9, nrow(sim)))
  for (i in 1:nrow(sim)){
    temp1 <- carbon_flux(std_Biomass = cnid$stdBiomass, numberInd = cnid$numeric_density, temp = rsst[whichcnid,"ann"],
                        prey_md = rzoomass[whichcnid, "md_zoomass"], prey_lg = rzoomass[whichcnid, "lg_zoomass"],
                        prey_prod_md = rzooprod[whichcnid, "md_zooprod"], prey_prod_lg = rzooprod[whichcnid, "lg_zooprod"],
                        esd_slope = sim$ESD_SLOPE[i], esd_intercept = sim$ESD_INTERCEPT[i],
                        cr_slope = sim$CR_SLOPE[i], cr_intercept = sim$CR_INTERCEPT[i],
                        resp_slope = sim$CNID_RESP_SLOPE[i], resp_intercept = sim$CNID_RESP_INTERCEPT[i],
                        coeff_temp = sim$coeff_cn_temp[i], 
                        coeff_DOC = sim$coeff_cn_DOC[i], coeff_repro = sim$coeff_cn_repro[i],
                        predation_ee = sim$coeff_cn_predEE[i],
                        AE=sim$ae_cn_ct[i], type="cnid")
     if(! all(is.na(temp1))){
       saveidx1 = c(saveidx1,i) }
  
    dcnid_allres[,,i] <- as.matrix(cbind(rd[whichcnid, LATLON], temp1))
  }
  
  print(length(saveidx1))
  colMeans(sim[saveidx1,])/colMeans(sim)

  
  
  saveidx2 = NULL
  ## Ctenophores (MC sim)
  dcten_allres <- array(NA, dim = c(length(whichcten), 9, nrow(sim)))
  for (i in 1:nrow(sim)){
    temp2 <- carbon_flux(std_Biomass = cten$stdBiomass, numberInd = cten$numeric_density, temp = rsst[whichcten,"ann"],
                        prey_md = rzoomass[whichcten, "md_zoomass"], prey_lg = rzoomass[whichcten, "lg_zoomass"],
                        prey_prod_md = rzooprod[whichcten, "md_zooprod"], prey_prod_lg = rzooprod[whichcten, "lg_zooprod"],
                        esd_slope = sim$ESD_SLOPE[i], esd_intercept = sim$ESD_INTERCEPT[i],
                        cr_slope = sim$CR_SLOPE[i], cr_intercept = sim$CR_INTERCEPT[i],
                        resp_slope = sim$CTEN_RESP_SLOPE[i], resp_intercept = sim$CTEN_RESP_INTERCEPT[i],
                        coeff_temp = sim$coeff_ct_temp[i], 
                        coeff_DOC = sim$coeff_ct_DOC[i], coeff_repro = sim$coeff_ct_repro[i],
                        predation_ee = sim$coeff_ct_predEE[i],
                        AE=sim$ae_cn_ct[i], type="cten")

    if(! all(is.na(temp2))){
      saveidx2 = c(saveidx2,i)  }

    dcten_allres[,,i] <- as.matrix(cbind(rd[whichcten, LATLON], temp2))
  }

  print(length(saveidx2))
  colMeans(sim[saveidx2,])/colMeans(sim)
  
  
  saveidx3 = NULL
  ## Chordates
  dchor_allres <- array(NA, dim = c(length(whichchor), 9, nrow(sim)))
  for (i in 1:nrow(sim)){
    temp3 <- carbon_flux(std_Biomass = chor$stdBiomass, numberInd = chor$numeric_density, temp = rsst[whichchor,"ann"],
                        prey_md = phyto[whichchor, "phytoC"], prey_lg = 0.33 * rzoomass[whichchor, "sm_zoomass"],
                        prey_prod_md = pp[whichchor, "pp"], prey_prod_lg = 0.33 * rzooprod[whichchor, "sm_zooprod"],
                        esd_slope = sim$ESD_SLOPE[i], esd_intercept = sim$ESD_INTERCEPT[i],
                        cr_slope = sim$CH_CR_SLOPE[i], cr_intercept = sim$CH_CR_INTERCEPT[i],
                        resp_slope = sim$CHOR_RESP_SLOPE[i], resp_intercept = sim$CHOR_RESP_INTERCEPT[i],
                        coeff_temp = sim$coeff_ch_temp[i], 
                        coeff_DOC = sim$coeff_ch_DOC[i], coeff_repro = sim$coeff_ch_repro[i],
                        predation_ee = sim$coeff_ch_predEE[i], 
                        AE = sim$ae_chor[i], type="chor")

    if(! all(is.na(temp3))){
      saveidx3 = c(saveidx3,i)  }
    
    dchor_allres[,,i] <- as.matrix(cbind(rd[whichchor, LATLON], temp3))
  }
  
  print(length(saveidx3))
  colMeans(sim[saveidx3,])/colMeans(sim)
  
  
  # check to make sure cnidarian + ctenophore consumption does not exceed mesozooplankton production
  check_idx = intersect(saveidx1,saveidx2)
  prey_prod = rzooprod[!duplicated(rzooprod),c("lat","lon","md_zooprod","lg_zooprod")]
  for (i in check_idx){
    len_dd = check_high_mesozoop_predation(dcnid = dcnid_allres[,,i], dcten = dcten_allres[,,i], 
                                               prey_prod = prey_prod)
    if(len_dd != 0){
      print(str_c('exceeded production: sim #', i))
      # partitioned = partition_cnid_cten_fluxes(dcnid = dcnid_allres[,,i], dcten = dcten_allres[,,i], 
      #                                          prey_prod_md_cnidx = rzooprod[whichcnid, "md_zooprod"],
      #                                          prey_prod_lg_cnidx = rzooprod[whichcnid, "lg_zooprod"])
      # dcnid_allres[,3:9,i] = partitioned['cnid']
      # dcten_allres[,3:9,i] = partitioned['cten']
      dcnid_allres[,3:9,i] = partitioned['cnid']
      dcten_allres[,3:9,i] = partitioned['cten']
      saveidx1 = saveidx1[which(saveidx1 != i)]
      saveidx2 = saveidx2[which(saveidx2 != i)]
    }
  }
  
  ## More culling of the simulation models 
  
  # CRITERIA 3: Simulations must result in biologically realistic Respiration to Ingestion ratios
  # Respiration to Ingestion ratios need to fall within normal ranges - these ranges are broad
  # but prevent respiration from being < 10% of all assimilated energy
  
  # pull out the respiration / ingestion ratios, and flux calculations
  # column 4 is respiration, column 3 is ingestion
  dcnidRI <- dcnid_allres[,4,saveidx1] / dcnid_allres[,3,saveidx1]
  dctenRI <- dcten_allres[,4,saveidx2] / dcten_allres[,3,saveidx2]
  dchorRI <- dchor_allres[,4,saveidx3] / dchor_allres[,3,saveidx3]
  
  CN_RI <- c(0.2, 0.8)
  CT_RI <- c(0.2, 0.8)
  CH_RI <- c(0.2, 0.8)
  
  # calculate the mean of the simulations' respiration to ingestion
  dcnidRImean <- colMeans(dcnidRI, na.rm=TRUE)
  dctenRImean <- colMeans(dctenRI, na.rm=TRUE)
  dchorRImean <- colMeans(dchorRI, na.rm=TRUE)
  
  # Accept or reject simulations based on those three criteria
  # find the intersection between rule 1 and 3
  idx_cn <- which(dcnidRImean > CN_RI[1] & dcnidRImean < CN_RI[2]); print(length(idx_cn)) # baseline: 247; high: 206; low: 275
  idx_ct <- which(dctenRImean > CT_RI[1] & dctenRImean < CT_RI[2]); print(length(idx_ct)) # baseline: 371; high: 293; low: 468
  idx_ch <- which(dchorRImean > CH_RI[1] & dchorRImean < CH_RI[2]); print(length(idx_ch)) # baseline: 6430; high: 3965; low: 8712
  
  idx_cn = saveidx1[idx_cn]
  idx_ct = saveidx2[idx_ct]
  idx_ch = saveidx3[idx_ch]
  
  
  # Keep just the accepted simluations
  dcnid_allres <- dcnid_allres[,,idx_cn]
  dcten_allres <- dcten_allres[,,idx_ct]
  dchor_allres <- dchor_allres[,,idx_ch]
  
  
  dimnames(dcnid_allres)[[2]] <- c("lat", "lon", "ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")
  dimnames(dcten_allres)[[2]] <- c("lat", "lon", "ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")
  dimnames(dchor_allres)[[2]] <- c("lat", "lon", "ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")


  
  # Posteriors = "Accepted Coefficient Values"
  simMelt <- melt(sim, id.vars="V1")
  simMelt <- simMelt[,-1]
  simMelt$id <- "sim"
  
  # grab the coefficients in the accepted simulations
  CNIDSIMS <- c("ESD_SLOPE", "ESD_INTERCEPT", "CR_SLOPE", "CR_INTERCEPT", "coeff_cn_temp", "CNID_RESP_SLOPE", 
                "CNID_RESP_INTERCEPT", "coeff_cn_DOC", "coeff_cn_repro", "coeff_cn_predEE", "ae_cn_ct") 
  cnid_posterior <- sim[idx_cn,CNIDSIMS]
  cnid_posteriorM <- melt(cnid_posterior)
  cnid_posteriorM$id <- "cnid"
  
  CTENSIMS <- c("ESD_SLOPE", "ESD_INTERCEPT", "CR_SLOPE", "CR_INTERCEPT", "coeff_ct_temp", "CTEN_RESP_SLOPE", 
                "CTEN_RESP_INTERCEPT", "coeff_ct_DOC", "coeff_ct_repro", "coeff_ct_predEE", "ae_cn_ct") 
  cten_posterior <- sim[idx_ct,CTENSIMS]
  cten_posteriorM <- melt(cten_posterior)
  cten_posteriorM$id <- "cten"
  
  CHORSIMS <- c("CH_CR_SLOPE", "CH_CR_INTERCEPT", "coeff_ch_temp", "CHOR_RESP_SLOPE", "CHOR_RESP_INTERCEPT", 
                "coeff_ch_DOC", "coeff_ch_repro", "coeff_ch_predEE", "ae_chor") 
  chor_posterior <- sim[idx_ch,CHORSIMS]
  chor_posteriorM <- melt(chor_posterior)
  chor_posteriorM$id <- "chor"
  
  posterior <- rbind(simMelt, cnid_posteriorM, cten_posteriorM, chor_posteriorM)
  
  posterior$variable <- factor(posterior$variable, levels=c("CR_SLOPE", "CR_INTERCEPT", "CH_CR_SLOPE", "CH_CR_INTERCEPT", 
                                        "ESD_SLOPE", "ESD_INTERCEPT", "CNID_RESP_SLOPE", "CNID_RESP_INTERCEPT", 
                                        "CTEN_RESP_SLOPE", "CTEN_RESP_INTERCEPT", "CHOR_RESP_SLOPE", "CHOR_RESP_INTERCEPT", 
                                        "ae_chor", "ae_cn_ct", "coeff_cn_temp", "coeff_ct_temp", "coeff_ch_temp", 
                                        "coeff_cn_DOC", "coeff_ct_DOC", "coeff_ch_DOC", 
                                        "coeff_cn_repro", "coeff_ct_repro", "coeff_ch_repro", 
                                        "coeff_cn_predEE", "coeff_ct_predEE", "coeff_ch_predEE"))
  posterior$id <- factor(posterior$id, levels=c("sim", "cnid", "cten", "chor"))

  
  col <- c("#000000", "#377eb8", "#4daf4a", "#e41a1c")
  labels <- c("Sim. Inputs", "Cnidarians", "Ctenophores", "Tunicates")
  p <- ggplot(posterior) + geom_density(aes(x=value, y=..density.., fill=id, colour=id), alpha=0.5, stat = "bin", bins=30) + 
    facet_wrap(~variable, scales="free", ncol=4) + scale_fill_manual("", values=col, labels=labels) + 
    scale_color_manual("",values=col, labels=labels) + theme_bw() + labs(x="", y="") + 
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid=element_blank())
  
  suppressMessages(ggsave(str_c("plots/",CASE,"/baseline/MCsim_posteriors.pdf"), plot = p, width=11, height=11, units="in"))
  
  # SAVE RESULTS
  dcnid <- dcnid_allres[,"flux",]
  dcten <- dcten_allres[,"flux",]
  dchor <- dchor_allres[,"flux",]

  # save the id values of the accepted sims
  save(idx_cn, file=str_c("data/",CASE,"/baseline/idx_cn.Rdata"))
  save(idx_ct, file=str_c("data/",CASE,"/baseline/idx_ct.Rdata"))
  save(idx_ch, file=str_c("data/",CASE,"/baseline/idx_ch.Rdata"))

  # all fluxes
  save(dcnid_allres, file=str_c("data/",CASE,"/baseline/results_cnid_allres_acceptedSims.Rdata"))
  save(dcten_allres, file=str_c("data/",CASE,"/baseline/results_cten_allres_acceptedSims.Rdata"))
  save(dchor_allres, file=str_c("data/",CASE,"/baseline/results_chor_allres_acceptedSims.Rdata"))

  # just export
  save(dcnid, file=str_c("data/",CASE,"/baseline/results_cnid_MLflux_acceptedSims.Rdata"))
  save(dcten, file=str_c("data/",CASE,"/baseline/results_cten_MLflux_acceptedSims.Rdata"))
  save(dchor, file=str_c("data/",CASE,"/baseline/results_chor_MLflux_acceptedSims.Rdata"))

  # accepted coeff values ("posteriors")
  save(cnid_posterior, file=str_c("data/",CASE,"/baseline/simulation_posteriors_cnid.Rdata"))
  # load(file=str_c("data/",CASE,"/baseline/simulation_posteriors_cnid.Rdata"))
  save(cten_posterior, file=str_c("data/",CASE,"/baseline/simulation_posteriors_cten.Rdata"))
  # load(file=str_c("data/",CASE,"/baseline/simulation_posteriors_cten.Rdata"))
  save(chor_posterior, file=str_c("data/",CASE,"/baseline/simulation_posteriors_chor.Rdata"))
  # load(file=str_c("data/",CASE,"/baseline/simulation_posteriors_chor.Rdata"))
  save(posterior, file=str_c("data/",CASE,"/baseline/simulation_posteriors.Rdata"))


}

# source('functions.R')
# load('data/biomes_chl_mixedlayer.Rdata')
# av <- read.csv("data/surfaceArea_volume_1-deg.csv")
# av <- join(av,biomes_df, by=c("lat","lon")) ; av$depth.m = av$vol / av$area.m2
# av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(vol=sum(x$vol)))})
# names(av_biomesum) <- c("biome", "vol")
# 
# # #CASE = "0-baseline"
# # load(str_c("data/",CASE,"/baseline/results_cnid_allres_acceptedSims_TEMP.Rdata"))
# # load(str_c("data/",CASE,"/baseline/results_cten_allres_acceptedSims_TEMP.Rdata"))
# # load(str_c("data/",CASE,"/baseline/results_chor_allres_acceptedSims_TEMP.Rdata"))
# 
# set.seed(0)
# SAMPLE=45
# dcnid_allres = dcnid_allres[,,sample(1:dim(dcnid_allres)[3], SAMPLE)]
# dcten_allres = dcten_allres[,,sample(1:dim(dcten_allres)[3], SAMPLE)]
# dchor_allres = dchor_allres[,,sample(1:dim(dchor_allres)[3], SAMPLE)]
# 
# mc_uores = list()
# for (i in 1:SAMPLE){
#   cn <- data.frame(dcnid_allres[,,i], taxon="cnid", stringsAsFactors = FALSE)
#   ct <- data.frame(dcten_allres[,,i], taxon="cten", stringsAsFactors = FALSE)
#   ch <- data.frame(dchor_allres[,,i], taxon="chor", stringsAsFactors = FALSE)  
#   uo_res <- rbind(cn,ct,ch)
#   uo_res <- join(uo_res, biomes_df, by=c("lat","lon"))
#   mc_uores[[i]] = uo_res
# }
# 
# 
# res <- data.frame()
# for (taxon in c("cnid","cten","chor")){
#   tres <- calculate_biome_mean_UO(mc_uores, taxon=taxon, biomes=biomes_df, n_biome=4, n_ens=SAMPLE, quantiles=FALSE, conf_int=TRUE, saveFiles=FALSE)
#   tres$taxon <- taxon
#   res <- rbind(res, tres)
#   print(str_c(CASE, ", ", taxon)); print(ddply(tres, ~var, function(x){return(sum(x$mean_tot))}))
# }
# 
# sum_res <- ddply(res, ~biome+var, function(x){
#   d = colSums(x[,c("mean", "lower", "upper", "vol", "mean_tot", "lower_tot", "upper_tot")])
#   d = as.data.frame(t(d))
#   d$vol = x$vol[1]
#   d$taxon = "sum"
#   return(d)
# })
# print(str_c(CASE, ", sum")); print(ddply(sum_res, ~var, function(x){return(sum(x$mean_tot))}))
# 
# cn=ddply(res[res$taxon=="cnid",], ~var, function(x){return(sum(x$mean_tot))})
# ct=ddply(res[res$taxon=="cten",], ~var, function(x){return(sum(x$mean_tot))})
# ch=ddply(res[res$taxon=="chor",], ~var, function(x){return(sum(x$mean_tot))})
# s=ddply(sum_res, ~var, function(x){return(sum(x$mean_tot))})
# ds=join_all(list(cn,ct,ch,s),by = "var")
# 
# names(ds)=c("var","cnidarians","ctenophores","tunicates","sum")
# ds
# (ds[1,]-ds[2,])
# ds[5,]/(ds[1,]-ds[2,])
