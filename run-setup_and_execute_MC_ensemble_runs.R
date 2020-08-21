#!/usr/bin/env Rscript

#-------------------------------------------------------
#
#     Set up for analysis of export flux error
#       - From the Monte Carlo accepted parameters
#       - Bootstrap 100 runs 
#
#-------------------------------------------------------

library("plyr")
library("stringr")

`%ni%` <- Negate(`%in%`)

# # import arguments from Unix
args <- commandArgs(trailingOnly=TRUE)
C_ID <- as.numeric(args[1])
BATCH <- as.numeric(args[2])


# create directory
#C_ID = 1
CASES = c("0-baseline","1-high_biomass","2-low_biomass")
CASE <<- CASES[C_ID]

dir.create(str_c("data/",CASE,"/MC_ExportError"), showWarnings = FALSE)
n = 100 #number of ensembles to run

# each BATCH will take SET_NUM runs
# if you set SET_NUM to 100, that would imply 1 BATCH
SET_NUM <- 100
print(str_c("Running case: ", CASE))
print(str_c("Running batch: ", BATCH, " of ", n/SET_NUM))

NEWSETUP <- TRUE
ENSEMBLE_RUN <- TRUE


## { Import Upper Ocean results, set up 100 sims ---------------------------------
# set up simulations folders (only need to be done once)?

if(NEWSETUP){
  print("Setting up new ensembles")
  # What needs to happen here is to pull out results from all the accepted simulations
  # Bootstrap a random selection of 100 runs (greatly overrepresents accepted sims for chordates) - not sure how else to do this though
  # Run the export flux on each of the 100 runs
  
  # Import all MC accepted simulations
  # Load
  load(str_c("data/",CASE,"/baseline/results_cnid_allres_acceptedSims.Rdata"))
  load(str_c("data/",CASE,"/baseline/results_cten_allres_acceptedSims.Rdata"))
  load(str_c("data/",CASE,"/baseline/results_chor_allres_acceptedSims.Rdata"))
  
  # quickly examine the structure of each file
  str(dcnid_allres)
  # num [1:2557, 1:9, 1:3141] -62.5 -62.5 -61.5 -51.5 -50.5 -45.5 -43.5 -42.5 -42.5 -42.5 ...
  # - attr(*, "dimnames")=List of 3
  # ..$ : NULL
  # ..$ : chr [1:9] "lat" "lon" "ingestion" "respiration" ...
  # ..$ : NULL
  dimnames(dcnid_allres)[2]
  # [1] "lat"          "lon"          "ingestion"    "respiration"  "DOC"         
  # [6] "reproduction" "predation"    "flux"         "egestion"        
  str(dcten_allres)
  # num [1:774, 1:9, 1:3553] -66.5 -65.5 -65.5 -64.5 -64.5 -64.5 -63.5 -63.5 -62.5 -62.5 ...
  # - attr(*, "dimnames")=List of 3
  # ..$ : NULL
  # ..$ : chr [1:9] "lat" "lon" "ingestion" "respiration" ...
  # ..$ : NULL
  str(dchor_allres)
  # num [1:2881, 1:9, 1:40] -65.5 -63.5 -62.5 -60.5 -60.5 -59.5 -59.5 -59.5 -58.5 -58.5 ...
  # - attr(*, "dimnames")=List of 3
  # ..$ : NULL
  # ..$ : chr [1:9] "lat" "lon" "ingestion" "respiration" ...
  # ..$ : NULL
  ## --> Looks like 3-D array, 1st dimension: lat/lon grid, 2nd dimension: all the fluxes, and 3rd dimension: accepted sims

  # draw n random (with replacement) accepted sims, set a seed
  set.seed(0)
  
  draw <- data.frame(cnid=sample(seq_len(dim(dcnid_allres)[3]), n, replace = FALSE), 
                     cten=sample(seq_len(dim(dcten_allres)[3]), n, replace = FALSE), 
                     chor=sample(seq_len(dim(dchor_allres)[3]), n, replace = FALSE)) 
  
  
  # original dataframes are as follows:
  # cf <- read.csv(str_c("data/",CASE,"/baseline/results_allMedian_ML_flux.csv"))
  # eg <- read.csv(str_c("data/",CASE,"/baseline/results_egestion_ML_flux.csv"))
  
  # set up dataframes for the values from the bootstrapped accepted sims
  # these will directly go into the code for export flux
  # two dataframes, one for the export flux (cf) and one for the egestion (eg)
  # this loop takes approx. 5 minutes to run
  
  for(k in 1:n){
    print(k)
    
    # cnidarians
    CN_SIM_N <- draw$cnid[k]
    cn <- data.frame(dcnid_allres[,,CN_SIM_N], taxon="cnid")
    cn_cf <- as.data.frame(dcnid_allres[,c("lat", "lon", "flux"), CN_SIM_N])
    cn_eg <- as.data.frame(dcnid_allres[,c("lat", "lon", "egestion"), CN_SIM_N])
    names(cn_cf)[3] <- names(cn_eg)[3] <- "Cnidaria"
    
    # ctenophores
    CT_SIM_N <- draw$cten[k]
    ct <- data.frame(dcten_allres[,,CT_SIM_N], taxon="cten")
    ct_cf <- as.data.frame(dcten_allres[,c("lat", "lon", "flux"), CT_SIM_N])
    ct_eg <- as.data.frame(dcten_allres[,c("lat", "lon", "egestion"), CT_SIM_N])
    names(ct_cf)[3] <- names(ct_eg)[3] <- "Ctenophora"
    
    # chordates
    CH_SIM_N <- draw$chor[k]
    ch <- data.frame(dchor_allres[,,CH_SIM_N], taxon="chor")
    ch_cf <- as.data.frame(dchor_allres[,c("lat", "lon", "flux"), CH_SIM_N])
    ch_eg <- as.data.frame(dchor_allres[,c("lat", "lon", "egestion"), CH_SIM_N])
    names(ch_cf)[3] <- names(ch_eg)[3] <- "Chordata"
    
    # join together
    cf <- join_all(list(cn_cf, ct_cf, ch_cf), by=c("lat", "lon"), type = "full") # join_all is a plyr function
    eg <- join_all(list(cn_eg, ct_eg, ch_eg), by=c("lat", "lon"), type = "full")
    uo_res <- rbind(cn, ct, ch)
  
    
    # create directories and save
    dir.create(str_c("data/",CASE,"/MC_ExportError/",k), showWarnings = FALSE)
    save(cf, file=str_c("data/",CASE,"/MC_ExportError/",k,"/results_MCsim_biomass_flux.Rdata")) # upper ocean flux
    save(eg, file=str_c("data/",CASE,"/MC_ExportError/",k,"/results_MCsim_egestion_flux.Rdata"))
    save(uo_res, file=str_c("data/",CASE,"/MC_ExportError/",k,"/results_MCsim_uo_res.Rdata"))
  }
  
}

# }

## { Run simulations for export flux ---------------------------------------------------

if (ENSEMBLE_RUN){
  for (j in 1:SET_NUM){
    
    # set run number
    k <- (BATCH-1) * SET_NUM + j


    print(paste("starting ensemble run", k, "at", Sys.time()))

    ptm <- proc.time() # starting time

    # load
    load(file=str_c("data/",CASE,"/MC_ExportError/",k,"/results_MCsim_biomass_flux.Rdata"))
    load(file=str_c("data/",CASE,"/MC_ExportError/",k,"/results_MCsim_egestion_flux.Rdata"))

    # set as global variables
    cf <<- cf
    eg <<- eg

    # set up DIRNAME
    DIRNAME <<- str_c("MC_ExportError/",k)

    # run
    source("run-export_flux.R")

    print(proc.time() - ptm) # end time
  }
}



# }


