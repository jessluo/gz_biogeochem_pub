#
#           Upper Ocean C-flux model-specific functions
#           1) carbon flux model
#           2) egestion flux 
#           3) a set of functions to fix the negative fluxes
#
#------------------------------------------------------------------------


# reading packages
library(ggplot2)
library(plyr)
library(reshape2)
library(oce)
library(gridExtra)
library(stringr)

Ea = -0.65
BOLTZMANN_K = 8.6173303E-5 #eV/K

##{ Carbon flux specific functions -------------------------------------
carbon_flux <- function(std_Biomass, numberInd, temp, prey_md, prey_lg, prey_prod_md, prey_prod_lg,
                        esd_slope, esd_intercept, cr_slope, cr_intercept, resp_slope, resp_intercept, 
                        coeff_temp=0.96, coeff_DOC=0.25, coeff_repro=0.2, predation_ee = 0.45, 
                        AE=0.8, gge=0.35, type=c("cnid","cten","chor")){
  
  OVERCONSUMPTION_THRESHOLD = 0.7
  NEGATIVEPRODUCTION_THRESHOLD = 0.4

  c <- std_Biomass
  
  NAdf = data.frame(array(NA,dim=c(length(c),7)))
  names(NAdf) = c("ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")
  
  # calculate the equiv. spherical diameter
  logESD <- esd_slope * log10(c) + esd_intercept
  
  # calculate the prey field
  # units of mg C m^-3
  PREY <- prey_md + prey_lg
  
  TEMP <- coeff_temp * temp
  
  # # GROWTH RATES
  # umax_corr <- -0.34 * logESD + 10.38
  # # convert umax_corr to umax
  # umax <- 10^(umax_corr) * exp((Ea/BOLTZMANN_K) / (273 + TEMP)) 
  # # half saturation const, vary by order of 5
  # k <- seq(from=1, to=50, by=2)
  # i <- 10
  # 
  # growth <- AE * umax * (PREY/(k[i]+PREY)) * c
  
  # CLEARANCE RATES (rates are in L ind-1 d-1)
  
  if(type %in% c("cnid","cten")){cr_c=c/1000} else {cr_c=c} # for cnid/cten biomass in equation is g C
    
  clearance_rate <- 10^(cr_intercept + cr_slope * log10(cr_c))
  # temperature correction
  clearance_rate <- clearance_rate * exp((Ea/BOLTZMANN_K) / (273 + TEMP))
  # convert to m^3 ind-1 d-1
  clearance_rate <- clearance_rate / 1000
  
  # INGESTION RATE
  consumption <- numberInd * (PREY) * clearance_rate
  
  # Cannot consume more than prey production
  TOTALPREYPROD <- prey_prod_md + prey_prod_lg
  dd <- which(consumption > TOTALPREYPROD)
  # To select for the appropriate clearance rates & prey capture parameters, 
  #    only allow for consumption > prey production in a certain % of cells
  #    we are calling this the overconsumption threshold amount
  if(length(dd) > 0){
    if(length(dd)/length(c) > OVERCONSUMPTION_THRESHOLD){
      return(NAdf)
    } else {
      consumption[dd] <- TOTALPREYPROD[dd]
    }
  }
  
  # ASSIMILATED CONSUMPTION = INGESTION
  ingestion <- AE * consumption
  egestion <- (1-AE) * consumption
  
  # RESPIRATION RATE
  # NB: cnidarians and ctenophores are going to use log(ESD)
  # NB: chordates: going to use log(Biomass)
  if(type=="chor"){
    respiration <- resp_slope * log10(c) + resp_intercept
  } else {
    respiration <- resp_slope * logESD + resp_intercept
  }
  
  respiration <- 10^respiration 
  # in units of ml O2 ind-1 hr-1 ==> must convert to Carbon
  respiration <- respiration  * exp((Ea/BOLTZMANN_K) / (273 + TEMP))
  # stoichiometric conversion: 1 ml O2 == 0.4287 mg C with a Respiratory Quotient of 0.8
  # JYL rechecked 06-04-2019:
  # 1 ml O2 * 0.04462 mmol O2/ml O2 [gas volume of O2] * 0.8 CO2/O2 * 12.011 mg C/mmolC = 0.4287
  respiration  <- respiration  * 0.4287
  # in mg C hr-1 so convert to mg C day-1
  respiration <- respiration * 24
  respiration <- numberInd * respiration
  
  # production
  production <- ingestion - respiration
  
  idx_negs = which(production < 0)
  world_fraction_neg_production = length(idx_negs)/length(production)
  
  if (world_fraction_neg_production > NEGATIVEPRODUCTION_THRESHOLD){
    return(NAdf)
  }
  
  # if negatives are not more than the negative production threshold #
  respiration[idx_negs] = ingestion[idx_negs]
  production[idx_negs] <- 0
  
  # RI_ratio = mean(respiration/ingestion, na.rm=TRUE)
  # if(RI_ratio < 0.2){
  #   return(NAdf)
  # }
  
  # Route the rest of the production to exudation, reproduction, predation, and non-predation mortality
  route_sum = predation_ee + coeff_DOC + coeff_repro
  if(route_sum > 1){
    predation_ee = predation_ee / route_sum
    coeff_DOC = coeff_DOC / route_sum
    coeff_repro = coeff_repro / route_sum
  }
  
  predation <- predation_ee * production
  exudation <- coeff_DOC * production
  reproductive_loss <- coeff_repro * production
  
  mortality <- (1-route_sum) * production
  
  # put it all together into a data frame
  allres <- data.frame(ingestion = ingestion,
                       respiration = respiration,
                       DOC = exudation,
                       reproduction = reproductive_loss,
                       predation = predation,
                       flux = mortality,
                       egestion = egestion)
  
  # return flux in g C
  allres <- allres / 1000
  
  # return flux in g C yr-1
  allres <- allres * 365
  
  return(allres)
}



check_high_mesozoop_predation <- function(dcnid, dcten, prey_prod){
  # convert to data frame
  dcnid_df = as.data.frame(dcnid)
  dcten_df = as.data.frame(dcten)

  names(dcten_df) <- names(dcnid_df) <- c("lat", "lon", "ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")

  if(ncol(prey_prod) > 3){
    prey_prod$prey_prod <- apply(prey_prod[,3:ncol(prey_prod)], 1, sum)
  } else {prey_prod$prey_prod <- prey_prod[,3]}

  dcnid_df <- plyr::join(dcnid_df, prey_prod[,c("lat","lon","prey_prod")], by=c("lat","lon"))
  # dcten_df <- plyr::join(dcten_df, prey_prod[,c("lat","lon","prey_prod")], by=c("lat","lon"))
  # 
  # dcnid_df$totconsumption = dcnid_df$ingestion + dcnid_df$egestion
  # dcten_df$totconsumption = dcten_df$ingestion + dcten_df$egestion
  
  idx_cnct <- join(dcnid_df[,c("lat", "lon")], dcten_df[,c("lat", "lon")], type="inner", by=c("lat","lon"))
  dcn <- join(idx_cnct, dcnid_df, by=c("lat","lon"))
  dct <- join(idx_cnct, dcten_df, by=c("lat","lon"))

  TOTAL_PREY_PROD = dcn$prey_prod

  total_consumption = dcn$ingestion + dcn$egestion + dct$ingestion + dct$egestion

  dd = which(total_consumption > TOTAL_PREY_PROD)

  # if(length(dd) == 0){
  #   return(list(dcnid=dcnid, dcten=dcten))
  # } else {
  #   return(list(dcnid=NA, dcten=NA))
  # }

  return(length(dd)/nrow(idx_cnct))
}

# partition_cnid_cten_fluxes <- function(dcnid, dcten, prey_prod_md_cnidx, prey_prod_lg_cnidx){
#   # convert to data frame
#   dcnid_df = as.data.frame(dcnid)
#   dcten_df = as.data.frame(dcten)
#   
#   names(dcten_df) <- names(dcnid_df) <- c("lat", "lon", "ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")
#   
#   dcnid_df$prey_prod_md = prey_prod_md_cnidx
#   dcnid_df$prey_prod_lg = prey_prod_lg_cnidx
#   
#   idx_cnct <- join(dcnid_df[,c("lat", "lon")], dcten_df[,c("lat", "lon")], type="inner", by=c("lat","lon"))
#   dcn <- join(idx_cnct, dcnid_df, by=c("lat","lon"))
#   dct <- join(idx_cnct, dcten_df, by=c("lat","lon"))
#   
#   TOTAL_PREY_PROD = dcn$prey_prod_md + dcn$prey_prod_lg
#   
#   total_consumption = dcn$ingestion + dcn$egestion + dct$ingestion + dct$egestion
#   
#   dd = which(total_consumption > TOTAL_PREY_PROD)
# }

# log.and.plot <- function(data, vars, plotname=NULL, plot.limits=NULL){
#   # function to calculate the log of the biomass and then to plot it and save it
#   latlon <- data[,LATLON]
#   # plotting loop
#   for (i in 1:length(vars)){
#     var <- vars[i]
#     temp <- data[,which(names(data)==var)]
#     temp <- as.data.frame(log1p(temp))
#     names(temp) <- var
#     temp <- cbind(latlon, temp)
#     p <- plotting(var=vars[i], data=temp, plot.title = paste0("Log Biomass in ", var), plot.limits = plot.limits)
#     #  p <- p + scale_colour_gradientn(na.value="white")
#     ggsave(paste0("plots/", plotname, i, "_", vars[i], ".pdf"), p)
#   }
# }

 

# }
