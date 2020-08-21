library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)

set.seed(123)

SIM_NUM <- 30000
sim <- as.data.frame(matrix(seq(1:SIM_NUM), ncol = 1, nrow=SIM_NUM))
# generate 10,000 random values for all the different coefficients

# 8a. Distribution of values for slopes/intercepts of C -> ESD conversions
ESD_SLOPE <- 0.327
ESD_SLOPE_SD <- 0.02147 * sqrt(34) * 10
ESD_INTERCEPT <- -0.0981
ESD_INTERCEPT_SD <- 0.05039 * sqrt(34) * 10

sim$ESD_SLOPE <- rnorm(SIM_NUM, mean = ESD_SLOPE, sd = ESD_SLOPE_SD)
sim$ESD_INTERCEPT <- rnorm(SIM_NUM, mean = ESD_INTERCEPT, sd = ESD_INTERCEPT_SD)

# 8b. Distribution of values for slopes/intercepts of CR calculations
CR_SLOPE <- 0.814
CR_SLOPE_SD <- 0.027 * sqrt(364)
CR_INTERCEPT <- 15.807
CR_INTERCEPT_SD <- 0.065 * sqrt(364)
CH_CR_SLOPE <- 1.20397
CH_CR_SLOPE_SD <- 0.02389 * sqrt(562) 
CH_CR_INTERCEPT <- 12.21788
CH_CR_INTERCEPT_SD <- 0.02206 * sqrt(562) 

## MODIFY
## CR_SLOPE <- CR_SLOPE / 20

sim$CR_SLOPE <- rnorm(SIM_NUM, mean = CR_SLOPE, sd = CR_SLOPE_SD)
sim$CR_INTERCEPT <- rnorm(SIM_NUM, mean = CR_INTERCEPT, sd = CR_INTERCEPT_SD)
sim$CH_CR_SLOPE <- rnorm(SIM_NUM, mean = CH_CR_SLOPE, sd = CH_CR_SLOPE_SD)
sim$CH_CR_INTERCEPT <- rnorm(SIM_NUM, mean = CH_CR_INTERCEPT, sd = CH_CR_INTERCEPT_SD)

# # 3. coeff_temp for cnidarians and ctenophores
# rt <- read.csv("data/woa_annual_temp_by_depth.csv")
# # ratio of temperature at depth X to surface temperature from the WOA data
# depths <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125, 150, 175, 200)
# temp_values <- adply(depths, 1, function(x){
#   col_index <- which(as.numeric(str_sub(names(rt)[3:104], 2, -1)) == x)
#   coeff_temp <- mean(rt[,col_index+2]/rt$X0, na.rm=T)
#   return(coeff_temp)
# })
# temp_values$X1 <- depths
# names(temp_values) <- c("depths", "coeff_temp")
# 
# # model a "continuous" distribution
# m <- loess(coeff_temp~depths, data=temp_values)
# save(m, file="data/WOA_depth_temperature_regression_model.Rdata")

load(file="data/WOA_depth_temperature_regression_model.Rdata")
# Model the cnid, cten, chor temperatures according to their depth distribution
# # from the JEDI data base pulled out the average distribution of the three grps (above 200 m)
# d <- read.csv("raw_data/JeDI_BCODMO_output_072014.csv", stringsAsFactors=FALSE)
# d$depth <- as.numeric(d$depth)
# # choose just the quantitative and presence/absence data. Does not change much when using quantitative vs quant + P/A. 
# mean(d[d$rank_phylum=="Cnidaria" & d$data_type %in% c("quantitative", "presence/absence") & d$depth <=200, "depth"], na.rm=TRUE)
# # --> 20.30
# mean(d[d$rank_phylum=="Ctenophora" & d$data_type %in% c("quantitative", "presence/absence") & d$depth <=200, "depth"], na.rm=TRUE)
# # --> 20.95
# mean(d[d$rank_phylum=="Chordata" & d$data_type %in% c("quantitative", "presence/absence") & d$depth <=200, "depth"], na.rm=TRUE)
# # --> 64.48
# cnid: 20, cten: 20, chor: 64
# seem to be heavily long tailed, so generated log-normal distributions (mean=1, sd=0.25) to fit to those means
# that log-normal distribution has roughly mean = 3, so multiplier coefficients are: cnid=15, cten=10, chor=20
depths <- rlnorm(SIM_NUM, mean=1, sd=0.25)
sim$coeff_cn_temp <- predict(m, newdata=data.frame(depths=depths*7))
sim$coeff_ct_temp <- predict(m, newdata=data.frame(depths=depths*7))
sim$coeff_ch_temp <- predict(m, newdata=data.frame(depths=depths*20))

# 5. a distribution curve for respiration slopes and intercepts
# NB: CHORDATE respiration rates calculated from biomass values
# CNID_RESP_SLOPE <- 2.5883
# CNID_RESP_SLOPE_SD <- 0.1199 * sqrt(33) * 10
# CNID_RESP_INTERCEPT <- 9.2158
# CNID_RESP_INTERCEPT_SD <- 0.0670 * sqrt(33) * 10
# CTEN_RESP_SLOPE <- 2.207
# CTEN_RESP_SLOPE_SD <- 0.6278 * sqrt(11) * 10
# CTEN_RESP_INTERCEPT <- 9.4416
# CTEN_RESP_INTERCEPT_SD <- 0.3321 * sqrt(11) * 10

CNID_RESP_SLOPE <- 2.5883
CNID_RESP_SLOPE_SD <- 0.1199 * sqrt(33) * 10
CNID_RESP_INTERCEPT <- 9.2158
CNID_RESP_INTERCEPT_SD <- 0.0670 * sqrt(33) * 10
CTEN_RESP_SLOPE <- 2.207
CTEN_RESP_SLOPE_SD <- 0.6278 * sqrt(11) * 10
CTEN_RESP_INTERCEPT <- 9.4416
CTEN_RESP_INTERCEPT_SD <- 0.3321 * sqrt(11) * 10

CHOR_RESP_SLOPE <- 1.05199
CHOR_RESP_SLOPE_SD <- 0.03297 * sqrt(192) * 5
CHOR_RESP_INTERCEPT <- 9.365
CHOR_RESP_INTERCEPT_SD <- 0.01837 * sqrt(192) * 5

sim$CNID_RESP_SLOPE <- rnorm(SIM_NUM, mean = CNID_RESP_SLOPE, sd = CNID_RESP_SLOPE_SD)
sim$CNID_RESP_INTERCEPT <- rnorm(SIM_NUM, mean = CNID_RESP_INTERCEPT, sd = CNID_RESP_INTERCEPT_SD)
sim$CTEN_RESP_SLOPE <- rnorm(SIM_NUM, mean = CTEN_RESP_SLOPE, sd = CTEN_RESP_SLOPE_SD)
sim$CTEN_RESP_INTERCEPT <- rnorm(SIM_NUM, mean = CTEN_RESP_INTERCEPT, sd = CTEN_RESP_INTERCEPT_SD)
sim$CHOR_RESP_SLOPE <- rnorm(SIM_NUM, mean = CHOR_RESP_SLOPE, sd = CHOR_RESP_SLOPE_SD)
sim$CHOR_RESP_INTERCEPT <- rnorm(SIM_NUM, mean = CHOR_RESP_INTERCEPT, sd = CHOR_RESP_INTERCEPT_SD)


# 6. coeff_DOC values
coeff_cn_DOC <- seq(0, 0.2, 0.005)
coeff_ct_DOC <- seq(0.15, 0.45, 0.005)
coeff_ch_DOC <- seq(0.1, 0.4, 0.005)

length <- c(cn=length(coeff_cn_DOC), ct = length(coeff_ct_DOC), ch = length(coeff_ch_DOC))

probs=c(seq(1, floor(length[1]/2)+(length[1]%%2), 1), seq(floor(length[1]/2), 1, -1))
sim$coeff_cn_DOC <- sample(coeff_cn_DOC, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_cn_DOC)

probs=c(seq(1, floor(length[2]/2)+(length[2]%%2), 1), seq(floor(length[2]/2), 1, -1))
sim$coeff_ct_DOC <- sample(coeff_ct_DOC, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_ct_DOC)

probs=c(seq(1, floor(length[3]/2)+(length[3]%%2), 1), seq(floor(length[3]/2), 1, -1))
sim$coeff_ch_DOC <- sample(coeff_ch_DOC, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_ch_DOC)

# 7. predation_EE. Ecotropic efficiency. around 0.45 mean
coeff_cn_predEE <- seq(0.20, 0.55, 0.005)
coeff_ct_predEE <- seq(0.4, 0.9, 0.005)
coeff_ch_predEE <- seq(0.25, 0.55, 0.005)

length <- c(cn=length(coeff_cn_predEE), ct = length(coeff_ct_predEE), ch = length(coeff_ch_predEE))

probs=c(seq(1, floor(length[1]/2)+(length[1]%%2), 1), seq(floor(length[1]/2), 1, -1))
sim$coeff_cn_predEE <- sample(coeff_cn_predEE, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_cn_predEE)

probs=c(seq(1, floor(length[2]/2)+(length[2]%%2), 1), seq(floor(length[2]/2), 1, -1))
sim$coeff_ct_predEE <- sample(coeff_ct_predEE, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_ct_predEE)

probs=c(seq(1, floor(length[3]/2)+(length[3]%%2), 1), seq(floor(length[3]/2), 1, -1))
sim$coeff_ch_predEE <- sample(coeff_ch_predEE, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_ch_predEE)

# 8. coeff_reproduction values
# ctenophores reproduce sexually, release eggs / sperm into water column. Probably releases more DOC with reproduction than cnidarians

coeff_cn_repro <- seq(0, 0.08, 0.002)
coeff_ct_repro <- seq(0.08, 0.16, 0.002)
coeff_ch_repro <- seq(0.08, 0.16, 0.002)

length <- c(cn=length(coeff_cn_repro), ct = length(coeff_ct_repro), ch = length(coeff_ch_repro))

probs=c(seq(1, floor(length[1]/2)+(length[1]%%2), 1), seq(floor(length[1]/2), 1, -1))
sim$coeff_cn_repro <- sample(coeff_cn_repro, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_cn_repro)

probs=c(seq(1, floor(length[2]/2)+(length[2]%%2), 1), seq(floor(length[2]/2), 1, -1))
sim$coeff_ct_repro <- sample(coeff_ct_repro, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_ct_repro)

probs=c(seq(1, floor(length[3]/2)+(length[3]%%2), 1), seq(floor(length[3]/2), 1, -1))
sim$coeff_ch_repro <- sample(coeff_ch_repro, SIM_NUM, replace=TRUE, prob=probs)
#hist(sim$coeff_ch_repro)

# REBALANCE ROUTING
cn_routing <- sim[,c("coeff_cn_predEE","coeff_cn_DOC","coeff_cn_repro")]
ct_routing <- sim[,c("coeff_ct_predEE","coeff_ct_DOC","coeff_ct_repro")]
ch_routing <- sim[,c("coeff_ch_predEE","coeff_ch_DOC","coeff_ch_repro")]

rebalance_routing <- function(data, colnames=c("pred","DOC","repro")){
  data=as.data.frame(data)
  names(data) = colnames
  # step 1, reduce reproduction, prioritize DOC and predation
  rsd1 = rowSums(data[,c("pred","DOC")])
  
  d1gt1 = which(rsd1 >= 1)
  if(length(d1gt1) > 0){
    data[d1gt1, "repro"] = 0
    data[d1gt1, "pred"] = data[d1gt1, "pred"] / (rsd1[d1gt1]+0.002)
    data[d1gt1, "DOC"] = data[d1gt1, "DOC"] / (rsd1[d1gt1]+0.002)
  }

  # step 2, reduce predation 
  rsd2 = rowSums(data)
  d2gt1 = which(rsd2 > 1)
  if (length(d2gt1) > 0) {
    data[d2gt1,"pred"] = data[d2gt1,"pred"] - (rsd2[d2gt1]-1+0.002)
  }
  
  return(data)
}

sim[,c("coeff_cn_predEE","coeff_cn_DOC","coeff_cn_repro")] = rebalance_routing(cn_routing)
sim[,c("coeff_ct_predEE","coeff_ct_DOC","coeff_ct_repro")] = rebalance_routing(ct_routing)
sim[,c("coeff_ch_predEE","coeff_ch_DOC","coeff_ch_repro")] = rebalance_routing(ch_routing)

cn_falls = 1- rowSums(sim[,c("coeff_cn_predEE","coeff_cn_DOC","coeff_cn_repro")])
ct_falls = 1- rowSums(sim[,c("coeff_ct_predEE","coeff_ct_DOC","coeff_ct_repro")])
ch_falls = 1- rowSums(sim[,c("coeff_ch_predEE","coeff_ch_DOC","coeff_ch_repro")])

mean(cn_falls)
mean(ct_falls)
mean(ch_falls)

# 11. assimilation efficiency
# log-normally distributed. mean of ~55% but skewed to the lower range
ae_chor <- rnorm(SIM_NUM+10000, mean=0.5, sd=0.2) # mean = 2.84
ae_chor <- ae_chor[ae_chor <= 0.9 & ae_chor >= 0.1]
sim$ae_chor <- sample(ae_chor, SIM_NUM, replace=FALSE)

ae_cn_ct <- rnorm(SIM_NUM+10000, 0.8, 0.1)
sim$ae_cn_ct <- sample(ae_cn_ct[which(ae_cn_ct < 1)], SIM_NUM, replace=FALSE)


write.csv(sim, "data/MonteCarloSimulations.csv", row.names=F)

