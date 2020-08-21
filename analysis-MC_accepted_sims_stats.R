#!/usr/bin/env Rscript
#
#
#    Extract Monte Carlo mean accepted values
#     and number of accepted values
#
# -------------------------------------------


# reading packages
library(plyr)
library(reshape2)
library(stringr)

# read original MC sims
sim <- read.csv("data/MonteCarloSimulations.csv")
pre_mean = data.frame(value=apply(sim, 2, mean))
pre_mean$variable = row.names(pre_mean)
pre_mean$case = "orig"

sim_cn <- sim_ct <- sim_ch <- data.frame()

n_accepted=as.data.frame(array(NA,dim = c(3,4)))
n_accepted[,1]=c('cnid','cten','chor')

for (C_ID in 1:3){
  CASES = c("0-baseline","1-high_biomass","2-low_biomass")
  CASE = CASES[C_ID]

  load(str_c("data/",CASE,"/baseline/simulation_posteriors_cnid.Rdata"))
  load(str_c("data/",CASE,"/baseline/simulation_posteriors_cten.Rdata"))
  load(str_c("data/",CASE,"/baseline/simulation_posteriors_chor.Rdata"))
  
  sim_cn_case <- as.data.frame(t(colMeans(cnid_posterior)))
  sim_ct_case <- as.data.frame(t(colMeans(cten_posterior)))
  sim_ch_case <- as.data.frame(t(colMeans(chor_posterior)))
  
  sim_cn_case$case = CASE
  sim_ct_case$case = CASE
  sim_ch_case$case = CASE
  
  sim_cn <- rbind(sim_cn, sim_cn_case)
  sim_ct <- rbind(sim_ct, sim_ct_case)
  sim_ch <- rbind(sim_ch, sim_ch_case)
  
  n_cn=dim(cnid_posterior)[1]
  n_ct=dim(cten_posterior)[1]
  n_ch=dim(chor_posterior)[1]
  
  n_accepted[,C_ID+1] = c(n_cn, n_ct, n_ch)
  
}

sim_cn = data.frame(c(melt(sim_cn, id.vars = 'case'), taxa='cnid'))
sim_ct = data.frame(c(melt(sim_ct, id.vars = 'case'), taxa='cten'))
sim_ch = data.frame(c(melt(sim_ch, id.vars = 'case'), taxa='chor'))

meansims = rbind(sim_cn, sim_ct, sim_ch)

orig=join(meansims[meansims$case=="0-baseline",c("variable","taxa")],pre_mean)
meansims = rbind(meansims, orig)
#write.csv(meansims, 'data/MonteCarlo_AcceptedSims_mean.csv', row.names=FALSE)

# dcast test
meansims$case = factor(meansims$case, levels=c("orig","2-low_biomass","0-baseline","1-high_biomass"))
dc = dcast(meansims, variable~taxa+case)
write.csv(dc, 'data/MonteCarlo_AcceptedSims_mean.csv', row.names=FALSE)

# number of accepted
names(n_accepted) = c("taxa", CASES)
write.csv(n_accepted, 'data/MonteCarlo_NumberofAcceptedSims.csv', row.names=FALSE)

##



  