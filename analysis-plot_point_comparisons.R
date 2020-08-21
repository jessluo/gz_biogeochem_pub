#!/usr/bin/env Rscript
#
#
#   Plot single point comparisons
#
#
# -------------------------------------------

# load libraries ----------------------------
library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")
source('functions.R')

CASES = c("0-baseline","1-high_biomass","2-low_biomass")

CONF_INT=TRUE
QUANTILES=FALSE

## { Read in single point model data ---------------------
singleptinfo <- read.csv("data/site_comparisons_info.csv")
singleptinfo$biome = factor(singleptinfo$biome, levels=c("COAST","LC","HCPS","HCSS"))
singleptinfo$site = as.character(singleptinfo$site)
singleptinfo = arrange(singleptinfo, taxa, type, biome, site)

plotarrange = c("atl_wiebe_cf", "bats_cf", "sargasso_sea_cf", "ivory_coast", "nw_med", "stn_M_carcass", "tasman_sea", 
                "atl_madin_eg", "atl_wiebe_eg", "atl_caron_eg", "bats_eg", "sargasso_sea_eg", "ne_pacific", 
                "stn_M_fecalpellets", "southern_ocean_atl", "west_antarctica", "lazarev_sea", 
                "chesapeake_bay", "norway_fjord", "sea_of_japan", "gulf_of_oman")

d = data.frame()

for (CASE in CASES){
  singlept = read.csv(str_c("data/",CASE,"/site_comparisons_allpoints.csv"))
  tmp = melt(singlept, measure.vars = names(singlept), variable.name = "site", value.name="flux_gC.m2.y")
  tmp = tmp[complete.cases(tmp),]
  tmp = join(tmp, singleptinfo, by="site")
  tmp$case=CASE
  d = rbind(d, tmp)
}

# } 

## { Calculate min/max and quantiles / 95% confidence intervals ------
minmax = ddply(d, ~site, function(x){
  min=min(x$flux_gC.m2.y)
  max=max(x$flux_gC.m2.y)
  return(data.frame(min,max))
})

print (minmax)

# site          min        max
# 1      chesapeake_bay 4.822571e-05  0.5283815
# 2        norway_fjord 3.568640e-05  1.2149753
# 3        sea_of_japan 4.489364e-04  2.9590836
# 4        gulf_of_oman 0.000000e+00  1.6321514
# 5              nw_med 1.203817e-03  4.4289016
# 6         ivory_coast 0.000000e+00  2.1871902
# 7          tasman_sea 3.997744e-04  4.6254443
# 8       stn_M_carcass 2.210351e-04  2.4946805
# 9  stn_M_fecalpellets 2.258508e-04 12.3607571
# 10    west_antarctica 3.979841e-03 28.5092138
# 11        lazarev_sea 3.498851e-03 12.1954666
# 12 southern_ocean_atl 4.108478e-05 14.3664228
# 13         ne_pacific 3.599207e-08 11.7162756
# 14       atl_wiebe_cf 1.602874e-03 10.3926439
# 15       atl_wiebe_eg 8.696498e-03 41.4571525
# 16       atl_madin_eg 7.493082e-02 57.6746016
# 17       atl_caron_eg 2.786393e-02 55.7661698
# 18            bats_cf 6.786471e-05  2.9045361
# 19            bats_eg 1.671171e-10 12.9024902
# 20    sargasso_sea_cf 6.191009e-04  4.2427273
# 21    sargasso_sea_eg 9.303376e-03 18.2626819

conf_int = ddply(d, ~site, function(x){
  if(CONF_INT)  res = return_bci(x$flux_gC.m2.y)
  if(QUANTILES)  res = return_quant(x$flux_gC.m2.y)
  return(data.frame(lower=res[1],upper=res[2]))
}, .progress="text")

print(conf_int)

# site        lower       upper
# 1      chesapeake_bay 1.434737e-02  0.15595452
# 2        norway_fjord 8.672355e-03  0.17529128
# 3        sea_of_japan 7.707006e-02  1.14521397
# 4        gulf_of_oman 0.000000e+00  0.16298465
# 5              nw_med 1.249342e-01  0.81696233
# 6         ivory_coast 8.825412e-02  0.56382430
# 7          tasman_sea 1.289690e-01  0.85735771
# 8       stn_M_carcass 2.417814e-02  0.18244621
# 9  stn_M_fecalpellets 1.725991e-02  0.72177883
# 10    west_antarctica 3.000739e-01  4.15399984
# 11        lazarev_sea 3.676372e-01  2.84025169
# 12 southern_ocean_atl 4.320564e-03  0.06390801
# 13         ne_pacific 1.472996e-05  0.06890449
# 14       atl_wiebe_cf 1.105694e-01  0.99650659
# 15       atl_wiebe_eg 5.640951e-01  4.28939786
# 16       atl_madin_eg 1.960454e+00 11.17073039
# 17       atl_caron_eg 7.804380e-01  5.56083676
# 18            bats_cf 1.953406e-02  0.24771915
# 19            bats_eg 1.338830e-08  0.87414170
# 20    sargasso_sea_cf 5.590469e-02  0.49126110
# 21    sargasso_sea_eg 3.178218e-01  2.15638743

res = join(minmax, conf_int)
write.csv(res, "data/site_comparisons.csv", row.names=FALSE)
# } 

## { Plot up the comparions with data ---------------------
tmp <- unique(d[,c("site", "taxa", "type", "case", "biome")])
tmp$case <- "observations"
tmp$flux_gC.m2.y <- 0
d <- rbind(d, tmp)

obs = rbind(data.frame(site='atl_madin_eg', min=1.8, max=33),
            data.frame(site='atl_caron_eg', min=3.65E-3, max=0.0255),
            data.frame(site='west_antarctica', min=0.36, max=1.8),
            data.frame(site='lazarev_sea', min=7.92, max=7.92),
            data.frame(site='southern_ocean_atl', min=0.031, max=2.1),
            data.frame(site='atl_wiebe_eg', min=3.1, max=50),
            data.frame(site='atl_wiebe_cf', min=1.31, max=1.31),
            data.frame(site='ne_pacific', min=3.17, max=3.17),
            data.frame(site='stn_M_fecalpellets', min=0, max=1.75),
            data.frame(site='stn_M_carcass', min=0, max=0.1825),
            data.frame(site='ivory_coast', min=1.0, max=22.0),
            data.frame(site='nw_med', min=6E-4, max=0.016),
            data.frame(site='tasman_sea', min=16, max=16),
            data.frame(site='bats_eg', min=0, max=1.726),
            data.frame(site='bats_cf', min=0.0035, max=0.0425),
            data.frame(site='sargasso_sea_eg', min=0, max=1.726),
            data.frame(site='sargasso_sea_cf', min=0.0035, max=0.0425),
            data.frame(site='gulf_of_oman', min=1.5, max=78),
            data.frame(site='norway_fjord', min=6.36, max=11.4),
            data.frame(site='chesapeake_bay', min=1.7E-4, max=0.112),
            data.frame(site='sea_of_japan', min=0.562, max=14.32))

obs <- join(obs, unique(d[,c("site", "taxa", "type")]))
obsplot = rbind(data.frame(c(obs, case='observations')),
                data.frame(c(obs, case='0-baseline')), 
                data.frame(c(obs, case='1-high_biomass')),
                data.frame(c(obs, case='2-low_biomass')))
obsplot[obsplot$case != "observations", c("min", "max")] <- 0

# set levels for plotting
obsplot$case = factor(obsplot$case, levels=c( "2-low_biomass","0-baseline", "1-high_biomass", "observations"))
d$case = factor(d$case, levels=c( "2-low_biomass","0-baseline", "1-high_biomass", "observations"))

obsplot$site = factor(obsplot$site, levels=plotarrange)
d$site = factor(d$site, levels=plotarrange)

p <- ggplot() + geom_boxplot(aes(x=site, y=flux_gC.m2.y, fill=case, linetype=case), outlier.alpha = 0.2, data=d) + 
  geom_errorbar(aes(x=site, ymin=min, ymax=max, color=case, alpha=case), width=0.7, position=position_dodge(1.25), size=1, data=obsplot)+
  facet_grid(.~taxa+type, scales = "free_x", space = "free_x") + 
  scale_fill_manual(values=c("#4daf4a", "#377eb8", "#984ea3", "white")) + 
  scale_linetype_manual(values=c("solid", "solid", "solid", "blank")) +
  scale_color_manual(values=c("#4daf4a", "#377eb8", "#984ea3", "#e41a1c")) + scale_alpha_manual(values=c(0,0,0,1)) + 
  scale_y_sqrt(expression(paste('POC Flux (g C ',m^-2,' ',y^{-1},')')), breaks=c(1, 4, 9, 16, 25, 36, 49, 64), limits=c(0,78)) + xlab('') +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=0.2, vjust=0.2))

ggsave(filename='plots/site_comparisons.pdf', p, height=6, width=12)

# } 

