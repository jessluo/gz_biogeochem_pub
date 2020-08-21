#!/usr/bin/env Rscript
#
#     Generate sinking/export flux
#       plots for all three cases put together
#     requires outputs from analysis-generate_ensemble_means.R
#
#
# ------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("gridExtra")
library("stringr")
library("viridis")
library("RColorBrewer")

source("functions.R")
Spectral <- brewer.pal(8, "Spectral")

# load biomes
load('data/biomes_chl_mixedlayer.Rdata')
N_BIOME=4

# import area/volume
av <- read.csv("data/surfaceArea_volume_1-deg.csv")
av <- join(av,biomes_df, by=c("lat","lon"))
av$depth.m = av$vol / av$area.m2 #vol in m3
av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(vol=sum(x$vol)))})
names(av_biomesum) <- c("biome", "vol")
av_biomesum <- rbind(av_biomesum, data.frame(biome='TOTAL', vol=sum(av_biomesum[av_biomesum$biome != 'TOTAL','vol'])))
av_biomesum$volfrac = av_biomesum$vol/av_biomesum[av_biomesum$biome=="TOTAL","vol"]

av_biomedepth <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(depth.m=mean(x$depth.m)))})


# Import Export to depth data for all three biomes
exportres = data.frame()
uores = data.frame()

for (C_ID in 1:3){
  CASES = c("0-baseline","1-high_biomass","2-low_biomass")
  CASE = CASES[C_ID]
  
  load(file=str_c("data/",CASE,"/baseline/allexport.Rdata")) # units are all g C m^-3 y^-1
  load(file=str_c("data/",CASE,"/baseline/upperocean_ens_mean.Rdata"))
  
  allexport$case <- CASE
  ens_mean$case <- CASE
  
  exportres <- rbind(exportres, allexport)
  uores <- rbind(uores, ens_mean)
}

exportres$case = factor(exportres$case, levels=c("2-low_biomass", "0-baseline", "1-high_biomass"))
exportres <- join(exportres, biomes_df, by=c("lat","lon"))

## { Carcass fraction --------------------------------

carc = exportres
carc = ddply(carc[carc$taxon=="sum",], ~case+taxon+biome, function(x){
  data.frame(export=mean(x$export), export_cf=mean(x$export_cf),
             seq_d=mean(x$seq_d), seq_d_cf=mean(x$seq_d_cf), 
             seafloor=mean(x$seafloor), seafloor_cf=mean(x$seafloor_cf))
})
carc$export_cf_frac = carc$export_cf / carc$export
carc$seq_d_cf_frac = carc$seq_d_cf / carc$seq_d
carc$seafloor_cf_frac = carc$seafloor_cf / carc$seafloor
carc = carc[,c("case","biome","export_cf_frac","seq_d_cf_frac","seafloor_cf_frac")]

carcm = melt(carc, id.vars = c("case","biome"))
carcm = join(carcm, av_biomesum)
carcm$value_scaled=carcm$value * carcm$volfrac
carcm = ddply(carcm, ~case+variable, function(x){sum(x$value_scaled)})

write.csv(carcm, "data/carcass_fraction_by_depth.csv", row.names=FALSE)

## seafloor flux at depths > 50 m
deepseafl = join(exportres, av)
deepseafl = deepseafl[deepseafl$depth > 50,]
deepseafl = deepseafl[deepseafl$taxon=="sum",c("lat","lon","biome","case","seafloor")]
deepseafl = ddply(deepseafl, ~case+biome, function(x){data.frame(seafloor=mean(x$seafloor))})
deepseafl = join(deepseafl, av_biomesum)
deepseafl$seafloor_vol = deepseafl$seafloor * deepseafl$vol
deepseafl = ddply(deepseafl, ~case, function(x){data.frame(seafloor_vol=sum(x$seafloor_vol))})
# }
## { Barplots ---------------------------------------
# create barplots
col_names = names(exportres)[which(names(exportres) %ni% c("lat","lon","taxon","case","biome"))]

# take the mean within biomes
plotdf <- ddply(exportres, ~case+taxon+biome, function(x){
  tmpd = apply(x[,col_names], 2, mean)
})


plotdf = melt(plotdf, id.vars = c("case","biome","taxon"), measure.vars = c("export_cf","export_eg"), variable.name = "loc")
plotdf = join(plotdf, av_biomedepth)
plotdf$value = plotdf$value * plotdf$depth.m

# set levels
plotdf$taxon = factor(plotdf$taxon, levels=c("cnid", "cten", "chor", "sum"))
plotdf$biome = factor(plotdf$biome, levels=c("COAST", "HCPS", "HCSS", "LC"))

# New facet label names 
expLabels <- c("Carcasses", "Fecal Matter")
names(expLabels) <- unique(plotdf$loc) #c("export_cf", "export_eg")
col <- c("#66c2a5", "#fc8d62", "#8da0cb")

p <- ggplot(plotdf[plotdf$taxon != "sum",]) + geom_bar(aes(x=case, y=value, fill=taxon), stat = 'identity') + 
  facet_grid(~biome+loc, labeller=labeller(loc=expLabels)) + theme_bw() + labs(y='Flux (g C m-2)',x="") + 
  scale_fill_manual("Taxon", values = col, labels=c("Cnidarians","Ctenophores","Tunicates")) +
  scale_x_discrete(labels=c("Low","Base", "High")) + scale_y_continuous(limits=c(0,11)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggsave('plots/barchart_export_depth_1deg_allcases.pdf', p, height=4, width=10)


# seafloor
plotdf <- ddply(exportres, ~case+taxon+biome, function(x){
  tmpd = apply(x[,col_names], 2, mean)
})

plotdf = melt(plotdf, id.vars = c("case","biome","taxon"), measure.vars = c("seafloor_cf","seafloor_eg"), variable.name = "loc")
plotdf = join(plotdf, av_biomedepth)
plotdf$value = plotdf$value * plotdf$depth.m

# set levels
plotdf$taxon = factor(plotdf$taxon, levels=c("cnid", "cten", "chor", "sum"))
plotdf$biome = factor(plotdf$biome, levels=c("COAST", "HCPS", "HCSS", "LC"))

# New facet label names 
expLabels <- c("Carcasses", "Fecal Matter")
names(expLabels) <- unique(plotdf$loc) #c("export_cf", "export_eg")
col <- c("#66c2a5", "#fc8d62", "#8da0cb")

p <- ggplot(plotdf[plotdf$taxon != "sum",]) + geom_bar(aes(x=case, y=value, fill=taxon), stat = 'identity') + 
  facet_grid(~biome+loc, labeller=labeller(loc=expLabels)) + theme_bw() + labs(y='Flux (g C m-2)',x="") + 
  scale_fill_manual("Taxon", values = col, labels=c("Cnidarians","Ctenophores","Tunicates")) +
  scale_x_discrete(labels=c("Low","Base", "High")) + scale_y_continuous(limits=c(0,11)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggsave('plots/barchart_seafloor_1deg_allcases.pdf', p, height=4, width=10)

# }


## { Latitudinal averages ----------------------------
exportres <- join(exportres, av[,c("lat","lon","depth.m")])
uores <- join(uores, av[,c("lat","lon","depth.m")])


# set factors as such so you don't get high biomass line over the baseline line in the plot
# just for aesthetics
exportres$case = factor(exportres$case, levels=c("2-low_biomass","1-high_biomass","0-baseline"))
uores$case = factor(uores$case, levels=c("2-low_biomass","1-high_biomass","0-baseline"))

## e-ratio by latitude
lat_eratio = data.frame()
for (C_ID in 1:3){
  CASES = c("0-baseline","1-high_biomass","2-low_biomass")
  CASE = CASES[C_ID]
  
  data = read.csv(str_c("data/",CASE,"/latitudinal_mean_seafloor.csv"), as.is=TRUE)
  data$case=CASE
  
  lat_eratio = rbind(lat_eratio, data)
}
lat_eratio$case = factor(lat_eratio$case, levels=c("2-low_biomass","1-high_biomass","0-baseline"))

p1 <- ggplot(lat_eratio) + geom_line(aes(x=lat, y=mean_seafloor_flux, color=case)) + labs(x="", y="") + theme_bw() + 
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), breaks=c("2-low_biomass", "0-baseline", "1-high_biomass")) + 
  coord_flip() + theme(legend.position = c(0.6, 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = str_c("plots/lat_seafloor_flux_1deg.pdf"), p1, width=2.5, height=4)

p2 <- ggplot(lat_eratio) + geom_line(aes(x=lat, y=e_ratio, color=case)) + labs(x="", y="") + theme_bw() + 
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), breaks=c("2-low_biomass", "0-baseline", "1-high_biomass")) + 
  coord_flip() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = str_c("plots/lat_seafloor_e_ratio_1deg.pdf"), p2, width=2.5, height=4)


# supplemental plots: lat averages by taxa & biome 

lat_averages <- ddply(exportres[exportres$taxon!="sum",], ~lat+case+taxon+biome, function(x){
  export_cf_flux.gCm2 = mean(x$export_cf * x$depth, na.rm=TRUE)
  export_eg_flux.gCm2 = mean(x$export_eg * x$depth, na.rm=TRUE)
  seq_d_cf_flux.gCm2 = mean(x$seq_d_cf * x$depth, na.rm=TRUE)
  seq_d_eg_flux.gCm2 = mean(x$seq_d_eg * x$depth, na.rm=TRUE)
  seafloor_cf_flux.gCm2 = mean(x$seafloor_cf * x$depth, na.rm=TRUE)
  seafloor_eg_flux.gCm2 = mean(x$seafloor_eg * x$depth, na.rm=TRUE)
  return(data.frame(export_cf_flux.gCm2, export_eg_flux.gCm2, seq_d_cf_flux.gCm2, 
                    seq_d_eg_flux.gCm2, seafloor_cf_flux.gCm2, seafloor_eg_flux.gCm2))
})

## plotting manipulations
lat_averages <- melt(lat_averages, id.vars = c("lat","case","taxon","biome"))
lat_averages$depth = str_sub(lat_averages$variable, 1, -14)
lat_averages$depth = factor(lat_averages$depth, levels=c("export", "seq_d", "seafloor"))
lat_averages$type = str_sub(lat_averages$variable, -12, -1)

# facet labels
expLabels <- c("Carcasses", "Fecal Matter")
names(expLabels) <- c("cf_flux.gCm2", "eg_flux.gCm2")
expLabels1 <- c("100 m", "1000 m", "Seafloor")
names(expLabels1) <- c("export", "seq_d", "seafloor")

p1 = ggplot(lat_averages[lat_averages$taxon=="cnid",]) + 
  geom_line(aes(x=lat, y=value, color=case), size=0.3) + 
  labs(x="", y="POC Flux (g C m-2)", title="Cnidarians - POC Export") + theme_bw() + 
  facet_grid(biome~depth+type, labeller=labeller(type=expLabels, depth=expLabels1)) +
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), 
                     breaks=c("2-low_biomass", "0-baseline", "1-high_biomass"), labels=c("Low", "Base", "High")) + 
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/cnid_latitudinal_means_biome.pdf", p1, height=10, width=7.5)

p2 = ggplot(lat_averages[lat_averages$taxon=="cten",]) + 
  geom_line(aes(x=lat, y=value, color=case), size=0.3) + 
  labs(x="", y="POC Flux (g C m-2)", title="Ctenophores - POC Export") + theme_bw() + 
  facet_grid(biome~depth+type, labeller=labeller(type=expLabels, depth=expLabels1)) +
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), 
                     breaks=c("2-low_biomass", "0-baseline", "1-high_biomass"), labels=c("Low", "Base", "High")) + 
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/cten_latitudinal_means_biome.pdf", p2, height=10, width=7.5)

p3 = ggplot(lat_averages[lat_averages$taxon=="chor",]) + 
  geom_line(aes(x=lat, y=value, color=case), size=0.3) + 
  labs(x="", y="POC Flux (g C m-2)", title="Tunicates - POC Export") + theme_bw() + 
  facet_grid(biome~depth+type, labeller=labeller(type=expLabels, depth=expLabels1)) +
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), 
                     breaks=c("2-low_biomass", "0-baseline", "1-high_biomass"), labels=c("Low", "Base", "High")) + 
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/chor_latitudinal_means_biome.pdf", p3, height=10, width=7.5)


uo_lat_averages <- ddply(uores[uores$taxon!="sum",], ~lat+case+taxon+biome, function(x){
  totingestion.gCm2 = mean((x$ingestion + x$egestion) * x$depth, na.rm=TRUE)
  egestion.gCm2 = mean(x$egestion * x$depth, na.rm=TRUE)
  respiration.gCm2 = mean(x$respiration * x$depth, na.rm=TRUE)
  alldoc.gCm2 = mean((x$DOC + x$reproduction) * x$depth, na.rm=TRUE)
  predation.gCm2 = mean(x$predation * x$depth, na.rm=TRUE)
  
  x[which(x$flux < 0),] = 0
  flux.gCm2 = mean(x$flux * x$depth, na.rm=TRUE)
  return(data.frame(totingestion.gCm2, egestion.gCm2, respiration.gCm2, alldoc.gCm2, 
                    predation.gCm2, flux.gCm2))
})

uo_lat_averages <- melt(uo_lat_averages, id.vars = c("lat","case","taxon","biome"))
# facet labels
expLabels2 <- c("Ingestion", "Egestion", "Respiration", "Excretion", "Predation", "Mortality")
names(expLabels2) <- c("totingestion.gCm2", "egestion.gCm2", "respiration.gCm2", "alldoc.gCm2", 
                      "predation.gCm2", "flux.gCm2")

p1 <- ggplot(uo_lat_averages[uo_lat_averages$taxon=="cnid",]) + 
  geom_line(aes(x=lat, y=value, color=case), size=0.3) + 
  labs(x="", y="Carbon Flux (g C m-2)", title="Cnidarians - Upper Ocean Fluxes") + theme_bw() + 
  facet_grid(biome~variable, labeller=labeller(variable=expLabels2)) +
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) +
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), 
                     breaks=c("2-low_biomass", "0-baseline", "1-high_biomass"), labels=c("Low", "Base", "High")) + 
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/cnid_UO_latitudinal_means_biome.pdf", p1, height=10, width=7.5)

p2 = ggplot(uo_lat_averages[uo_lat_averages$taxon=="cten",]) + 
  geom_line(aes(x=lat, y=value, color=case), size=0.3) + 
  labs(x="", y="Carbon Flux (g C m-2)", title="Ctenophores - Upper Ocean Fluxes") + theme_bw() + 
  facet_grid(biome~variable, labeller=labeller(variable=expLabels2)) +
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), 
                     breaks=c("2-low_biomass", "0-baseline", "1-high_biomass"), labels=c("Low", "Base", "High")) + 
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/cten_UO_latitudinal_means_biome.pdf", p2, height=10, width=7.5)

p3 = ggplot(uo_lat_averages[uo_lat_averages$taxon=="chor",]) + 
  geom_line(aes(x=lat, y=value, color=case), size=0.3) + 
  labs(x="", y="Carbon Flux (g C m-2)", title="Tunicates - Upper Ocean Fluxes") + theme_bw() + 
  facet_grid(biome~variable, labeller=labeller(variable=expLabels2)) +
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + 
  scale_y_continuous(trans="sqrt", breaks=c(10,100,300), limits=c(0,500)) +
  scale_color_manual("Case", values = c('#ef8a62','#67a9cf','#030303'), 
                     breaks=c("2-low_biomass", "0-baseline", "1-high_biomass"), labels=c("Low", "Base", "High")) + 
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/chor_UO_latitudinal_means_biome.pdf", p3, height=10, width=7.5)


# }
