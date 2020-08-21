#!/usr/bin/env Rscript
#
#
#       Analysis/Plotting to evaluate GZ POC flux
#         & compare with COBALT sinkng POC flux
#
#
# ---------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("gridExtra")
library("stringr")
library("viridis")
library("RColorBrewer")

source("functions.R")
Spectral <- brewer.pal(8, "Spectral")

args <- commandArgs(trailingOnly=TRUE)
C_ID <- as.numeric(args[1])
# C_ID = 1
CASES = c("0-baseline","1-high_biomass","2-low_biomass")
CASE = CASES[C_ID]

## { Read and load in data -----------------------------------------------------------
# import functions and constants
OCEANS_SURFACE_AREA <- 3.619e+14
OCEANS_SURFACE_VOL <- 6.708e+16

# laod cobalt export flux
load(file="data/gfdl_POCflux_100m_1000m_seafloor.Rdata")

# load coastline
coastline.world <- read.csv("data/gshhg_world_c.csv")

# convert to standard units of g C m^-2 y-1
cobalt_det[, c("det_100", "det_1000", "det_btm")] <- cobalt_det[, c("det_100", "det_1000", "det_btm")] / 1000 * 365

# load gridded results data
load(file=str_c("data/",CASE,"/baseline/allexport.Rdata")) # note that these units are all g C m^-3 y^-1

# load npp, units of mg C m^-3 d^-1
load("data/gz_model_inputs/pp.Rdata")
names(pp)[3] <- "npp"
pp$npp <- pp$npp / 1000 * 365 # convert to g C m^-3 y^-1

# load biomes
load('data/biomes_chl_mixedlayer.Rdata')

# add in area and volume
# import area/volume
av <- read.csv("data/surfaceArea_volume_1-deg.csv")
av$depth <- av$vol / av$area.m2
av$depth100 <- av$depth
av$depth100 <- ifelse(av$depth100 > 100, 100, av$depth100)
av$vol100 <- av$depth100 * av$area.m2
av <- join(av,biomes_df, by=c("lat","lon"))
av <- av[,c("lat","lon","area.m2","vol","vol100","depth","depth100","biome")]

av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(c(sum(x$vol), sum(x$vol100), mean(x$depth), mean(x$depth100)))})
names(av_biomesum) <- c("biome", "vol",  "vol100", "depth", "depth100")

allexport <- join(join(allexport, pp), av)
allexport[is.na(allexport$biome),"depth"] <- 1
allexport[is.na(allexport$biome),"area.m2"] <- 1e8
allexport[is.na(allexport$biome),"biome"] <- 'COAST'
# } 

## { Pull out summary stats for 1 degree ------------------------------------------------------

# test out the sum of export flux 
mean(cobalt_det$det_100, na.rm=TRUE) * OCEANS_SURFACE_AREA
# --> 6.488302e+15 g C y^-1

# NB: allexport variable does not allow for negative fluxes to be propagated downwards
exportcols <- c("lat", "lon", "taxon", "export", "export_cf", "export_eg", 
                "seq_d", "seq_d_cf", "seq_d_eg", "seafloor", "seafloor_cf", "seafloor_eg", 
                "area.m2", "depth", "depth100", "vol", "vol100")

exportres <- allexport[,c(exportcols, "biome", "npp")]
exportres$e_ratio <- ((exportres$export * exportres$depth100) / (exportres$npp * exportres$depth))

# add in GFDL/COBALT model export flux
exportres <- join(exportres, cobalt_det)

# set up a data frame for the comparison ratios
ratios <- exportres[,c("lat", "lon", "taxon", "e_ratio")]

# calculate the ratios between jelly export and model (COBALT) export at different depths
ratios$jelly_v_cobalt_100 <- (exportres$export * exportres$depth) / exportres$det_100
ratios$jelly_v_cobalt_1000 <- (exportres$seq_d * exportres$depth) / exportres$det_1000
ratios$jelly_v_cobalt_seafl <- (exportres$seafloor * exportres$depth) / exportres$det_btm
  
# 1. PIECHART
# calculate fecal pellet and carcass flux past 100 m
FPC = c("fecal_pellets","carcasses")
plotdf <- ddply(exportres, ~taxon+biome, function(x){
  area <- mean(x$area.m2,na.rm=TRUE)
  fecal_pellets = mean(x$seafloor_eg, na.rm=TRUE) * area
  carcasses = mean(x$seafloor_cf, na.rm=TRUE) * area
  return(data.frame(fecal_pellets,carcasses))
})
# get group-specific fraction of sum
plotdf <- ddply(plotdf, ~biome, function(x){
  for (taxon in c("cnid","cten","chor")){
    x[x$taxon==taxon,FPC] <- x[x$taxon==taxon,FPC] / x[x$taxon=="sum",FPC]}
  return(x[x$taxon != "sum",])
})

piechart_df <- ddply(plotdf, ~taxon, function(x){return(colSums(x[,FPC]))})
piechart_df <- melt(piechart_df, id.vars="taxon")
piechart_df$taxon <- factor(piechart_df$taxon, levels=c("cnid", "cten", "chor"))

blank_theme <- theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(),
        axis.text=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
p = ggplot(piechart_df, aes(x=1, y=value, fill=taxon, alpha=variable)) + geom_bar(stat="identity") + 
  coord_polar("y") + blank_theme + scale_alpha_discrete("type", range = c(0.5, 0.8))
ggsave(filename = str_c("plots/",CASE,"/export/piechart_1deg.pdf"), p, width=4.5, height=4.5)


# flux to seafloor
plotdf <- exportres[exportres$taxon=="sum", c("lat", "lon", "depth", "biome", "seafloor", "e_ratio")]
blues <- brewer.pal(9, "Blues")


p1 <- ggplot(mapping=aes(x=lon, y=lat)) + 
  geom_point(aes(colour=seafloor * depth), size=0.32, data=plotdf, shape=15) + labs(x="", y="") + 
  scale_x_continuous(breaks=c(-90, 0, 90), expand = c(0,0)) +
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) + 
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                     panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
  scale_color_gradientn("seafloor", colours = c(rev(Spectral), "red"), na.value="grey80", limits=c(0,40))
ggsave(filename = str_c("plots/",CASE,"/export/seafloor_flux_1deg.pdf"), p1, width=8.3, height=3.85)


# latitudional averages
lat_seafloor <- ddply(plotdf, ~biome+lat, function(x){
  seafloor_flux = x$seafloor * x$depth
  mean_seafloor_flux = mean(seafloor_flux, na.rm=TRUE)
  e_ratio = mean(x$e_ratio, na.rm=TRUE)
  return(data.frame(mean_seafloor_flux, e_ratio))
})

lat_seafloor <- ddply(lat_seafloor, ~lat, function(x){
  mean_seafloor_flux = mean(x$mean_seafloor_flux)
  e_ratio = mean(x$e_ratio)
  return(data.frame(mean_seafloor_flux, e_ratio))
})

write.csv(lat_seafloor, file = str_c("data/",CASE,"/latitudinal_mean_seafloor.csv"), row.names=FALSE)

p1 <- ggplot(lat_seafloor) + geom_line(aes(x=lat, y=mean_seafloor_flux)) + labs(x="", y="") + theme_bw() + 
  scale_x_continuous(limits=c(-90,90), breaks = c(-60, -30, 0, 30, 60)) + coord_flip()
ggsave(filename = str_c("plots/",CASE,"/export/lat_seafloor_flux_1deg.pdf"), p1, width=2.5, height=4)

# e-ratio
p1 <- ggplot(mapping=aes(x=lon, y=lat)) + 
  geom_point(aes(colour=e_ratio), size=0.32, data=ratios[ratios$taxon=="sum",], shape=15) + labs(x="", y="") + 
  scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) + 
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) + 
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                     panel.grid.minor=element_blank(), panel.grid.major=element_blank())

ggsave(filename = str_c("plots/",CASE,"/export/jellyC_export_ratio_nolegend_1deg.pdf"), 
       plot = p1 + scale_color_gradientn("e-ratio", colours = plasma(20), na.value="grey80", limits=c(0,0.45), guide="none"), width=7, height=4)
ggsave(filename = str_c("plots/",CASE,"/export/jellyC_export_ratio_1deg.pdf"), 
       plot = p1 + scale_color_gradientn("e-ratio", colours = plasma(20), na.value="grey80", limits=c(0,0.45)), width=8.5, height=4)


plotdf <- melt(ratios[ratios$taxon=="sum", which(names(ratios) != "taxon")], id.vars=c("lat", "lon"))

p <- ggplot(mapping=aes(x=lon, y=lat)) + 
  geom_point(aes(colour=value), size=0.35, data=plotdf[plotdf$variable != "e_ratio",], shape=15) + labs(x="", y="") + 
  scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) + 
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) + 
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                     panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
  scale_color_gradientn("", colours = c("white", viridis(12),rev(heat.colors(16))[7:16]), na.value="grey80", limits=c(0,48)) + 
  facet_grid(variable~., labeller=as_labeller(c(`jelly_v_cobalt_100`="100 m", `jelly_v_cobalt_1000`="1000 m", `jelly_v_cobalt_seafl`="seafloor"))) 

ggsave(filename = str_c("plots/",CASE,"/export/jelly_v_cobalt_all_1deg.pdf"), plot = p, width=8.425, height=11.275)

# }

## { Transfer efficiency to seafloor ---------------------------------------------------------------------
# calculate benthic transfer efficiency
cobalt_det$trans_eff <- cobalt_det$det_btm / cobalt_det$det_100
cobalt_det$ts_1000 <- cobalt_det$det_1000 / cobalt_det$det_100

ratios$j_trans_eff <- exportres$seafloor / exportres$export

ratios <- join(ratios, cobalt_det[,c("lat", "lon", "trans_eff")])

# calculate a factor difference
ratios$trans_eff_fctdiff <- (ratios$j_trans_eff / ratios$trans_eff)


p1 <- ggplot(mapping=aes(x=lon, y=lat)) + 
  geom_point(aes(colour=j_trans_eff), size=0.35, data=ratios[ratios$taxon=="sum",], shape=15) + labs(x="", y="") + 
  scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) + 
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) + 
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                    panel.grid.minor=element_blank(), panel.grid.major=element_blank())

ggsave(filename = str_c("plots/",CASE,"/export/jelly_benthic_transfer_efficiency_nolegend_1deg.pdf"), 
       plot=p1 + scale_color_gradientn("", colours = plasma(20), na.value="grey80", guide="none"), width=7, height=4)
ggsave(filename = str_c("plots/",CASE,"/export/jelly_benthic_transfer_efficiency_1deg.pdf"), 
       plot=p1 + scale_color_gradientn("", colours = plasma(20), na.value="grey80"), width=8.6, height=4)


p2 <- ggplot(mapping=aes(x=lon, y=lat)) + 
  geom_point(aes(colour=trans_eff_fctdiff), size=0.35, data=ratios[ratios$taxon=="sum",], shape=15) + labs(x="", y="") + 
  scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) + 
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) + 
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                     panel.grid.minor=element_blank(), panel.grid.major=element_blank())

ggsave(filename = str_c("plots/",CASE,"/export/jelly_v_cobalt_benthic_transfer_efficiency_nolegend_1deg.pdf"), 
       plot=p2 + scale_color_gradientn("", colours = viridis(20), na.value="grey80", guide="none"), width=7, height=4)
ggsave(filename = str_c("plots/",CASE,"/export/jelly_v_cobalt_benthic_transfer_efficiency_1deg.pdf"), 
       plot=p2 + scale_color_gradientn("", colours = viridis(20), na.value="grey80"), width=8.4, height=4)


## Calculate global means, weighted by biome ##
# 1. e-ratio
av_biomesum$percent_of_ocean <- av_biomesum$vol / sum(av_biomesum$vol)
av_biomesum$percent_of_ocean_100 <- av_biomesum$vol100 / sum(av_biomesum$vol100)

biomesum <- ddply(exportres[exportres$taxon=="sum",],~biome,function(x){
  export_flux = mean(x$export * x$depth, na.rm=TRUE) # originally had it as mean(x$export * x$depth100) but e-ratio was half of what I expected
  npp_zint200m = mean(x$npp * x$depth, na.rm=TRUE)
  e_ratio = export_flux / npp_zint200m
  return(data.frame(export_flux,npp_zint200m,e_ratio))
})
biomesum <- cbind(biomesum, percent_of_ocean=av_biomesum$percent_of_ocean)


ratios <- join(ratios, biomes_df)
biomeratio <- ddply(ratios, ~biome, function(x){
  jelly_v_cobalt_100 = mean(x$jelly_v_cobalt_100, na.rm=TRUE)
  jelly_v_cobalt_1000 = mean(x$jelly_v_cobalt_1000, na.rm=TRUE)
  jelly_v_cobalt_seafl = mean(x$jelly_v_cobalt_seafl, na.rm=TRUE)
  j_trans_eff = mean(x$j_trans_eff, na.rm=TRUE)
  trans_eff_fctdiff = mean(x$trans_eff_fctdiff, na.rm=TRUE)
  return(data.frame(jelly_v_cobalt_100,jelly_v_cobalt_1000,jelly_v_cobalt_seafl,j_trans_eff,trans_eff_fctdiff))
})
biomeratio <- cbind(biomeratio, percent_of_ocean=av_biomesum$percent_of_ocean)

print("Biome-weighted GZ-model vs COBALT export at 100 m:")
print(sum(biomeratio$jelly_v_cobalt_100 * biomeratio$percent_of_ocean))
# --> 0.4063026

print("Biome-weighted GZ-model vs COBALT export at 1000 m:")
print(sum(biomeratio$jelly_v_cobalt_1000 * biomeratio$percent_of_ocean))
# --> 1.661407

print("Biome-weighted GZ-model vs COBALT export at the seafloor:")
print(sum(biomeratio$jelly_v_cobalt_seafl * biomeratio$percent_of_ocean))
# --> 2.11825

print("Biome-weighted GZ-model e-ratio:")
print(sum(biomesum$e_ratio * biomesum$percent_of_ocean))
# --> 0.07032388

print("Biome-weighted GZ-model transfer efficiency to depth:")
print(sum(biomeratio$j_trans_eff * biomeratio$percent_of_ocean))
# --> 0.2838655

print("Biome-weighted GZ-model vs COBALT transfer efficiency factor difference:")
print(sum(biomeratio$trans_eff_fctdiff * biomeratio$percent_of_ocean))
# --> 5.108108

# }