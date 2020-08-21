#!/usr/bin/env Rscript
#
#     Analyze baseline results for the 1-degree grid
#     C-flux model, Jessica Luo
#
#----------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("gridExtra")
library("oce")

# import functions and constants
OCEANS_SURFACE_AREA <- 3.619e+14
OCEANS_SURFACE_VOL <- 6.708e+16

source("functions.R")
Spectral <- brewer.pal(8, "Spectral")

CASES = c("0-baseline","1-high_biomass","2-low_biomass")

## { --- Import standard files -------
# import coastline
coastline.world <- read.csv("data/gshhg_world_c.csv")
load("data/coastline.grid_1deg.Rdata") # 1 deg grid

# load biomes
load('data/biomes_chl_mixedlayer.Rdata')
N_BIOME=4

# import area/volume
av <- read.csv("data/surfaceArea_volume_1-deg.csv")
av <- join(av,biomes_df, by=c("lat","lon"))
av$depth.m = av$vol / av$area.m2
av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(vol=sum(x$vol)))})
av_biomesum <- rbind(av_biomesum, data.frame(biome='TOTAL', vol=sum(av_biomesum[av_biomesum$biome != 'TOTAL','vol'])))
av_biomesum$volfrac = av_biomesum$vol/av_biomesum[av_biomesum$biome=="TOTAL","vol"]

av_biomedepth <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(depth.m=mean(x$depth.m)))})
sum(av_biomesum$vol)

# flux variables
vars <- c("ingestion", "respiration", "DOC", "reproduction", "predation", "flux", "egestion")

# } 

for (C_ID in 1:3){
  # C_ID=1
  CASE = CASES[C_ID]
  
  ## { Generate plots for the input biomass ---------------------------------
  
  # 1 degree grid
  rd <- read.csv(str_c("data/",CASE,"/jelly_biomass_1_deg_grid.csv"), as.is=TRUE)
  rd$rank_phylum <- factor(rd$rank_phylum, levels=c("Cnidaria", "Ctenophora", "Chordata"))
  rd$Biomass.gCm3 <- rd$Biomass.mgCm3 / 1000
  
  p1 <- ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(aes(colour=Biomass.gCm3), size=0.3, data=rd, shape=15) + 
    geom_point(colour="grey80", size=0.2, data=rd[which(rd$Biomass.gCm3==0),], shape=15) + 
    labs(x="", y="") +
    scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
    theme_bw() + scale_color_gradientn(expression(bold(paste("SSB \n(g C ", m^bold("-3"), ")"))), 
                                       colours=rev(Spectral), na.value="white", trans="log10") +
    facet_grid(rank_phylum~.) +  theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                                      panel.grid.minor=element_blank(), panel.grid.major=element_blank())
  ggsave(str_c("plots/",CASE,"/baseline/baseline_biomass_1deg.pdf"), plot = p1, width=6.9, height=9)

  
  ## Biomass data; sum by biome
  rd_biome = join(rd, av, by=c("lat", "lon"))

  rd_biome$weights=0
  # generate weights by biome
  gen_weights_by_biome <- function(x, biome=c("COAST", "LC", "HCSS", "HCPS"), taxon=c("Cnidaria", "Ctenophora", "Chordata"), volfrac=1){
    totn=nrow(x[x$biome==biome & x$rank_phylum==taxon,])
    x[x$biome==biome & x$rank_phylum==taxon,"weights"]= (1/totn) * volfrac
    return(x)
  }
  for (b in c("COAST", "LC", "HCSS", "HCPS")){
    for (t in c("Cnidaria", "Ctenophora", "Chordata")){
      rd_biome <- gen_weights_by_biome(rd_biome, biome=b, taxon=t, volfrac=av_biomesum[av_biomesum$biome==b,"volfrac"])
    }
  }
  
  # complete plot
  p2 <- ggplot(rd_biome) + geom_boxplot(aes(x=biome,y=Biomass.gCm3), alpha=0.5) + facet_grid(rank_phylum~.) + 
    scale_y_continuous(trans="log10") + labs(x="") + theme_bw()
  
  require(grid)
  r <- rectGrob(gp=gpar(col="white"))
  
  pdf(str_c("plots/",CASE,"/baseline/baseline_biomass_1deg_with_boxplot.pdf"), width=11, height=9)
  grid.arrange(p1, r, p2, ncol=2, layout_matrix = rbind(c(1,2),c(1,3)), heights=c(0.02,1), widths=c(1,0.5))
  dev.off()
  
  # calculate stats for table
  biomass_biome = ddply(rd_biome, ~rank_phylum+biome, function(x){
    biomass_gm_mean.mgCm3 = gm_mean(x$Biomass.mgCm3)
    biomass_gm_std.mgCm3 = gm_sd(x$Biomass.mgCm3)
    ngrid=nrow(x)
    return(data.frame(biomass_gm_mean.mgCm3, biomass_gm_std.mgCm3, ngrid))
  })

  # total
  set.seed(0)
  biomass_biome_sum = ddply(rd_biome, ~rank_phylum, function(x){
    biomass_gm_mean.mgCm3 = weighted.geomean(x$Biomass.mgCm3, x$weights)
    
    min_nbiome=min(table(x$biome))
    xx = ddply(x, ~biome, function(xi, n=min_nbiome){
      Biomass.mgCm3=sample(xi$Biomass.mgCm3, size = n, replace = FALSE)
      return(data.frame(Biomass.mgCm3))
    })
    
    biomass_gm_std.mgCm3 = gm_sd(xx$Biomass.mgCm3)
    
    ngrid=nrow(x)#sum(x$ngrid)
    biome='TOTAL'
    return(data.frame(biome, biomass_gm_mean.mgCm3, biomass_gm_std.mgCm3, ngrid))
  })
  
  min_nbiome=min(table(rd_biome$biome, rd_biome$rank_phylum))
  global_sd_sample = ddply(rd_biome, ~rank_phylum+biome, function(x, n=min_nbiome){
    Biomass.mgCm3=sample(x$Biomass.mgCm3, size = n, replace = TRUE)
    return(data.frame(Biomass.mgCm3))
  })
  
  global_sd = gm_sd(global_sd_sample$Biomass.mgCm3)
  
  biomass_tot <- data.frame(rank_phylum="ALL", biome="TOTAL", 
                            biomass_gm_mean.mgCm3=sum(biomass_biome_sum$biomass_gm_mean.mgCm3),
                            biomass_gm_std.mgCm3=global_sd, 
                            ngrid=nrow(dcast(rd, lat+lon~rank_phylum, value.var="Biomass.mgCm3")))

  biomass_biome <- rbind(biomass_biome, biomass_biome_sum, biomass_tot)
  biomass_biome <- join(biomass_biome, av_biomesum, by="biome")
  biomass_biome$tot_biomass.gC <- biomass_biome$biomass_gm_mean.mgCm3 / 1000 * biomass_biome$vol
  biomass_biome$sd_biomass.gC <- biomass_biome$biomass_gm_std.mgCm3 / 1000 * biomass_biome$vol
  
  write.csv(biomass_biome, str_c("data/",CASE,"/baseline/biomass_1deg_biomes.csv"), row.names=FALSE)
  # }

  
  ## { Load upper ocean ensemble mean results -------------------------
  # must run analysis-generate_ensemble_means.R
  load(str_c("data/",CASE,"/baseline/upperocean_ens_mean.Rdata"))

  ens_sum <- ens_mean[ens_mean$taxon=="sum",]
  
  ens_sum$allfood <- ens_sum$ingestion + ens_sum$egestion
  ens_sum$allpocproduction <- ifelse(ens_sum$flux >= 0, yes=ens_sum$flux, no=0) + ens_sum$egestion
  # }
  
  ## { Plotting global maps ------------------------------------------------
  
  # assign a non-zero value to all food for plotting
  #ens_sum[which(ens_sum$allfood ==0),"allfood"] <- 1e-7
  plotd <- ens_sum
  plotd[,vars] <- ens_sum[,vars]/ens_sum$allfood * 100 #convert to percentage
  plotd <- join(plotd, av[,c("lat","lon","depth.m")])
  plotd$allfood_gCm2 = plotd$allfood * plotd$depth.m
  #plotd <- melt(plotd, id.vars = c("lat","lon","taxon"))
  
  ingestp <- ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(aes(colour=allfood_gCm2), size=0.3, data=plotd, shape=15) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) + labs(x="", y="") +
    scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_color_gradientn(expression(paste("C_I\n(g C ", m^"-2", " ", y^"-1", ")")), colours=rev(Spectral), na.value="grey75", trans="log10")
  
  egestp <- ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(aes(colour=egestion), size=0.3, shape=15, data=plotd) +
    labs(x="", y="") + scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
    theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_color_gradientn("% of Input", colours=brewer.pal(8, "Purples"), limits=c(0, 50), na.value="grey80")
  
  respp <- ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(colour="black", size=0.2, shape=15, data=plotd[plotd$respiration > 100,]) +
    geom_point(aes(colour=respiration), size=0.3, shape=15, data=plotd) +
    labs(x="", y="") + scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
    theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
    scale_color_gradientn("% of Input", colours=brewer.pal(8, "Reds"), limits=c(0, 100), na.value="grey80")
  
  fluxp <- ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(aes(colour=flux), size=0.3, shape=15, data=plotd) +
    geom_point(colour='black', size=0.3, shape=15, data=plotd[which(plotd$flux < 0),]) +
    labs(x="", y="") + scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
    theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_color_gradientn("% of Input", colours=brewer.pal(8, "Blues"), limits=c(0,55), na.value="grey80")
  
  predp <- ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(aes(colour=predation), size=0.3, shape=15, data=plotd) +
    labs(x="", y="") + scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
    theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_color_gradientn("% of Input", colours=brewer.pal(8, "Greens"), limits=c(0, 60), na.value="grey80")
  
  ggsave(plot=ingestp, file=str_c("plots/",CASE,"/baseline/1deg_1_base_fluxes_ingestion_SUM.pdf"), width=5.7, height=2.5)
  ggsave(plot=egestp, file=str_c("plots/",CASE,"/baseline/1deg_2_base_fluxes_egestion_pct_SUM.pdf"), width=5.4, height=2.5)
  ggsave(plot=respp, file=str_c("plots/",CASE,"/baseline/1deg_3_base_fluxes_respiration_pct_SUM.pdf"), width=5.4, height=2.5)
  ggsave(plot=fluxp, file=str_c("plots/",CASE,"/baseline/1deg_4_base_fluxes_carcasses_pct_SUM.pdf"), width=5.4, height=2.5)
  ggsave(plot=predp, file=str_c("plots/",CASE,"/baseline/1deg_5_base_fluxes_predation_pct_SUM.pdf"), width=5.4, height=2.5)
  
  # }
  
  
  # Plot production to biomass ratio
  d = ens_mean[which(ens_mean$taxon != "sum"),]
  d$production.gCm3 = d$ingestion - d$respiration
  rd$taxon=tolower(str_sub(rd$rank_phylum, 1, 4)) 
  d <- join(d,rd[,c("taxon","lat","lon","Biomass.gCm3")], by=c("taxon","lat","lon"))
  d$p_b = d$production.gCm3/d$Biomass.gCm3
  
  # get rid of NA values and extraneous point
  d <- d[!is.na(d$p_b),]
  #d <- d[which(d$p_b < 3e3),]
  
  d$taxon <- factor(d$taxon, levels=c("cnid","cten","chor"))

  gm_mean_pb = ddply(d, ~taxon, function(x){data.frame(p_b=gm_mean(x$p_b/365*100))})
  
  # calculate a biome-scaled mean P/B & range
  range_pb = ddply(d, ~taxon, function(x){t(range(x$p_b))})
  names(range_pb)=c("taxon","lower_range","upper_range")
  
  # res=data.frame()
  # for (t in c("cnid","cten","chor")){
  #   for (b in c("COAST", "LC", "HCSS", "HCPS")){
  #     x=d[d$taxon==t & d$biome==b,"p_b"]
  #     x=x/365*100
  #     bci=return_bci(x)
  #     tmp=data.frame(taxon=t,biome=b,lower=bci[1],upper=bci[2], stringsAsFactors=FALSE)
  #     res=rbind(res,tmp)
  #   }
  # }
  # res

  
  mean_pb = ddply(d, ~taxon+biome, function(x){data.frame(p_b=mean(x$p_b/90*100))})
  mean_pb = join(mean_pb, av_biomesum)
  mean_pb$p_b_scaled = mean_pb$p_b * mean_pb$volfrac
  mean_pb_biomescaled = join(ddply(mean_pb, ~taxon, function(x){data.frame(p_b=sum(x$p_b_scaled))}),
                             range_pb)
  print(mean_pb_biomescaled)
  
  
  # plotting
  TaxaLabels <- c("Cnidarians", "Ctenophores", "Tunicates")
  names(TaxaLabels) <- c("cnid","cten","chor")
  
  # ggplot(d, aes(x=Biomass.gCm3, y=production.gCm3)) + 
  #   geom_point(alpha=0.5) + geom_smooth(method="lm") + facet_wrap(~taxon, scales="free") + theme_bw() +
  #   scale_x_continuous("Biomass (g C m-3)", trans="log10") + scale_y_continuous("Production (g C m-3)", trans="log10")

  p1 = ggplot() + 
    geom_point(aes(x=Biomass.gCm3, y=p_b), alpha=0.5, size=1, data=d) + 
    geom_hline(aes(yintercept=p_b), color="red", linetype=2, data=gm_mean_pb) + 
    facet_grid(taxon~., labeller=labeller(taxon=TaxaLabels)) + theme_bw() +
    scale_x_continuous(expression(paste("Biomass (g C ", m^"-3", ")")), trans="log10") + 
    scale_y_continuous("Production/Biomass (P/B)", trans="log10")
  
  p2 = ggplot(mapping=aes(x=lon, y=lat)) +
    geom_point(aes(colour=p_b), size=0.11, shape=15, data=d) +
    labs(x="", y="") + scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
    scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
    geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) +
    theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    facet_grid(taxon~., labeller=labeller(taxon=TaxaLabels)) +  theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), 
                                       panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + 
    scale_color_gradientn("P/B", colours=rev(Spectral), trans="log10", na.value="grey80")
  #ggsave(str_c("plots/",CASE,"/baseline/production_to_biomass.pdf"), plot = p2, width=6.9, height=9)
  
  require(grid)
  r <- rectGrob(gp=gpar(col="white"))
  
  pdf(str_c("plots/",CASE,"/baseline/production_to_biomass.pdf"), width=10, height=7.5)
  grid.arrange(r, p1, p2, ncol=2, layout_matrix = rbind(c(1,3),c(2,3)), heights=c(0.03,1), widths=c(0.5,1))
  dev.off()
  
  
  # lm(log10(production.gCm3)~log10(Biomass.gCm3), data=d[d$Biomass.gCm3 != 0,])
  # lm(log10(production.gCm3)~log10(Biomass.gCm3), data=d[d$taxon=="cnid",])
  # lm(log10(production.gCm3)~log10(Biomass.gCm3), data=d[d$taxon=="cten",])
  # lm(log10(production.gCm3)~log10(Biomass.gCm3), data=d[d$taxon=="chor" & d$Biomass.gCm3 != 0,])
  
  
}


