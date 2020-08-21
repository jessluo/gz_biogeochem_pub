#!/usr/bin/env Rscript
#
#     Combine biome results
#
#
# -------------------------------------------

# load libraries ----------------------------
library("plyr")
library("reshape2")
library("ggplot2")
library("stringr")

CASES = c("0-baseline","1-high_biomass","2-low_biomass")

#### Pull together values with quantiles for Fig 2. ####
## Upper Ocean ##
ds = data.frame()

for (C_ID in 1:3){
  # C_ID = 1
  
  CASE = CASES[C_ID]
  upper_ocean = read.csv(str_c("data/",CASE,"/mc_res/Biomes_UO_sum_MeanResults.csv"), as.is=TRUE)

  upper_ocean <- rbind(data.frame(dcast(upper_ocean, taxon~var, value.var='lower_tot', fun=sum), value='lower'),
                       data.frame(dcast(upper_ocean, taxon~var, value.var='mean_tot', fun=sum), value='mean'),
                       data.frame(dcast(upper_ocean, taxon~var, value.var='upper_tot', fun=sum), value='higher'))

  upper_ocean$case=CASE

  if(C_ID==1){
    ds = upper_ocean
  } else {
    ds = rbind(ds, upper_ocean)
  }
}  

ds$totIngestion = ds$ingestion + ds$egestion
ds$totDOC = ds$reproduction + ds$DOC

summed <- ddply(ds, ~case+value, function(x){
  colSums(x[,which(sapply(x, class) == "numeric")])
})

ds <- rbind(ds, data.frame(summed, taxon="sum"))
ds$taxon=factor(ds$taxon, levels=c("cnid","cten","chor","sum"))
ds = arrange(ds,case,value,taxon)
ds <- ds[,c("case","value","taxon","totIngestion","ingestion","egestion","respiration","totDOC","DOC","reproduction","predation","flux")]
write.csv(ds, "data/uo_results_CIs_cases_combined.csv", row.names=FALSE)

## Export to Depth ##
dd = data.frame()
for (C_ID in 1:3){
  #C_ID = 1
  
  CASE = CASES[C_ID]
  three_d = read.csv(str_c("data/",CASE,"/mc_res/Biomes_Depth_Combined_Fluxes_CIs.csv"), as.is=TRUE)
  
  three_d <- ddply(three_d, ~depth+taxon+speed+value, function(x){
    return(data.frame(flux_tot=sum(x$flux_tot)))
    })
  
  three_d$case=CASE
  
  if(C_ID==1){
    dd = three_d
  } else {
    dd = rbind(dd, three_d)
  }
}  

summed <- ddply(dd, ~depth+speed+case+value, function(x){
  return(data.frame(flux_tot=sum(x$flux_tot)))
})
dd <- rbind(dd, data.frame(summed, taxon="sum"))

dd$taxon=factor(dd$taxon, levels=c("cnid","cten","chor","sum"))
dd$speed=factor(dd$speed, levels=c("all","slow","fast"))
dd$value=factor(dd$value, levels=c("mean","low","high"))
dd$depth=factor(dd$depth, levels=c("Export","Sequestration","Seafloor"))

dd = arrange(dd, case, depth, taxon, speed, value)
write.csv(dd, "data/depth_results_CIs_cases_combined.csv", row.names=FALSE)


#### Global Biome Means ####
# pull together all means for table S6.
db = data.frame()

for (C_ID in 1:3){
  # C_ID = 1
  
  CASE = CASES[C_ID]
  upper_ocean = read.csv(str_c("data/",CASE,"/mc_res/Biomes_UO_sum_MeanResults.csv"), as.is=TRUE)
  three_d = read.csv(str_c("data/",CASE,"/mc_res/Biomes_Depth_Combined_Fluxes_CIs.csv"), as.is=TRUE)
  
  upper_ocean <- rbind(data.frame(dcast(upper_ocean, biome+taxon~var, value.var='mean_tot'), CI='mean'),
                       data.frame(dcast(upper_ocean, biome+taxon~var, value.var='lower_tot'), CI='lower'),
                       data.frame(dcast(upper_ocean, biome+taxon~var, value.var='upper_tot'), CI='upper'))
  
  colvars=c("depth","biome","taxon","flux_tot")
  three_d_m <- dcast(three_d[three_d$value=="mean", colvars], biome+taxon~depth, value.var = "flux_tot")
  three_d_l <- dcast(three_d[three_d$value=="low" & three_d$speed=="slow", colvars], biome+taxon~depth, value.var = "flux_tot")
  three_d_h <- dcast(three_d[three_d$value=="high" & three_d$speed=="fast", colvars], biome+taxon~depth, value.var = "flux_tot")
  
  three_d <- rbind(data.frame(three_d_m, CI='mean'), data.frame(three_d_l, CI='lower'), data.frame(three_d_h, CI='upper'))
    
  joined = join(upper_ocean, three_d, by=c("biome","taxon","CI"))
  joined$case=CASE

  
  if(C_ID==1){
    db = joined
  } else {
    db = rbind(db, joined)
  }
}  

db$totIngestion = db$ingestion + db$egestion

taxa_global = ddply(db, ~case+CI+taxon, function(x){
  vals=colSums(x[,which(sapply(x, class) == "numeric")])
  return(data.frame(biome="SUM",t(vals)))
})

biome_sum = ddply(db, ~case+CI+biome, function(x){
  vals=colSums(x[,which(sapply(x, class) == "numeric")])
  return(data.frame(taxon="all",t(vals)))
})

global = ddply(db, ~case+CI, function(x){
  vals=colSums(x[,which(sapply(x, class) == "numeric")])
  return(data.frame(biome="SUM",taxon="all",t(vals)))
})

db <- rbind(db, taxa_global, biome_sum, global)
db$taxon <- factor(db$taxon, levels=c("cnid","cten","chor","all"))
db$biome <- factor(db$biome, levels=c("COAST","HCPS","HCSS","LC","SUM"))

db <- rbind(arrange(db[db$CI=="mean",], case, CI, taxon, biome),
            arrange(db[db$CI!="mean",], case, taxon, biome, CI))
db <- arrange(db, case)

db <- db[,c("case", "CI", "taxon", "biome","totIngestion", "egestion",
            "respiration","reproduction","DOC","predation","flux",
            "Export", "Seafloor", "Sequestration")]

write.csv(db, "data/all_results_CIs_cases_combined.csv", row.names=FALSE)

# calculate non-coastal transfer efficiency
# import area/volume
av <- read.csv("data/surfaceArea_volume_1-deg.csv")
av <- join(av,biomes_df, by=c("lat","lon"))
av$depth.m = av$vol / av$area.m2 #vol in m3
av_biomesum <- ddply(av[!is.na(av$biome),],~biome,function(x){return(data.frame(vol=sum(x$vol)))})
names(av_biomesum) <- c("biome", "vol")
av_biomesum <- rbind(av_biomesum, data.frame(biome='TOTAL', vol=sum(av_biomesum[av_biomesum$biome != 'TOTAL','vol'])))
av_biomesum$volfrac = av_biomesum$vol/av_biomesum[av_biomesum$biome=="TOTAL","vol"]

# non-coastal
nc = db[db$biome %ni% c("COAST","SUM") & db$taxon=="all",]
# transfer efficiency
nc <- ddply(nc, ~case+biome+CI, function(x){
  export_to_seq=x$Sequestration/x$Export
  export_to_seafl=x$Seafloor/x$Export
  return(data.frame(export_to_seq,export_to_seafl))
})

nc <- join(nc,av_biomesum[,c("biome","volfrac")], by="biome")
nc <- ddply(nc, ~case+CI, function(x){
  export_to_seq = sum(x$export_to_seq * x$volfrac)
  export_to_seafl = sum(x$export_to_seafl * x$volfrac)
  return(data.frame(export_to_seq,export_to_seafl))
})
####
