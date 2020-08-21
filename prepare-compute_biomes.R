#
#
#   Compute Biomes
#     Use 3-biome formulation as in Stock et al. 2014
#     In 2-degree regular grid
#
#
#----------------------------------------------------


library("plyr")
library("reshape2")
library("ncdf4")
library("matlab")
library("oce")
library("ggplot2")
library("stringr")

## ------------------ FUNCTIONS ------------------------
source(file = "functions.R")

seafloor.depth <- read.csv("data/woa_seafloor_depth.csv")
seafloor.depth$loc <- ifelse(seafloor.depth$depth <= 200, "coast", "ocean") 

# load data
# Mixed Layer Depth data from de Boyer Montegut et al. 2004
# MLD defined in density terms with a fixed threshold criterion (0.03 kg/m3)
nc <- nc_open("raw_data/mld_DR003_c1m_reg1.0.nc")
#nc <- nc_open("raw_data/Argo_mixedlayers_monthlyclim_03192017.nc")
mld <- ncvar_get(nc, "mld") # dims : [lon,lat,time]
mld[mld>1e+07] <- NA #mld[mld==1e9] <- NA

lat_1x <- seq(-89.5,89.5,by=1)
lon_1x <- seq(-179.5,179.5,by=1)

nc <- nc_open("raw_data/S19972472010334.L3m_CU_CHL_chl_ocx_regularx1.nc")
chl <- ncvar_get(nc, "chl_ocx") # dims : [lon,lat]

max_mld <- apply(mld, c(1,2), max, na.omit=TRUE) # dims : [lon,lat]

biomes = array(NA, dim = c(length(lon_1x),length(lat_1x)))

# low chlorophyll
biomes[which(chl <= 0.125)] = 2#"LC"

# high chlorophyll, seasonally stratified
biomes[which(chl > 0.125 & max_mld > 75)] = 3#"HCSS"

# high chlorophyll, permanently stratified
biomes[which(chl > 0.125 & max_mld <= 75)] = 4#"HCPS"


biomes_df <- as.data.frame(biomes)
names(biomes_df) <- str_c("X",str_replace(as.character(lat_1x),'-','_'))
biomes_df <- data.frame(lon=lon_1x, biomes_df)

biomes_df <- melt(biomes_df, id.vars = "lon", variable.name = "lat", value.name = "biome_num")
biomes_df$lat <- as.character(biomes_df$lat)
biomes_df$lat <- str_replace(biomes_df$lat, "X", "")
biomes_df$lat <- as.numeric(str_replace(biomes_df$lat, "_", "-"))

# coastal biome, defined as max depth <= 200 m
biomes_df <- join(biomes_df,seafloor.depth)
biomes_df[which(biomes_df$loc=="coast"),"biome_num"] <- 1 #COAST
#  1   2   3
# 23 547 777

maxmld_df <- as.data.frame(max_mld)
names(maxmld_df) <- str_c("X",str_replace(as.character(lat_1x),'-','_'))
maxmld_df <- data.frame(lon=lon_1x, maxmld_df)
maxmld_df <- melt(maxmld_df, id.vars = "lon", variable.name = "lat", value.name = "max_mld")
maxmld_df$lat <- biomes_df$lat

# plotting('max_mld', maxmld_df)


# add missing values from the dataset, mostly coastal values
rd = read.csv("data/0-baseline/jelly_biomass_1_deg_grid.csv", as.is=TRUE)
biomes_df$id = row.names(biomes_df)
temp <- join(rd, biomes_df)
missingbiome <- temp[which(is.na(temp$biome)),]
plotting("Biomass.mgCm3", missingbiome) # mostly coastal
biomes_df[which(biomes_df$id %in% missingbiome$id),"biome_num"] <- 1

# add names to biome numbers
biomes_df <- join(biomes_df, data.frame(biome_num=1:4,biome=c("COAST","LC","HCSS","HCPS")))
biomes_df$biome <- factor(biomes_df$biome, levels=c("COAST","LC","HCSS","HCPS"))
biomes_df <- biomes_df[,c("lon","lat","biome")]

save(biomes_df, file = "data/biomes_chl_mixedlayer.Rdata")

coastline.world <- read.csv("data/gshhg_world_c.csv")

p <- ggplot(mapping=aes(x=lon, y=lat)) +
  geom_point(aes(colour=biome), size=0.5, data=biomes_df, shape=15) + labs(x="", y="") +
  scale_x_continuous(breaks=c(-90, 0, 90), expand=c(0,0)) +
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60), expand=c(0,0)) +
  geom_polygon(data=coastline.world, fill="grey60", size=0.1) + scale_shape_discrete(solid=T) + scale_color_discrete("Biome",na.value="grey40") + 
  theme_bw() + theme(plot.margin=unit(c(1,0.75,0.25,0.75), "lines"), panel.grid.minor=element_blank(), panel.grid.major=element_blank())

ggsave(filename = "plots/biomes_1degree.pdf",p,height=4.7,width=9.8)

