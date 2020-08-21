#
#
#     Calculate Salp respiration allometric relationship
#
###########################################################

library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)

## { Respiration --------------------------------------------------------------
t <- read.csv("raw_data/salps_respiration.csv", stringsAsFactors=FALSE)
t <- t[t$order=="Salpida",]

# convert all to C weight in mg
# C = a L^b

# bootstrapping  

tb <- ddply(t, ~genus_species+generation, function(x){
  if(!is.na(x$C_weight_mg[1])){
    C_weight_mg <- rnorm(n=x$n, mean = x$C_weight_mg, sd=x$C_weight_error)
    x <- x[rep(seq_len(nrow(x)), x$n), ]
    x$C_weight_mg <- C_weight_mg
  }
  if(is.na(x$C_weight_mg[1])){
    low_range <- as.numeric(str_sub(x$length_range_mm[1], 1, 2))
    high_range <- as.numeric(str_sub(x$length_range_mm[1], 4, -1))
    length_mm <- runif(sum(x$n), low_range, high_range)
    C_weight_mg <- x$a_length_to_weight[1] * length_mm ^ x$b_length_to_weight[1]
    
    x <- x[rep(seq_len(nrow(x)), x$n), ]
    x$length_mm <- length_mm
    x$C_weight_mg <- C_weight_mg
  } 

  x$rate <- x$a_weight_to_rate * x$C_weight_mg ^ x$b_weight_to_rate
  
  return(x)
})

tb <- tb[,c("phylum", "order", "genus_species", "generation","temp_mean", "length_mm", "C_weight_mg", "rate", "rate_units", "biomass_units", "n")]
tb <- tb[tb$C_weight_mg > 0,]

# 44.66 micromol O2 = 1 mL O2
# molar volume at STP (0-deg C) for O2 = 22.4 L/mol
tb$rate <- tb$rate / 44.66
tb$rate_units <- "mL O2 ind-1 h-1"

tb$temp_correction <- exp(-7542.938 / (273 + tb$temp_mean)) 
tb$rate_temp_corr <- tb$rate/tb$temp_correction
tb$log_resp <- log10(tb$rate_temp_corr)
tb$log_biomass <- log10(tb$C_weight_mg)

plot(log_resp ~ log_biomass, data=tb)

#write.csv(tb, "data/salps_respiration_bootstrapped_calculated.csv", row.names=F)

m <- lm(log_resp~log_biomass, data=tb)
summary(m)
capture.output(summary(m), file="plots/allometric_relationships/salps_respiration_regression.txt")


label <- data.frame(x=c(1,2), value=c(paste("y =", round(m$coefficients[1], 2), "+", round(m$coefficients[2],2), "* x"), paste("r^2 =", round(summary(m)$adj.r.squared,3))))
p <- ggplot(data=tb, aes(x=log_biomass, y=log_resp)) + geom_point(alpha=0.8) + geom_smooth(method="lm", linetype=5, size=1.5, se=FALSE) + labs(x="Log B (mg C)", y=expression(paste("Log R (mL ", O[2]," ", ind^-1, " ", h^-1, " / ", exp^(-Ea/kT), ")", sep=" "))) + theme_bw() #+ geom_text(aes(x=-0.9, y=10.5-label$x/4, label=label$value, hjust=0)); p
ggsave(p, file="plots/Salps_Resp_regression.pdf", height=5, width=6.5)

d <- read.csv("raw_data/salps_weight-specific-respiration.csv", stringsAsFactors = FALSE)
q10 = 2.3
d$temp_correction = q10^(-d$temp_mean/10)
d$weight_specific_umolO2.mgChr_0deg = d$Weight.specific * d$temp_correction
d$weight_specific_mgC.mgCd_0deg = d$weight_specific_umolO2.mgChr_0deg * 1/(44.62) * 24 * 0.5359 * 0.8 #RQ can also be 0.95
mean(d$weight_specific_mgC.mgCd_0deg)
# --> 0.03041577

# }


## { Clearance Rate -----------------------------------------------------
d <- read.csv("raw_data/salps_clearancerate.csv", stringsAsFactors=FALSE)
d <- d[d$length_range_mm != "",]

db <- ddply(d, ~genus_species+generation, function(x){
  length_range <- str_split_fixed(x$length_range_mm, "-", 2)
  C_weight_mg <- NULL
  length <- NULL
  
  for (i in 1:nrow(x)){
    # bootstrapping
    low_range <- as.numeric(length_range[i,1])
    high_range <- as.numeric(length_range[i,2])
    length_mm <- runif(x$n[i], low_range, high_range)
    
    # calculating a carbon weight from the length
    C_weight_mg_temp <- x$a_length_to_weight[i] * length_mm ^ x$b_length_to_weight[i]
    
    C_weight_mg <- c(C_weight_mg, C_weight_mg_temp)
    length <- c(length, length_mm)
  }  

  x <- x[rep(seq_len(nrow(x)), x$n), ]
  x$length_mm <- length
  x$C_weight_mg <- C_weight_mg 
  
  x$rate <- x$a_weight_to_rate * x$length_mm ^ x$b_weight_to_rate
  
  return(x)
})

db <- db[,c("phylum", "order", "genus_species", "generation","temp_mean", "length_mm", "C_weight_mg", "rate", "rate_units", "biomass_units", "n")]

db$temp_correction <- exp(-7542.938 / (273 + db$temp_mean))
db$rate_temp_corr <- db$rate/db$temp_correction
db$rate_temp_corr <- db$rate_temp_corr / 1000 * 24 # convert from mL hr-1 to L d-1
db$log_cr <- log10(db$rate_temp_corr)
db$log_biomass <- log10(db$C_weight_mg)

plot(log_cr ~ log_biomass, data=db)

write.csv(db, "data/salps_CR_bootstrapped_calculated.csv", row.names=F)

m <- lm(log_cr ~ log_biomass, data=db)
summary(m)
capture.output(summary(m), file="plots/allometric_relationships/salps_CR_regression.txt")

label <- data.frame(x=c(1,2), value=c(paste("y =", round(m$coefficients[1], 2), "+", round(m$coefficients[2],2), "* x"), paste("r^2 =", round(summary(m)$adj.r.squared,3))))
p <- ggplot(data=db, aes(x=log_biomass, y=log_cr)) + geom_point(alpha=0.8) + geom_smooth(method="lm", linetype=5, size=1.5, se=FALSE) + labs(x="Log B (mg C)", y=expression(paste("Log CR (L ", d^-1, " / ", exp^(-Ea/kT), ")"))) + theme_bw()
# + geom_text(aes(x=-1.9, y=14.5-label$x/2, label=label$value, hjust=0))
ggsave(p, file="plots/Salps_CR_regression.pdf", height=5, width=6.5)


# }