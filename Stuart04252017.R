
setwd("C:/Users/stuart/Google Drive/DataForAnalysis/grasshoppers/SexCombined")

data <- read.csv("C1_1959-2015_eggDiapause.csv")
data <- read.csv("B1_1959-2015_eggDiapause.csv")
data <- read.csv("CHA_1959-2012_eggDiapause.csv")

# Check for missing data
sum(is.na(data$ordinal))
sum(is.na(data$in6))
sum(is.na(data$year))

# Check species names
data$species <- gsub(" $", "", data$species) # Remove any spaces accidentally
                                             # added at end of species names
sort(unique(data$species))

# Subset to common species (B1 and CHA only)
# B1
comsps <- c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus boulderensis",
            "Melanoplus dawsoni", "Melanoplus packardii")
data <- data[data$species %in% comsps,]
# CHA
comsps <- c("Aeropedellus clavatus", "Hesperotettix viridis", "Melanoplus bivittatus",
            "Melanoplus dawsoni", "Melanoplus sanguinipes")
data <- data[data$species %in% comsps,]

# Check years
unique(data$year)
# CHA has no 2006 and no 2013-2015

# Subset to adults
data <- data[data$in6 > 0,]

#==============================================
# Calculating total overlap by species and year
#==============================================

# Pianka's overlap metric is calculated using the proportions of total individuals
# of each species observed in the entire season that were present at each sampling
# period. The proportions of species j and k present in sampling period i are
# written as pij and pik respectively. Overlap is calculated from these proportions
# as follows:
# overlap = sum(pij * pik) / sqrt(sum(pij^2) * sum(pik^2))
# Need to create vectors of pij and pik before calculating

# The code below calculates total overlap experienced by each species in each
# year at site C1

# Get list of species and years
years <- unique(data$year)
sps <- sort(unique(data$species))

# Set up data storage
overlap <- matrix(NA, nrow = length(years), ncol = length(sps))
colnames(overlap) <- sps
rownames(overlap) <- years

# Loop through years
for(j in 1:length(years)) {
  # Subset to year
  dat <- data[data$year == years[j],]
  
  # Create matrix of abundances for all species at each sampling event
  sumByOrd <- tapply(dat$in6, list(dat$species, dat$ordinal), FUN = sum)
  # Replace NAs with 0
  sumByOrd[is.na(sumByOrd)] <- 0
  # Calculate total individuals of each species observed in entire season
  totals <- apply(sumByOrd, MARGIN = 1, FUN = sum)
  
  # Loop through species calculating overlap for each one
  for(i in 1:nrow(sumByOrd)) {
    # Calculate proportion of focals observed in entire season present at each
    # sampling event (if focal species was detected in this year)
    if(totals[i] > 0) {
      pij <- as.vector(sumByOrd[i,] / totals[i])
      # Calculate same proportions for other species combined
      sumOthers <- apply(sumByOrd[-i,], MARGIN = 2, FUN = sum)
      pik <- sumOthers / sum(sumOthers)
      # Calculate overlap and store result
      overlap[j,i] <- sum(pij * pik) / sqrt(sum(pij^2) * sum(pik^2))
      
    } # Close if statement
  } # Close species loop
} # Close years loop

# Look at variation in overlap across years
plot(overlap[,5])

#======================================
# Calculating measures of seasonal temp
#======================================

# There are a number of different ways to summarize temperature in each year.
# We will use the number of GDDs accumulated because this has a mechanistic link
# to development rate of these species. We could look at the total number of
# GDDs accumulated in the season, but it might also be that the variation in
# GDD accumulation over the season is important. To look at the timing of warming
# we can look at the number of GDDs accumulated midway through the season.

setwd("C:/Users/stuart/Google Drive/DataForAnalysis/climate")
allClim <- read.csv("AlexanderClimateAll.csv")

# Subset to site, relevant years, and dates relevant to grasshopper season
# (March 1 to Aug 31 = ordinal 60 to 243)
clim <- allClim[allClim$Site == "NOAA" & allClim$Year %in% years &
                  allClim$Julian > 59 & allClim$Julian < 244,] # Remember to change site!

# Check for missing temperature data and remove
clim[is.na(clim$Min) | is.na(clim$Max),c("Year", "Julian")]
# C1 - 8 days missing in 2013 - probably no need to worry about this small number, no 
# climate data from 2015
# B1 - only 1 day missing in 2012, so nothing to worry about, no climate data for 2013
# CHA - no missing days
clim <- clim[!is.na(clim$Min) & !is.na(clim$Max),]

# Calculate GDDs accumulated on each day
setwd("C:/Users/stuart/Google Drive/University of Washington/Buckley lab/Phenological overlap project/HopperPhenology")
source("degreedays.R")
clim$dd <- apply(clim[,c("Min","Max")], MARGIN = 1, FUN = degree.days.mat, 
                 LDT= 12)

# Calculate cumulative number of GDDs accumulated by each day in each year
# (ignoring GDDs accumulated before March 1st)
library(dplyr)
clim <- clim %>% arrange(Year, Julian) %>% group_by(Year) %>% 
  mutate(cdd = cumsum(dd))

# Based on Fig 2 from Nufio et al (2010 PLoS one) it seems like the GDD 
# accumulation by ordinal dates 160, 200, and 240 might reveal relevant 
# differences between years.
ords <- c(160, 200, 240)
tempMeas <- clim[clim$Julian %in% ords, c("Year", "Julian", "cdd")]
plot(x = tempMeas$Year, y = tempMeas$cdd, col = as.factor(tempMeas$Julian))

# Now calculate, for each year, number of GDDs accumulated between ordinal
# dates 141 and 170, and between ordinal dates 171 and 200 (inclusive)
dd141_170 <- clim[clim$Julian > 140 & clim$Julian <= 170,] %>% 
  group_by(Year) %>% summarize(dds = sum(dd))

dd171_200 <- clim[clim$Julian > 170 & clim$Julian <= 200,] %>%
  group_by(Year) %>% summarize(dds = sum(dd))

# Create a data frame of the different temperature metrics
Temp <- data.frame(dd141_170$Year,
                   tempMeas[tempMeas$Julian == 160, "cdd"],
                   tempMeas[tempMeas$Julian == 200, "cdd"],
                   tempMeas[tempMeas$Julian == 240, "cdd"],
                   dd141_170$dds,
                   dd171_200$dds)
names(Temp) <- c("Year", "GDDa160", "GDDa200", "GDDa240", "GDD141_170",
                 "GDD171_200")

#================================================
# Relating overlap to various temperature metrics
#================================================

# Combine overlap and temp data
# C1
overlapTemp <- cbind(Temp, overlap[1:nrow(overlap)-1,])

# B1
overlapTemp <- cbind(Temp, overlap)

# CHA
overlapTemp <- cbind(Temp, overlap)

# Add whether year was initial or resurvey
overlapTemp$period <- as.factor(c("Initial", "Initial", 
                                  rep("Resurvey", times = nrow(overlapTemp) - 2)))

par(mar = c(0,0,0,0), oma = c(4,4,1,3))
par(mfcol = c(5,5))

plot(y = overlapTemp[,7], x = overlapTemp[,2], xaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,8], x = overlapTemp[,2], xaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,9], x = overlapTemp[,2], xaxt = "n", col = overlapTemp$period)
mtext("Overlap", side = 2, line = 2, cex = 1)
plot(y = overlapTemp[,10], x = overlapTemp[,2], xaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,11], x = overlapTemp[,2], col = overlapTemp$period)
mtext("Up to ord 160", side = 1, line = 2, cex = 0.6)

plot(y = overlapTemp[,7], x = overlapTemp[,3], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,8], x = overlapTemp[,3], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,9], x = overlapTemp[,3], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,10], x = overlapTemp[,3], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,11], x = overlapTemp[,3], yaxt = "n", col = overlapTemp$period)
mtext("Up to ord 200", side = 1, line = 2, cex = 0.6)

plot(y = overlapTemp[,7], x = overlapTemp[,4], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,8], x = overlapTemp[,4], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,9], x = overlapTemp[,4], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,10], x = overlapTemp[,4], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,11], x = overlapTemp[,4], yaxt = "n", col = overlapTemp$period)
mtext("Up to ord 240", side = 1, line = 2, cex = 0.6)

plot(y = overlapTemp[,7], x = overlapTemp[,5], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,8], x = overlapTemp[,5], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,9], x = overlapTemp[,5], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,10], x = overlapTemp[,5], xaxt = "n", yaxt = "n", col = overlapTemp$period)
plot(y = overlapTemp[,11], x = overlapTemp[,5], yaxt = "n", col = overlapTemp$period)
mtext("Ords 141-170", side = 1, line = 2, cex = 0.6)

plot(y = overlapTemp[,7], x = overlapTemp[,6], xaxt = "n", yaxt = "n", col = overlapTemp$period)
#mtext("M.bouderensis", side = 4, line = 0, cex = 0.5)
plot(y = overlapTemp[,8], x = overlapTemp[,6], xaxt = "n", yaxt = "n", col = overlapTemp$period)
#mtext("M.sanguinipes", side = 4, line = 0, cex = 0.5)
plot(y = overlapTemp[,9], x = overlapTemp[,6], xaxt = "n", yaxt = "n", col = overlapTemp$period)
#mtext("M.fasciatus", side = 4, line = 0, cex = 0.5)
plot(y = overlapTemp[,10], x = overlapTemp[,6], xaxt = "n", yaxt = "n", col = overlapTemp$period)
#mtext("C.pellucida", side = 4, line = 0, cex = 0.5)
plot(y = overlapTemp[,11], x = overlapTemp[,6], yaxt = "n", col = overlapTemp$period)
#mtext("C.abdominalis", side = 4, line = 0, cex = 0.5)
mtext("Ords 171-200", side = 1, line = 2, cex = 0.6)
