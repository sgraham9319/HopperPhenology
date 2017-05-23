#######################################
# Calculating pairwise overlap and GDDs
#######################################

# Load packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(zoo)

# Source degree days file
setwd("C:/Users/stuart/Google Drive/University of Washington/Buckley lab/Phenological overlap project/HopperPhenology")
source("degreedays.R")

# Load climate data
setwd("C:/Users/stuart/Google Drive/DataForAnalysis/climate")
allClim <- read.csv("AlexanderClimateAll.csv")

# Load grasshopper data
setwd("C:/Users/stuart/Google Drive/DataForAnalysis/grasshoppers/SexCombined")
C1 <- read.csv("C1_1959-2015_eggDiapause.csv")
B1 <- read.csv("B1_1959-2015_eggDiapause.csv")
CHA <- read.csv("CHA_1959-2012_eggDiapause.csv")

#========================
# Format grasshopper data
#========================

# Add site column and combine sites
C1$site <- "C1"
B1$site <- "B1"
CHA$site <- "CHA"
data <- rbind(CHA, B1, C1)

# Check for missing data
sum(is.na(data$ordinal))
sum(is.na(data$in6))
sum(is.na(data$year))

# Fix typos in species names
data$species <- gsub(" $", "", data$species)
data$species <- gsub("Melanoplus bouderensis", "Melanoplus boulderensis", data$species)
sort(unique(data$species))

# Subset to common species
comsps <- c("Aeropedellus clavatus", "Melanoplus boulderensis", "Camnula pellucida",
            "Melanoplus sanguinipes", "Chloealtis abdominalis", "Melanoplus dawsoni")
data <- data[data$species %in% comsps,]

# Subset to adults
data <- data[data$in6 > 0,]

#=============================
# Calculating pairwise overlap
#=============================

# Define sites, years, and species in sample
sites <- unique(data$site)
years <- sort(unique(data$year))
sps <- unique(data$species)

# Reorder species vector before defining all combinations
sps <- sps[c(1,6,5,2,4,3)]
spcmbs <- t(combn(sps, 2))

# Set up data storage
#overlap <- matrix(NA, nrow = nrow(spcmbs), ncol = length(years))
overlap <- array(NA, dim = c(nrow(spcmbs), length(years), length(sites)))
colnames(overlap) <- years

# Make species column a factor
data$species <- as.factor(data$species)

# Calculate pairwise overlap
for(site in 1:length(sites)) {
  
  # Subset to 1 site
  site.sub <- data[data$site == sites[site],]
  
  # Loop through years
  for(year in 1:length(years)) {
    
    # Subset to 1 year
    dat <- site.sub[site.sub$year == years[year],]
    
    # Create matrix of abundances for all species at each sampling event
    sumByOrd <- tapply(dat$in6, list(dat$species, dat$ordinal), FUN = sum)
    
    # Replace NAs with 0
    sumByOrd[is.na(sumByOrd)] <- 0
    
    # Loop through species combinations
    for(cmb in 1:nrow(spcmbs)) {
      
      # Calculate abundance of species 1 at each sampling period as proportion of
      # total abundance across season
      pij <- as.vector(sumByOrd[spcmbs[cmb,1],] / sum(sumByOrd[spcmbs[cmb,1],]))
      
      # Calculate same proportions for species 2
      pik <- as.vector(sumByOrd[spcmbs[cmb,2],] / sum(sumByOrd[spcmbs[cmb,2],]))
      
      # Calculate overlap and store result
      overlap[cmb, year, site] <- sum(pij * pik) / sqrt(sum(pij^2) * sum(pik^2))
      
    } # Close species combinations loop
  } # Close years loop
} # Close sites loop

# Replace NaN with NA
overlap[is.nan(overlap)] <- NA

# Convert to data frame
longOver <- as.data.frame(rbind(overlap[,,1], overlap[,,2], overlap[,,3]))
longOver$sp1 <- as.character(rep(spcmbs[,1], times = 3))
longOver$sp2 <- as.character(rep(spcmbs[,2], times = 3))
longOver$combID <- rep(1:15, times = 3)
longOver$site <- rep(sites, each = nrow(spcmbs))

# Reshape data for ggplot2
longOver <- melt(longOver, id.vars = c("site", "combID", "sp1", "sp2"),
                 variable.name = "year",
                 value.name = "overlap")

# Add column for resurvey/initial
longOver$period <- "resurvey"
longOver[longOver$year %in% c(1959, 1960), "period"] <- "initial"

# Plot pairwise overlap vs. year
ggplot(data = longOver, aes(x = year, y = overlap, col = site, shape = period)) +
  geom_point() +
  facet_grid(sp1 ~ sp2, drop = TRUE) +
  theme_bw()

#========================================
# Adding temperature data to overlap data
#========================================

# Subset climate data to sites and years for which we have grasshopper data
sites <- c("B1", "C1", "NOAA")
years <- c(1959, 1960, 2006:2015)
allClim <- droplevels(allClim[allClim$Site %in% sites & allClim$Year %in% years,])

# Subset to ordinal dates relevant for grasshopper season (March 1 to Aug 31 = ordinal 60 to 243)
allClim <- allClim[allClim$Julian > 59 & allClim$Julian < 244,]

# Check that all ordinal dates are included in the climate data
ords <- 60:243
results <- matrix(NA, nrow = length(sites), ncol = length(years))
rownames(results) <- sites
colnames(results) <- years
for(i in 1:length(sites)){
  for(j in 1:length(years)){
    clim <- allClim[allClim$Site == sites[i] & allClim$Year == years[j],]
    results[i,j] <- sum(ords %in% clim$Julian)/length(ords)
  }
}

# No climate data for C1 or NOAA (CHA) in 2015. Otherwise all ordinal dates are 
# present for all sites in all years

# Calculate GDDs accumulated on each day
clim <- allClim[!is.na(allClim$Min) & !is.na(allClim$Max),]
clim$dd <- apply(clim[,c("Min","Max")], MARGIN = 1, FUN = degree.days.mat, 
                 LDT= 12)
allClim <- merge(allClim, clim, all.x = T)

# Check dates for which GDDs could not be calculated
allClim[is.na(allClim$dd),]

# Longest run of days without GDDs is 3. Interpolation should not lead to
# major errors

# Re-order allClim
allClim <- allClim[order(allClim$Site, allClim$Year, allClim$Julian),]

# Obtain rows which need to be interpolated
x <- which(is.na(allClim$dd))

# Interpolate data
allClim$order <- 1:nrow(allClim)
allClim$dd <- na.approx(allClim$dd, allClim$order, na.rm = F)

# Check interpolation
allClim[x,]

# Mistake for B1 2012 ordinal 60 because the interpolation is done between the
# last day of 2011 and the second day of 2012. Change this value to zero manually
allClim[x[1],"dd"] <- 0

# Extract annual measures of GDD accumulation for each site. The measures are the
# GDDs accumulated in 4 periods: ordinals 60-135, 136-181, 182-212, and 213-243.
accumPeriod <- c("ords60_135", "ords136_181", "ords182_212", "ords213_243")
allClim$accumPeriod <- NA
allClim[allClim$Julian %in% 60:135, "accumPeriod"] <- accumPeriod[1]
allClim[allClim$Julian %in% 136:181, "accumPeriod"] <- accumPeriod[2]
allClim[allClim$Julian %in% 182:212, "accumPeriod"] <- accumPeriod[3]
allClim[allClim$Julian %in% 213:243, "accumPeriod"] <- accumPeriod[4]
GDDdata <- allClim %>% group_by(Site, Year, accumPeriod) %>%
  summarise(GDDs = sum(dd))

# Change to wide format
GDDdata <- dcast(GDDdata, Site + Year ~ accumPeriod, value.var = "GDDs")

# Add GDD data to overlap data
GDDdata$Site <- as.character(GDDdata$Site)
GDDdata[GDDdata$Site == "NOAA", "Site"] <- "CHA"
finalDat <- merge(longOver, GDDdata, by.x = c("site", "year"), by.y = c("Site", "Year"),
      all.x = TRUE)

#======================================================
# Which GDD accumulation window seems most interesting?
#======================================================

ggplot(data = finalDat, #aes(x = ords60_135,
#                       aes(x = ords136_181,
                        aes(x = ords182_212,
#                        aes(x = ords213_243, 
                            y = overlap, color = site)) +
  geom_point(aes(shape = period)) +
  facet_grid(sp1 ~ sp2, drop = TRUE) +
  theme_bw() + 
  geom_smooth(method=lm) # this code makes separate models for initial and resurvey

# ords 182-212 seems the most interesting

# 2009 at B1 has suspiciously low GDD accumulation in all windows. Nothing looks
# strange in the max and min data, it just seems to be a very cold year
