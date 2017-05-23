#################################
# Phenology vs. timing of warming
#################################

library(ggplot2)
library(dplyr)
library(reshape2)

#=================================
# Load and format grasshopper data
#=================================

# Load grasshopper data
setwd("C:/Users/stuart/Google Drive/DataForAnalysis/grasshoppers/SexCombined")
#data <- read.csv("C1_1959-2015_eggDiapause.csv")
#data <- read.csv("B1_1959-2015_eggDiapause.csv")
data <- read.csv("CHA_1959-2012_eggDiapause.csv")

# Check for missing data
sum(is.na(data$ordinal))
sum(is.na(data$in6))
sum(is.na(data$year))

# Remove any spaces accidentally added at end of species names
data$species <- gsub(" $", "", data$species)
sort(unique(data$species))

# Subset to adults
data <- data[data$in6 > 0,]

# Subset to common species (B1 and CHA only)
# B1
comsps <- c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus boulderensis",
            "Melanoplus dawsoni", "Melanoplus packardii")
data <- data[data$species %in% comsps,]

# CHA
comsps <- c("Aeropedellus clavatus", "Hesperotettix viridis", "Melanoplus bivittatus",
            "Melanoplus dawsoni", "Melanoplus sanguinipes")
data <- data[data$species %in% comsps,]

#======================
# Quantifying phenology
#======================

# Define species and years in sample
sps <- as.character(sort(unique(data$species)))
years <- sort(unique(data$year))

# Make species column a factor
data$species <- as.factor(data$species)

# Set up data storage
results <- array(NA, dim = c(length(sps), 3, length(years)))
dimnames(results)[[1]] <- sps
dimnames(results)[[2]] <- c("FirstApp", "PeakAbund", "SeasonLen")

# Calculate phenology metrics
for(year in 1:length(years)){ # Loop through years
  
  # Subset to 1 year
  dat <- data[data$year == years[year],]
  
  # Create matrix of abundances for all species at each sampling event
  sumByOrd <- tapply(dat$in6, list(dat$species, dat$ordinal), FUN = sum)
  
  # Replace NAs with 0
  sumByOrd[is.na(sumByOrd)] <- 0
  
  for(sp in 1:length(sps)) { # Loop through species
    if(sum(sumByOrd[sp,]) > 0) { # Check species was observed
      
      # Get date of first adult observation
      results[sp,"FirstApp",year] <- 
        as.numeric(names(which(sumByOrd[sp,] > 0)[1]))
      
      # Get date of peak adult abundance
      results[sp,"PeakAbund",year] <- 
        weighted.mean(as.numeric(colnames(sumByOrd)), w = sumByOrd[sp,])
      
      # Get season length
      results[sp,"SeasonLen",year] <- 
        max(as.numeric(names(which(sumByOrd[sp,] > 0)))) - results[sp,"FirstApp",year]
    } # Close if statement
  } # Close species loop
} # Close years loop

# Combine into long data frame
pheno <- as.data.frame(results[,,1])
for(yr in 2:length(years)){
  pheno <- rbind(pheno, results[,,yr])
}
pheno$species <- as.factor(rep(sps, times = length(years)))
pheno$year <- rep(years, each = length(sps))
rownames(pheno) <- NULL

#======================================
# Adding temperature data to data frame
#======================================

# Load climate data
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
clim <- clim[!is.na(clim$Min) & !is.na(clim$Max),]

# Calculate GDDs accumulated on each day
setwd("C:/Users/stuart/Google Drive/University of Washington/Buckley lab/Phenological overlap project/HopperPhenology")
source("degreedays.R")
clim$dd <- apply(clim[,c("Min","Max")], MARGIN = 1, FUN = degree.days.mat, 
                 LDT= 12)

# Calculate cumulative number of GDDs accumulated by each day in each year
# (ignoring GDDs accumulated before March 1st)
clim <- clim %>% arrange(Year, Julian) %>% group_by(Year) %>% 
  mutate(cdd = cumsum(dd))

# Extract different temperature metrics for each year
#ords <- c(160, 200, 240)
#tempMeas <- clim[clim$Julian %in% ords, c("Year", "Julian", "cdd")]
dd141_160 <- clim[clim$Julian > 140 & clim$Julian <= 160,] %>% 
  group_by(Year) %>% summarize(dds = sum(dd))
dd161_180 <- clim[clim$Julian > 160 & clim$Julian <= 180,] %>% 
  group_by(Year) %>% summarize(dds = sum(dd))
dd181_200 <- clim[clim$Julian > 180 & clim$Julian <= 200,] %>% 
  group_by(Year) %>% summarize(dds = sum(dd))
dd201_220 <- clim[clim$Julian > 200 & clim$Julian <= 220,] %>% 
  group_by(Year) %>% summarize(dds = sum(dd))
dd221_240 <- clim[clim$Julian > 220 & clim$Julian <= 240,] %>%
  group_by(Year) %>% summarize(dds = sum(dd))

# Create a data frame of the different temperature metrics
Temp <- data.frame(dd141_160$Year,
                   #tempMeas[tempMeas$Julian == 160, "cdd"],
                   #tempMeas[tempMeas$Julian == 200, "cdd"],
                   #tempMeas[tempMeas$Julian == 240, "cdd"],
                   dd141_160$dds,
                   dd161_180$dds,
                   dd181_200$dds,
                   dd201_220$dds,
                   dd221_240$dds)
names(Temp) <- c("Year", #"GDDa160", "GDDa200", "GDDa240", 
                 "Ords141-160", "Ords161-180", "Ords181-200",
                 "Ords201-220", "Ords221-240")

# Add temp data to phenology data
pheno <- pheno[pheno$year %in% Temp$Year,]
Temp <- Temp[match(pheno$year, Temp$Year),]
plotPheno <- cbind(pheno, Temp[,-which(names(Temp) == "Year")])

#==================
# Plotting the data
#==================

# Find out which are early and late season species
ggplot(data = plotPheno, aes(x = year, y = PeakAbund, col = species)) +
  geom_point()

# Put data in long format for plotting
pheno1 <- melt(plotPheno, id.vars = 1:5, value.name = "GDDsAccumulated")

# Add column for resurvey/initial
pheno1$period <- "resurvey"
pheno1[pheno1$year %in% c(1959, 1960), "period"] <- "initial"

# Re-order factor levels of species columns so that plot shows from
# left to right and top to bottom, early-season to late-season species
pheno1$species <- as.factor(pheno1$species)
# B1
#pheno1$species <- factor(pheno1$species, levels = c("Aeropedellus clavatus", "Melanoplus boulderensis", "Camnula pellucida", "Melanoplus packardii", "Melanoplus dawsoni"))
# C1
#pheno1$species <- factor(pheno1$species, levels = c("Melanoplus bouderensis", "Melanoplus fasciatus", "Melanoplus sanguinipes", "Camnula pellucida", "Chloealtis abdominalis"))
# CHA
pheno1$species <- factor(pheno1$species, levels = c("Aeropedellus clavatus", "Hesperotettix viridis", "Melanoplus bivittatus", "Melanoplus sanguinipes", "Melanoplus dawsoni"))

# Plot data
ggplot(data = pheno1, aes(x = GDDsAccumulated, y = SeasonLen, col = period)) +
  geom_point() +
  facet_grid(variable ~ species) +
  theme_bw()


