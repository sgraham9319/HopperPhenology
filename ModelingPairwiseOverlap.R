#############################################################
# Modeling pairwise overlap as a function of GDD accumulation
#############################################################

# This analysis used the data frame named "finalDat" which is the output of the
# R script named "PairwiseOverlapVsTemp.R"

# Define sites
sites <- sort(unique(finalDat$site))

# Create results storage
modelOutput <- array(NA, dim = c(length(unique(finalDat$combID)), 4, length(sites)))
colnames(modelOutput) <- c("slope", "SE", "tvalue", "pvalue")

for(site in 1:length(sites)){
  for(com in 1:length(unique(finalDat$combID))){
    
    # Subset to 1 site and 1 species combination
    dat <- finalDat[finalDat$site == sites[site] & finalDat$combID == com,]
    
    # Remove year 2009 for B1
    if(sites[site] == "B1"){
      dat <- dat[dat$year != 2009,]
    }
    
    # Only try to create model if data exists
    if(sum(!is.na(dat$overlap)) > 0){
      # Create linear model
      mod <- lm(overlap ~ ords182_212, data = dat)
      
      # Extract coefficients table
      modelOutput[com,,site] <- summary(mod)[["coefficients"]][2,]
    }
  }
}

# Combine all output in a data frame
combIDs <- finalDat[match(sort(unique(finalDat$combID)), finalDat$combID), c("combID", "sp1", "sp2")]
coefs <- rbind(modelOutput[,,1], modelOutput[,,2], modelOutput[,,3])
IDs <- rbind(combIDs, combIDs, combIDs)
modelResults <- cbind(IDs, coefs)
modelResults$site <- rep(sites, each = 15)

# View marginally significant results
View(modelResults[!is.na(modelResults$pvalue) & modelResults$pvalue < 0.1,])

# Write results to csv
#write.csv(modelResults, "PairwiseOverlapBySpsComb.csv", row.names = F)

#===========
# Full model
#===========

# Assigning life histories
# Levels = "early", "late", "mix"
# early = 1, 2, 6 late = 13, 14, 15 mix = 3,4,5,7,8,9,10,11,12
#lifeHistDat <- data.frame(1:15, c(rep("early", times = 2), rep("mix", times = 3),
#                                  "early", rep("mix", times = 6), rep("late", times = 3)))
#names(lifeHistDat) <- c("site", "lifeHist")
#
#test <- merge(finalDat, lifeHistDat, by.x = site)

# Full model
#mod <- lm(overlap ~ site * ords213_243 + lifeHist * site, data = test)