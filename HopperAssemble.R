
# Load packages
library(dplyr)
# see dplyr introduction at: https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html

# IMPORTANT: if plyr package is loaded after dplyr package, some functions of
# dplyr will work incorrectly but no warning will be given!

# Set working directory
setwd("C:/Users/stuart/Google Drive/DataForAnalysis/grasshoppers/SexCombined")

# Load data
a1 <- read.csv("A1_1958-2011_allSpecies.csv")
b1 <- read.csv("B1_1959-2015_eggDiapause.csv")
c1 <- read.csv("C1_1959-2015_eggDiapause.csv")
cha <- read.csv("CHA_1959-2012_eggDiapause.csv")


test <- a1[1:10,]
test[,"species"] <- as.factor(c(rep("Melanoplus dodgei", times = 3), 
                                rep("Trimerotropis cincta", times = 3), 
                                rep("Melanoplus sanguinipes", times = 3), 
                                "Melanoplus dodgei"))

a2 <- group_by(test, species)
summarise(a2, x = sum(in6))#, in5 = sum(in5), in4 = sum(in4), in3 = sum(in3),
          #in2 = sum(in2), in1 = sum(in1))
apply(test[5:10], MARGIN = 2, FUN = sum)
typeof(a2)
