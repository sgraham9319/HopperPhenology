#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

#--------------------------------------
#create path to data directory
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:/Users/stuart/Google Drive/DataForAnalysis/"

#load climate data
setwd(paste(fdir, "climate", sep= ""))

clim.a1= read.csv("Climate_A1_1953-2014.csv")
clim.b1= read.csv("Climate_B1_1953-2014.csv")
clim.c1= read.csv("Climate_C1_1953-2014.csv")
clim.d1= read.csv("Climate_D1_1953-2014.csv")
clim.noaa= read.csv("Climate_noaa_1953-2014.csv")

#add site column
clim.a1$Site= "A1"
clim.b1$Site= "B1"
clim.c1$Site= "C1"
clim.d1$Site= "D1"
clim.noaa$Site= "NOAA"

#combine ##currently omits flags
cols= c("Site", "Date", "Year", "Month", "Julian", "Max", "Min", "Mean")
clim= rbind(clim.a1[,cols], clim.b1[,cols], clim.c1[,cols], clim.d1[,cols], 
            clim.noaa[,cols])

#==================================
#add recent climate data from Cesar
#==================================

setwd(paste(fdir, "climate\\Recent\\", sep=""))  

#B1
clim.b1.2015= read.csv("B1_2015.csv")
#calculate Julian
tmp <- as.POSIXlt(clim.b1.2015$Date.Time..GMT.07.00, format = "%m/%d/%Y %H:%M")
clim.b1.2015$Julian= tmp$yday+1
#calc min and max
clim.b1.2015= clim.b1.2015 %>% group_by(Year, Julian) %>% 
  summarise(Date=Date.Time..GMT.07.00[1], Min= min(Temp_C), Max= max(Temp_C))
clim.b1.2015= as.data.frame(clim.b1.2015)
clim.b1.2015$Site= "B1"
clim.b1.2015$Mean= NA
clim.b1.2015$Month= NA
clim= rbind(clim, clim.b1.2015[,cols])

clim.b1.2009.tmax= read.csv("B1_2009_Tmax.csv")
clim.b1.2009.tmin= read.csv("B1_2009_Tmin.csv")
clim.b1.2010.tmax= read.csv("B1_2010_Tmax.csv")
clim.b1.2010.tmin= read.csv("B1_2010_Tmin.csv")
clim.b1.2012.tmax= read.csv("B1_2012_Tmax.csv")
clim.b1.2012.tmin= read.csv("B1_2012_Tmin.csv")

clim.b1.2009.tmax$Min= clim.b1.2009.tmin$Min
clim.b1.2010.tmax$Min= clim.b1.2010.tmin$Min
clim.b1.2012.tmax$Min= clim.b1.2012.tmin$Min

clim.b1.20122013= rbind(clim.b1.2009.tmax, clim.b1.2010.tmax, clim.b1.2012.tmax)
clim.b1.20122013$Site= "B1"
clim.b1.20122013$Date= NA
clim.b1.20122013$Mean= NA
clim= rbind(clim, clim.b1.20122013[,cols])

#C1
clim.c1.20092012= read.csv("C1_2009_2012.csv")
#add missing data from 2011
JulianNAs <- clim[clim$Site == "C1" & clim$Year == 2011 & is.na(clim$Max), "Julian"]
newDat <- clim.c1.20092012[clim.c1.20092012$Year == 2011 & clim.c1.20092012$Julian %in% JulianNAs,]
clim[clim$Site == "C1" & clim$Year == 2011 & clim$Julian %in% JulianNAs, "Max"] = newDat$Max
clim[clim$Site == "C1" & clim$Year == 2011 & clim$Julian %in% JulianNAs, "Min"] = newDat$Min
#add missing data from 2010
JulianNAs <- clim[clim$Site == "C1" & clim$Year == 2010 & is.na(clim$Max), "Julian"]
newDat <- clim.c1.20092012[clim.c1.20092012$Year == 2010 & clim.c1.20092012$Julian %in% JulianNAs,]
clim[clim$Site == "C1" & clim$Year == 2010 & clim$Julian %in% JulianNAs, "Max"] = newDat$Max
clim[clim$Site == "C1" & clim$Year == 2010 & clim$Julian %in% JulianNAs, "Min"] = newDat$Min
# add 2012 data
clim.c1.2012= clim.c1.20092012[clim.c1.20092012$Year > 2011,]
clim.c1.2012$Site= "C1"
clim= rbind(clim, clim.c1.2012[,cols])

#add 2014 data
clim.c1.2014= read.csv("C1_2014.csv")
clim.c1.2014$Site= "C1"
clim.c1.2014$Month= NA
clim= rbind(clim, clim.c1.2014[,cols])

#=============
#add LTER data
#=============

setwd(paste(fdir, "climate\\LTERdownload_Recent\\", sep= ""))  

#A1
a1= read.csv("a-1hobo.hourly.jm.data.csv") #2013-2014
#calc min and max
a1= a1 %>% group_by(year, jday) %>% summarise(Date= date[1], 
                                              Min= min(airtemp_avg), 
                                              Max= max(airtemp_avg))
a1= as.data.frame(a1)
a1$Site= "A1"
a1$Year= a1$year
a1$Julian= a1$jday
a1$Mean= NA
a1$Month= NA

clim= rbind(clim, a1[,cols])

#B1
b1= read.csv("b-1hobo.hourly.jm.data.csv") #2013-2014
#calc min and max
b1= b1 %>% group_by(year, jday) %>% summarise(Date= date[1], 
                                              Min= min(airtemp_avg), 
                                              Max= max(airtemp_avg))
b1= as.data.frame(b1)
b1$Site= "B1"
b1$Year= b1$year
b1$Julian= b1$jday
b1$Mean= NA
b1$Month= NA
b1201314= b1[b1$Year == 2013 | b1$Year == 2014,]

clim= rbind(clim, b1201314[,cols])

#C1
c1= read.csv("c-1tdayv.ml.data.csv")
#calculate Julian
tmp <- as.POSIXlt(c1$Date, format = "%m/%d/%Y")
c1$Julian= tmp$yday
c1$Site= "C1"
c1$Month= NA
#subset to missing years
c12013= c1[c1$Year == 2013,]
c12014= c1[c1$Year == 2014 & c1$Julian >= 335,]
clim= rbind(clim, c12013[,cols], c12014[,cols])

#D1
d1= read.csv("d-1cr23x-cr1000.daily.ml.data_up.csv")
d1$Site= "D1"
d1$Month= NA
#subset to missing years
d1= subset(d1, d1$Year > 2008)

clim= rbind(clim, d1[,cols])

#=====================
# Clean and check data
#=====================

#replace NaNs
clim[clim == "NaN"] = "NA"

#change temps to numeric
typeof(clim$Min)
clim[,"Min"]= as.numeric(clim[,"Min"])
clim[,"Max"]= as.numeric(clim[,"Max"])
clim[,"Mean"]= as.numeric(clim[,"Mean"])

#check counts of data
years= c(1958:1960, 2006:2015)

climsub= clim[clim$Year %in% years,]
#counts across sites years
climsub= climsub %>% group_by(Site, Year) %>% 
  summarise(Min= length(na.omit(Min)), Max= length(na.omit(Max)), 
            Mean= length(na.omit(Mean)))
climsub= as.data.frame(climsub)

#=========
#Plot data
#=========

#summer means
clim1= clim[which(clim$Julian > 59 & clim$Julian < 244),]
clim1= clim1 %>% group_by(Year, Site) %>% summarise(Min= mean(Min, na.rm= TRUE),
                                                    Max= mean(Max, na.rm= TRUE),
                                                    Mean= mean(Mean, na.rm= TRUE))

ggplot(data= clim1, aes(x= Year, y= Min, color= Site )) + geom_line() +
  theme_bw()

#===============
# Write out data
#===============

setwd(paste(fdir, "climate", sep= "")) 
write.csv(clim, "AlexanderClimateAll.csv", row.names = FALSE)
