#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

#--------------------------------------
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   

clim.a1= read.csv( "Climate_A1_1953-2014.csv" )
clim.b1= read.csv( "Climate_B1_1953-2014.csv" )
clim.c1= read.csv( "Climate_C1_1953-2014.csv" )
clim.d1= read.csv( "Climate_D1_1953-2014.csv" )
clim.noaa= read.csv( "Climate_noaa_1953-2014.csv" )

#add site
clim.a1$Site="A1"
clim.b1$Site="B1"
clim.c1$Site="C1"
clim.d1$Site="D1"
clim.noaa$Site="NOAA"

#combine ##currently omits flags
clim= rbind(clim.a1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.b1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.c1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.d1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.noaa[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#------
#load recent climate data from Cesar
setwd( paste(fdir, "climate\\Recent\\", sep="") )  

#columns: "Site","Date","Year","Month","Julian","Max","Min","Mean"

#A1
clim.a1.2011= read.csv( "A1_2011.csv" )
clim.a1.2011$Site="A1"
clim.a1.2011$Date=NA
clim= rbind(clim, clim.a1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#B1
clim.b1.2011= read.csv( "B1_2011.csv" )
clim.b1.2011$Site="B1"
clim.b1.2011$Date=NA
clim= rbind(clim, clim.b1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

clim.b1.2015= read.csv( "B1_2015.csv" )
#calculate Julian
tmp <- as.POSIXlt(clim.b1.2015$Date.Time..GMT.07.00, format = "%m/%d/%Y %H:%M")
clim.b1.2015$Julian= tmp$yday+1
#calc min and max
clim.b1.2015= clim.b1.2015 %>% group_by(Year,Julian) %>% summarise(Date=Date.Time..GMT.07.00[1],Min= min(Temp_C),Max= max(Temp_C) )
clim.b1.2015= as.data.frame(clim.b1.2015)
clim.b1.2015$Site="B1"
clim.b1.2015$Mean=NA
clim.b1.2015$Month=NA
clim= rbind(clim, clim.b1.2015[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

clim.b1.2009.tmax= read.csv( "B1_2009_Tmax.csv" )
clim.b1.2009.tmin= read.csv( "B1_2009_Tmin.csv" )
clim.b1.2010.tmax= read.csv( "B1_2010_Tmax.csv" )
clim.b1.2010.tmin= read.csv( "B1_2010_Tmin.csv" )
clim.b1.2012.tmax= read.csv( "B1_2012_Tmax.csv" )
clim.b1.2012.tmin= read.csv( "B1_2012_Tmin.csv" )
clim.b1.2013.tmax= read.csv( "B1_2012_Tmax.csv" )

clim.b1.2009.tmax$Min= clim.b1.2009.tmin$Min
clim.b1.2010.tmax$Min= clim.b1.2010.tmin$Min
clim.b1.2012.tmax$Min= clim.b1.2012.tmin$Min
clim.b1.2013.tmax$Min=NA

clim.b1.20122013= rbind(clim.b1.2009.tmax, clim.b1.2010.tmax, clim.b1.2012.tmax, clim.b1.2013.tmax)
clim.b1.20122013$Site="B1"
clim.b1.20122013$Date=NA
clim.b1.20122013$Mean=NA
clim= rbind(clim, clim.b1.20122013[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#C1
clim.c1.20092012= read.csv( "C1_2009_2012.csv" )
clim.c1.2014= read.csv( "C1_2014.csv" )

clim.c1.2012= subset(clim.c1.20092012, clim.c1.20092012$Year>2011)
clim.c1.2012$Site="C1"
clim= rbind(clim, clim.c1.2012[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

clim.c1.2014$Site="C1"
clim.c1.2014$Month=NA
clim= rbind(clim, clim.c1.2014[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#------
#load and process LTER data
setwd( paste(fdir, "climate\\LTERdownload_Recent\\", sep="") )  

#A1
a1= read.csv("a-1hobo.hourly.jm.data.csv") #2013-2014
#calc min and max
a1= a1 %>% group_by(year,jday) %>% summarise(Date=date[1],Min= min(airtemp_avg),Max= max(airtemp_avg) )
a1= as.data.frame(a1)
a1$Site="A1"
a1$Year=a1$year
a1$Julian=a1$jday
a1$Mean=NA
a1$Month=NA

clim= rbind(clim, a1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#B1
b1= read.csv("b-1hobo.hourly.jm.data.csv") #2013-2014
#calc min and max
b1= b1 %>% group_by(year,jday) %>% summarise(Date=date[1],Min= min(airtemp_avg),Max= max(airtemp_avg) )
b1= as.data.frame(b1)
b1$Site="B1"
b1$Year=b1$year
b1$Julian=b1$jday
b1$Mean=NA
b1$Month=NA

clim= rbind(clim, b1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#C1
c1= read.csv("c-1tdayv.ml.data.csv")
#calculate Julian
tmp <- as.POSIXlt(c1$Date, format = "%m/%d/%Y")
c1$Julian= tmp$yday
c1$Site="C1"
c1$Month=NA
#subset to missing years
c1= subset(c1, c1$Year %in% c(2013,2015))

clim= rbind(clim, c1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#D1
d1= read.csv("d-1cr23x-cr1000.daily.ml.data_up.csv")
d1$Site="D1"
d1$Month=NA
#subset to missing years
d1= subset(d1, d1$Year>2008)

clim= rbind(clim, d1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#=====================================
# Remove duplicates

clim$Site_Year_Julian= paste(clim$Site, clim$Year, clim$Julian, sep="")
# average across duplicates
clim1= clim %>% group_by(Site,Date,Year,Month,Julian) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )
clim1= as.data.frame(clim1)

#replace NaN
clim1[clim1 == "NaN"] = "NA"
#change temps to numeric
clim1[,"Min"]= as.numeric(clim1[,"Min"])
clim1[,"Max"]= as.numeric(clim1[,"Max"])
clim1[,"Mean"]= as.numeric(clim1[,"Mean"])

#------------------------
#Plot

#summer means
clim2= clim1[which(clim1$Julian>59 & clim1$Julian<244),]
clim2 = clim2 %>% group_by(Year,Site) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )

ggplot(data=clim2, aes(x=Year, y = Min, color=Site ))+geom_line() +theme_bw()

#=====================================
# WRITE OUT DATA
setwd( paste(fdir, "climate", sep="") )   
write.csv(clim1, "AlexanderClimateAll.csv")

#check counts of data
years= c(1958:1960, 2006:2015)

climsub=clim1[clim1$Year %in% years,]
#counts across sites years
climsub= climsub %>% group_by(Site,Year) %>% summarise(Min= length(na.omit(Min)),Max= length(na.omit(Max)),Mean= length(na.omit(Mean)) )
climsub= as.data.frame(climsub)

