#####################################################################
# start with clean slate: 
# rm(list = ls())

#####################################################################
#bring in libraries
library(CPUE)
library(stringr) #for str_trim
library(shiny)
library(reshape)
library(ncdf.tools)
library(ncdf4)
library(xtable)
library(shapefiles)
library(sp)
library(plyr)
library(dplyr)
library(tidyr)
library(markdown)
library(ggplot2)
library(DT)
library(scales)
# library(XLConnect)
library(proj4)a
library(plyr)
library(jsonlite)
library(lazyeval)
library(ncdf4)
library(binr)
library(tensorA)
library(tidyverse)
library(openxlsx)
library(ggthemes)
library(Hmisc)
library(hydroGOF)
library(corrgram)
library(plotrix)
library(ggmap)
library(rgdal)
library(gtable)
library(maps)
library(mapdata)
library(rgeos)
library(ggsubplot)
library(raster)
library(tree)

# Set paths and working directories

# Vidette:
thisPath<-"C:\\projects\\2015\\FIFI1501\\ATLANTISChathamRise\\ATLANTISChathamRiseRepository\\"

thisPath<-"C:\\projects\\2017\\ATLANTISChathamRiseRepository\\"

source(paste(thisPath,"R\\(0)Functions_general\\make.filename.R",sep=""))
source(paste(thisPath,"R\\(0)Functions_general\\assign.directories.R",sep=""))
DIR<-assign.directories(base=thisPath)
# tidyup
rm(assign.directories,make.filename)

#Source other functions needed
source(paste(DIR$'General functions',"get_first_number.R",sep=""))
source(paste(DIR$'General functions',"makeBlankPlot.R",sep=""))
source(paste(DIR$'General functions',"nonZeroMean.R",sep=""))
source(paste(DIR$'General functions',"formatShape.R",sep=""))
source(paste(DIR$'General functions',"myRounding.R",sep=""))
source(paste(DIR$'General functions',"multiplot.R",sep=""))

thisCex<-1.5

#set colours
myRed<-rgb("red"=0.8,"blue"=0.5,"green"=0)
myRed_trans<-rgb("red"=0.8,"blue"=0.5,"green"=0,"alpha"=0.3)
myBlue<-rgb("red"=0,"blue"=0.8,"green"=0.5)
myBlue_trans<-rgb("red"=0,"blue"=0.8,"green"=0.5,"alpha"=0.3)
myGreen<-rgb("red"=0.5,"blue"=0,"green"=0.8)
myGreen_trans<-rgb("red"=0.5,"blue"=0,"green"=0.8,"alpha"=0.3)
myOrange<-rgb("red"=1,"blue"=0,"green"=0.7)
myOrange_trans<-rgb("red"=1,"blue"=0,"green"=0.7,"alpha"=0.3)
myPurple<-rgb("red"=0.63,"blue"=0.94,"green"=0.12)
myPurple_trans<-rgb("red"=0.63,"blue"=0.94,"green"=0.12,"alpha"=0.3)
myGrey<-rgb("red"=0.5,"blue"=0.5,"green"=0.5)
myGrey_trans<-rgb("red"=0.5,"blue"=0.5,"green"=0.5,"alpha"=0.3)
myYellow<-rgb("red"=1,"blue"=0.01,"green"=1)
myYellow_trans<-rgb("red"=1,"blue"=0.01,"green"=1,"alpha"=0.3)
myAqua<-rgb("red"=0.30,"green"=0.89,"blue"=0.85)
myLightGrey<-rgb("red"=0.7,"green"=0.7,"blue"=0.7)
myVeryLightGrey<-rgb("red"=0.9,"green"=0.9,"blue"=0.9)
myLightBlue<-rgb("red"=0.7,"blue"=0.91,"green"=0.85)

myLightAqua<-rgb("red"=0.9,"green"=0.984,"blue"=0.98)
myDarkGreen<-rgb("green"=0.5,"blue"=0.25,"red"=0.125)
myMidnightGreen<-rgb("green"=0.2,"blue"=0.1,"red"=0.05)
myLightGreen<-rgb("red"=0.9,"blue"=0.62,"green"=0.99)

myGold<-rgb("red"=1, "green"=0.85, "blue"=0)
myDarkRed<-rgb("red"=0.7,"blue"=0.1,"green"=0)

mg_2_tonne<-2e-8; X_CN<-5.7
