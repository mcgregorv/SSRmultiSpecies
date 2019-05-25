#ratio of TL 4 and up over TL 3 biomass - outputs created in ecologicalIndicatorsFromSensRunsSETUPandSTORE.r
Version<-"C"
outFolder<-"SSRsens"
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep=""); basePath<-thisPath
dataOutPath<-paste(thisPath,"SSRsens\\SSRsensTracers\\CatchScenarios_V",Version,sep="")
plotPath<-paste(DIR$'Base',"ATLANTISmodels\\Figures\\Recruitment\\Paper2\\V",Version,sep="")
# plotPath<-paste(DIR$'Figures',"Recruitment\\postReview\\", sep="") ## write up

# thisOutList<-c("TLratio", "Qtimes", "TLtimes_df", "meanLengthByRun","biomOverCatch", "ratioPelBiom", "sensZG")
# save(list=thisOutList,file=paste(dataOutPath,"ecoIndicators",sep=""))
load(paste(dataOutPath,"ecoIndicators", sep=""))

TLarray<-Qtimes
lookup_df$scenario<-lookup_df$Scenario

allScenarios<-unique(lookup_df$Scenario); nscenarios<-length(allScenarios)
scenarioColors<-colorRampPalette(colors=c(myOrange,"red",myPurple,myBlue))(nscenarios)

thisMax<-max(TLarray,na.rm=TRUE)*1.1
colBySteepness<-colorRampPalette(colors=c(myGold, myGreen,myDarkGreen, "black"))(3)  
steepnessSens<-c(0.5,0.7,0.9)

QtimeYears<-seq(2010,2046); QtimeIndex<-QtimeYears-1865 
axisYears<-QtimeYears[seq(1,length(QtimeYears),by=5)]; axisYearsAt<-seq(1,length(QtimeYears), by=5)
snapshotYears<-c(2020,2025,2035,2045); snapshotAt<-match(snapshotYears, QtimeYears)

histYears<-seq(1970,2016); histYearsIndex<-histYears-1865
axisHistYears<-histYears[seq(1,length(histYears),by=5)]; axisHistYearsAt<-seq(1,length(histYears), by=5)
thisMin <- 0
preFishingMin<-min(Qtimes[,1:min(QtimeIndex)], na.rm=TRUE); preFishingMax<-max(Qtimes[,1:min(QtimeIndex)], na.rm=TRUE)
for(r in 1:2){
  for(h in 1:3){
    if(r==1){
      baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"No cap"
    } else{
      baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"Cap"
    }
    pdf(paste(plotPath, "KemptonsQScenario_h",h,"_s",r,".pdf", sep=""), height=3, width=4.5)
    par(mar=c(4,4,1.5,0.5))
    plot(Qtimes[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Kempton's Q")
    axis(at=axisYearsAt,labels = axisYears,side=1,las=2)
    polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
            border=NA, col=myGreen_trans)
    for(s in 1:nscenarios){
      if(s!=3){
        cat(lookup_df$scenario[baseCappedIndex][s],"--")
        points(Qtimes[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
      }
    }
    legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="bottomleft")
    abline(v=grep("2016",QtimeYears), col="red", lty=4)
    mtext(paste(capText,", h=", htext, sep=""), side=3,  adj=0)
    dev.off()
  }
}
thisCex<-1.5
## just h=0.5, for figure in text body
pdf(paste(plotPath, "KemptonsQScenario_h",h,".pdf", sep=""), height=4.5, width=15)
par(mar=c(4,4.5,1.5,0.5), mfrow=c(1,2))
h=1
r<-1
baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
capText<-"No cap"
  plot(Qtimes[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Kempton's Q", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  axis(at=axisYearsAt,labels = axisYears,side=1,las=2, cex.axis=thisCex, cex.lab=thisCex)
  polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
          border=NA, col=myGrey_trans)
  for(s in 1:nscenarios){
    if(s!=3){
      cat(lookup_df$scenario[baseCappedIndex][s],"--")
      points(Qtimes[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
    }
  }
  legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="bottomright", cex=thisCex)
  abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste("A: ",capText,", h=", htext, sep=""), side=3,  adj=0, font=2, cex=thisCex)
  
r<-2
  baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
  capText<-"Cap"
  plot(Qtimes[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Kempton's Q", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  axis(at=axisYearsAt,labels = axisYears,side=1,las=2, cex.axis=thisCex, cex.lab=thisCex)
  polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
          border=NA, col=myGrey_trans)
  for(s in 1:nscenarios){
    if(s!=3){
      cat(lookup_df$scenario[baseCappedIndex][s],"--")
      points(Qtimes[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
    }
  }
  legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="bottomright", cex=thisCex)
  abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste("B: ",capText,", h=", htext, sep=""), side=3,  adj=0, font=2, cex=thisCex)
dev.off()

## do historic period comparing h for each SSR

## do historic period comparing h for each SSR
catchScenIndex <- lookup_df$scenario=="All100catch"
## all one figure
pdf(paste(plotPath, "KemptonsQ_all.pdf", sep=""), height=15, width=13)
par(mar=c(4.5,4.5,2,1), mfrow=c(4,2))

for(r in 1:2){
  if(r==1){
    baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$scenario=="All100catch";  
    capText<-"No cap"
  } else{
    baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$scenario=="All100catch";  
    capText<-"Cap"
  }
  plot(TLarray[1,histYearsIndex],type="n", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Kempton's Q", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  axis(at=axisHistYearsAt,labels = axisHistYears,side=1,las=2, cex.axis=thisCex, cex.lab=thisCex)
  # polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
  #         border=NA, col=myGreen_trans)
  for(h in 1:3){
    points(TLarray[baseCappedIndex,histYearsIndex][h,],type="l", col=colBySteepness[h])
  }
  legend(legend=steepnessSens, lty=1, seg.len=3, bty="n", col=colBySteepness, x="bottomright", inset=0.01, title="Steepness", cex=thisCex)
  # abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste(capText, sep=""), side=3,  adj=0, cex=thisCex)

}

  for(h in 1:3){
    for(r in 1:2){
      
    if(r==1){
      baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"No cap"
    } else{
      baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"Cap"
    }
    plot(Qtimes[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Kempton's Q", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
    axis(at=axisYearsAt,labels = axisYears,side=1,las=2, cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
    polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
            border=NA, col=myGrey_trans)
    for(s in 1:nscenarios){
      if(s!=3){
        cat(lookup_df$scenario[baseCappedIndex][s],"--")
        points(Qtimes[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
      }
    }
    legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="bottomleft", cex=thisCex)
    abline(v=grep("2016",QtimeYears), col="red", lty=4)
    mtext(paste(capText,", h=", htext, sep=""), side=3,  adj=0, cex=thisCex)

  }
}

dev.off()

## separate figures
for(r in 1:2){
  if(r==1){
    baseCappedIndex<-lookup_df$SSR[catchScenIndex]=="bin";  
    capText<-"No cap"
  } else{
    baseCappedIndex<-lookup_df$SSR[catchScenIndex]!="bin";  
    capText<-"Cap"
  }
  pdf(paste(plotPath, "KemptonsQhist_s",r,".pdf", sep=""), height=3, width=4.5)
  par(mar=c(4,4,1.5,0.5))
  plot(TLarray[1,histYearsIndex],type="n", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Kempton's Q")
  axis(at=axisHistYearsAt,labels = axisHistYears,side=1,las=2)
  # polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
  #         border=NA, col=myGreen_trans)
  for(h in 1:3){
    points(TLarray[baseCappedIndex,histYearsIndex][h,],type="l", col=colBySteepness[h])
  }
  legend(legend=steepnessSens, lty=1, seg.len=3, bty="n", col=colBySteepness, x="bottomright", inset=0.01, title="Steepness")
  # abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste(capText, sep=""), side=3,  adj=0)
  dev.off()
}

# 
# 
# #create array of which ones fall outside the boundaries (they are all below, but check for above too)
# # dimensions are SSR, h, scenario, timestep
# summaryYears<-c(2020, 2025, 2030, 2040); summaryTimeIndex<-summaryYears-1865 ; nsy<-length(summaryYears)
# summaryArray<-array(NA, dim=c(2,nh, nscenarios, nsy))
# for(s in 1:2){
#   if(s==1){
#     thisData<-Qtimes[lookup_df$SSR=="bin" ,summaryTimeIndex]
#   } else{
#     thisData<-Qtimes[lookup_df$SSR!="bin" ,summaryTimeIndex]
#   }
# }


#####################################################
# thisHlines<-pretty(seq(0,thisYmax,length.out=5))
# ## snapshot at 20, then 50 years
# thisColors<-colorRampPalette(colors=c(myGreen, myBlue,"midnightblue"))(nTestScen)
# yearsToPlot<-c(2016,2020,2030,2040,2046); timeIndexes<-yearsToPlot-1865
# for(thisTime in timeIndexes){
#   par(lend=1)
#   plot(TLarray[lookup_df$SSR=="bin" & lookup_df$scenario=="All0catch",thisTime], type="n", lwd=5,xlab="Steepness", ylab="Q"
#        ,xaxt="n", ylim=c(0,thisYmax), xlim=c(0.5,5.5))
#   abline(h=thisHlines, col=myGrey_trans, lwd=2)
#   axis(at=seq(1,nh), labels = as.character(steepnessSens), side=1)
#   for(s in 1:nTestScen){
#     thisXshift<-(s/nTestScen)*0.4 - 0.25
#     thisY<-TLarray[lookup_df$SSR=="bin" & lookup_df$scenario==testScenarios[s],thisTime]
#     points(x=seq(1,nh)+thisXshift, y=thisY, col=thisColors[s], type="h", lwd=5)
#   }
#   mtext(paste(thisTime," years", sep=""), side=3, adj=1)
#   mtext("No cap", side=3, adj=0)
#   plot(TLarray[lookup_df$SSR=="bin" & lookup_df$scenario=="All0catch",thisTime], type="n", lwd=5,xlab="Steepness", ylab="Q"
#        ,xaxt="n", ylim=c(0,thisYmax), xlim=c(0.5,5.5))
#   abline(h=thisHlines, col=myGrey_trans, lwd=2)
#   axis(at=seq(1,nh), labels = as.character(steepnessSens), side=1)
#   for(s in 1:nTestScen){
#     thisXshift<-(s/nTestScen)*0.4 - 0.25
#     thisY<-TLarray[lookup_df$SSR!="bin" & lookup_df$scenario==testScenarios[s],thisTime]
#     points(x=seq(1,nh)+thisXshift, y=thisY, col=thisColors[s], type="h", lwd=5)
#   }
#   mtext("Cap", side=3, adj=0)
# }
# 
# 
# ##################################################################
# ## now to plot as heat map
# getColor<-function(x,thisMax){
#   thisCol<-"red"
#   if(!is.na(x) & x>=thisMin){
#     y<-round((x-thisMin)/(thisMax-thisMin),2)*100+1
#     thisCol<-thisColRamp[y]
#   }
#   if(x>thisMax){thisCol<-thisColRamp[length(thisColRamp)]}
#   return(thisCol)
# }
# plotGrid<-function(x,y){
#   thisX<-c(x-0.5,x-0.5,x+0.5,x+0.5); thisY<-c(y-0.5,y+0.5,y+0.5,y-0.5)
#   thisCol<-plotColour[x,y]
#   if(length(thisCol)>0){
#     polygon(x=thisX,y=thisY,col=thisCol,border=NA)
#   }
#   return(NULL)
# }
# thisColRamp<-colorRampPalette(colors=c(myLightAqua, myAqua,"midnightblue"))(102)[-1]
# 
# 
# 
# Qtimes_df<-apply(data.frame(Qtimes[,QtimeIndex]),c(1,2), FUN=function(x){as.double(as.character(x))})
# plotData<-t(Qtimes_df)
# thisMax<-max(plotData, na.rm=TRUE); thisMin<-min(plotData, na.rm=TRUE)
# 
# for(p in 1:dim(plotColor)[2]){
#   plot(rep(1,dim(plotColor)[2]), col=plotColour[,10], pch=20, cex=1.5)
#   mtext(p,side=3)
# }
# 
# plotColour<-apply(plotData,c(1,2),getColor,thisMax)
# tempDF<-data.frame(cbind("x"=rep(seq(1,dim(plotData)[1]),dim(plotData)[2]),"y"=sort(rep(seq(1,dim(plotData)[2]),dim(plotData)[1]))))
# pdf(paste(plotPath,"KemptonsQ_SSRsens.pdf",sep=""),height=5,width=8)
# par(mar=c(4.5,9,1.5,1))
# plot(x=seq(1,dim(plotData)[1]),y=rep(dim(plotData)[2],dim(plotData)[1]),type="n",xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,dim(plotData)[2]))
# axis(at=axisYearsAt,labels = axisYears,side=1,las=2)
# axis(at=seq(1,dim(plotData)[2]),labels=steepnessSens[lookup_df$h],side=2,las=1)
# par(las=1)
# axis(at=grep("All100catch",lookup_df$scenario), labels=c(rep("100",(nh * nSSR))), side=2, line=2.5)
# axis(at=grep("All50catch",lookup_df$scenario), labels=c(rep("50",(nh * nSSR))), side=2, line=2.5)
# axis(at=grep("All0catch",lookup_df$scenario), labels=c(rep("0",(nh * nSSR))), side=2, line=2.5)
# axis(at=grep("All150catch",lookup_df$scenario), labels=c(rep("150",(nh * nSSR))), side=2, line=2.5)
# 
# axis(at=grep("R0",lookup_df$SSR, invert=TRUE), labels=c(rep("",(nruns/2 -1)), "\nNot\ncapped"), side=2, line=5)
# axis(at=grep("R0",lookup_df$SSR), labels=c(rep("",(nruns/2 -1)), "Capped"), side=2, line=5)
# temp<-Map(plotGrid,x=tempDF$x,y=tempDF$y)
# 
# mtext("h", side=3, line=-0.5, adj=-0.05, font=2)
# mtext("Catch", side=3, line=-0.5, adj=-0.2, font=2)
# mtext("SSR", side=3, line=-0.5, adj=-0.3, font=2)
# abline(v=snapshotAt, col="red", lty=2, lwd=2)
# box()
# dev.off()
# 
# 
# legendText<-seq(1.5,5.5, by=1); legendCol<-unlist(lapply(legendText, getColor, thisMax=thisMax))
# pdf(paste(plotPath,"KemptonsQ_SSRsensLEGEND.pdf",sep=""),height=4,width=2)
# par(mar=c(0,0,0,0))
# makeBlankPlot()
# par(xpd=TRUE)
# legend(legend=legendText, col=legendCol, pch=15, pt.cex=1.5, x="center", title="Q", bty="n")
# dev.off()
# ########################################################################
# 
# 
# 
