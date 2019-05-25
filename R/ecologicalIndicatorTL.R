#ratio of pelagic biomass over total biomass 
Version<-"C"
outFolder<-"SSRsens"
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep=""); basePath<-thisPath
dataOutPath<-paste(thisPath,"SSRsens\\SSRsensTracers\\CatchScenarios_V",Version,sep="")
plotPath<-paste(DIR$'Base',"ATLANTISmodels\\Figures\\Recruitment\\Paper2\\V",Version,sep="")
# plotPath<-paste(DIR$'Figures',"Recruitment\\postReview\\", sep="") ## write up
colBySteepness<-colorRampPalette(colors=c(myGold, myGreen,myDarkGreen, "black"))(3)  
steepnessSens<-c(0.5,0.7,0.9)

# thisOutList<-c("TLratio", "Qtimes", "TLtimes_df", "meanLengthByRun","biomOverCatch", "ratioPelBiom", "sensZG")
# save(list=thisOutList,file=paste(dataOutPath,"ecoIndicators",sep=""))
load(paste(dataOutPath,"ecoIndicators", sep=""))
allScenarios<-unique(lookup_df$Scenario); nscenarios<-length(allScenarios)
scenarioColors<-colorRampPalette(colors=c(myOrange,"red",myPurple,myBlue))(nscenarios)

TLarray<-t(TLtimes_df); thisYmax<-max(TLarray, na.rm=TRUE)
thisYmax<-max(TLarray,na.rm=TRUE)*1.1

QtimeYears<-seq(2010,2046); QtimeIndex<-QtimeYears-1865 
axisYears<-QtimeYears[seq(1,length(QtimeYears),by=5)]; axisYearsAt<-seq(1,length(QtimeYears), by=5)
snapshotYears<-c(2020,2025,2035,2045); snapshotAt<-match(snapshotYears, QtimeYears)

histYears<-seq(1970,2016); histYearsIndex<-histYears-1865
axisHistYears<-histYears[seq(1,length(histYears),by=5)]; axisHistYearsAt<-seq(1,length(histYears), by=5)

thisMin<-min(TLarray, na.rm=TRUE); thisMax<-max(TLarray,na.rm=TRUE)
# thisMin<-21
preFishingMin<-min(TLarray[,1:min(QtimeIndex)], na.rm=TRUE); preFishingMax<-max(TLarray[,1:min(QtimeIndex)], na.rm=TRUE)
for(r in 1:2){
  for(h in 1:3){
    if(r==1){
      baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"No cap"
    } else{
      baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"Cap"
    }
    pdf(paste(plotPath, "MeanTLScenario_h",h,"_s",r,".pdf", sep=""), height=3, width=4.5)
    par(mar=c(4,4,1.5,0.5))
    plot(TLarray[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Mean trophic level")
    axis(at=axisYearsAt,labels = axisYears,side=1,las=2)
    polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
            border=NA, col=myGreen_trans)
    for(s in 1:nscenarios){
      if(s!=3){
        points(TLarray[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
      }
    }
    legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="topright")
    abline(v=grep("2016",QtimeYears), col="red", lty=4)
    mtext(paste(capText,", h=", htext, sep=""), side=3,  adj=0)
    dev.off()
  }
}

## do historic period comparing h for each SSR
for(r in 1:2){
  if(r==1){
    baseCappedIndex<-lookup_df$SSR=="bin"
    capText<-"No cap"
  } else{
    baseCappedIndex<-lookup_df$SSR!="bin";  
    capText<-"Cap"
  }
  pdf(paste(plotPath, "MeanTLhist_s",r,".pdf", sep=""), height=3, width=4.5)
  par(mar=c(4,4,1.5,0.5))
  plot(TLarray[1,histYearsIndex],type="n", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Mean trophic level")
  axis(at=axisHistYearsAt,labels = axisHistYears,side=1,las=2)
  # polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
  #         border=NA, col=myGreen_trans)
  for(h in 1:3){
    points(TLarray[baseCappedIndex,histYearsIndex][h,],type="l", col=colBySteepness[h])
  }
  legend(legend=steepnessSens, lty=1, seg.len=3, bty="n", col=colBySteepness, x="bottomleft", inset=0.01, title="Steepness")
  # abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste(capText, sep=""), side=3,  adj=0)
  dev.off()
}

################## COMBINE in one figure (version in paper)

## do historic period comparing h for each SSR
pdf(paste(plotPath, "MeanTL.pdf", sep=""), height=9, width=15)
par(mar=c(4,4.5,1.5,0.5), mfrow=c(2,2))

for(r in 1:2){
  if(r==1){
    baseCappedIndex<-lookup_df$SSR=="bin"
    capText<-"No cap"; thisLetter<-"A"
  } else{
    baseCappedIndex<-lookup_df$SSR!="bin";  
    capText<-"Cap"; thisLetter="B"
  }
  plot(TLarray[1,histYearsIndex],type="n", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Mean trophic level", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  axis(at=axisHistYearsAt,labels = axisHistYears,side=1,las=2, cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  # polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
  #         border=NA, col=myGreen_trans)
  for(h in 1:3){
    points(TLarray[baseCappedIndex,histYearsIndex][h,],type="l", col=colBySteepness[h])
  }
  legend(legend=steepnessSens, lty=1, seg.len=3, bty="n", col=colBySteepness, x="bottomleft", inset=0.01, title="Steepness", cex=thisCex)
  # abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste(thisLetter,": ",capText, sep=""), side=3,  adj=0, font=2, cex=thisCex)

}
for(r in 1:2){
  h=3
    if(r==1){
      baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"No cap"; thisLetter <-"C"
    } else{
      baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
      capText<-"Cap"; thisLetter <- "D"
    }
    plot(TLarray[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Mean trophic level", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
    axis(at=axisYearsAt,labels = axisYears,side=1,las=2, cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
    polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
            border=NA, col=myGrey_trans)
    for(s in 1:nscenarios){
      if(s!=3){
        points(TLarray[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
      }
    }
    legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="topright", cex=thisCex)
    abline(v=grep("2016",QtimeYears), col="red", lty=4)
    mtext(paste(thisLetter,": ",capText,", h=", htext, sep=""), side=3,  adj=0, font=2, cex=thisCex)

}
dev.off()

####################################################################
################## COMBINE for APPENDIX in paper #####################
#######################################################################
## do historic period comparing h for each SSR
pdf(paste(plotPath, "MeanTL_appendix.pdf", sep=""), height=15, width=13)
par(mar=c(4.5,4.5,2,1), mfrow=c(4,2))

for(r in 1:2){
  if(r==1){
    baseCappedIndex<-lookup_df$SSR=="bin"
    capText<-"No cap"; 
  } else{
    baseCappedIndex<-lookup_df$SSR!="bin";  
    capText<-"Cap"; 
  }
  plot(TLarray[1,histYearsIndex],type="n", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Mean trophic level", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  axis(at=axisHistYearsAt,labels = axisHistYears,side=1,las=2, cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  # polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
  #         border=NA, col=myGreen_trans)
  for(h in 1:3){
    points(TLarray[baseCappedIndex,histYearsIndex][h,],type="l", col=colBySteepness[h])
  }
  legend(legend=steepnessSens, lty=1, seg.len=3, bty="n", col=colBySteepness, x="bottomleft", inset=0.01, title="Steepness", cex=thisCex)
  # abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste(capText, sep=""), side=3,  adj=0, font=1, cex=thisCex)
  
}

  for(h in 1:3){
    for(r in 1:2){
  if(r==1){
    baseCappedIndex<-lookup_df$SSR=="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
    capText<-"No cap"; 
  } else{
    baseCappedIndex<-lookup_df$SSR!="bin" & lookup_df$h==h;  htext<-steepnessSens[h]
    capText<-"Cap"; 
  }
  plot(TLarray[baseCappedIndex,QtimeIndex][1,],type="l", ylim=c(thisMin, thisMax), xaxt="n", xlab="", ylab="Mean trophic level", cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  axis(at=axisYearsAt,labels = axisYears,side=1,las=2, cex=thisCex, cex.axis=thisCex, cex.lab=thisCex)
  polygon(x=c(seq(1,length(QtimeIndex)), rev(seq(1,length(QtimeIndex)))), y=c(rep(preFishingMin,length(QtimeIndex)), rep(preFishingMax, length(QtimeIndex))), 
          border=NA, col=myGrey_trans)
  for(s in 1:nscenarios){
    if(s!=3){
      points(TLarray[baseCappedIndex,QtimeIndex][s,],type="l", col=scenarioColors[s])
    }
  }
  legend(legend=c("Zero catch", "Status quo", "Half catch"), lty=1, seg.len=3, bty="n", col=scenarioColors[c(1:2,4)], x="topright", cex=thisCex)
  abline(v=grep("2016",QtimeYears), col="red", lty=4)
  mtext(paste(capText,", h=", htext, sep=""), side=3,  adj=0, font=1, cex=thisCex)
  }
}
dev.off()

# ##################################################
# thisHlines<-pretty(seq(thisYmin,thisYmax,length.out=5))
# 
# ## snapshot at 20, then 50 years
# thisColors<-colorRampPalette(colors=c(myGreen, myBlue,"midnightblue"))(nTestScen)
# thisTime<-20
# for(thisTime in seq(10,50,by=10)){
# par(lend=1)
# plot(TLarray[lookup_df$SSR=="bin" & lookup_df$scenario=="All0catch",thisTime], type="n", lwd=5,xlab="Steepness", ylab="Mean trophic level"
#      ,xaxt="n", ylim=c(thisYmin, thisYmax), xlim=c(0.5,5.5))
# abline(h=thisHlines, col=myGrey_trans, lwd=2)
# axis(at=seq(1,nh), labels = as.character(steepnessSens), side=1)
# for(s in 1:nTestScen){
#   thisXshift<-(s/nTestScen)*0.4 - 0.25
#   thisY<-TLarray[lookup_df$SSR=="bin" & lookup_df$scenario==testScenarios[s],thisTime]
#   points(x=seq(1,nh)+thisXshift, y=thisY, col=thisColors[s], type="h", lwd=5)
# }
# mtext(paste(thisTime," years", sep=""), side=3, adj=1)
# mtext("No cap", side=3, adj=0)
# plot(TLarray[lookup_df$SSR=="bin" & lookup_df$scenario=="All0catch",thisTime], type="n", lwd=5,xlab="Steepness", ylab="Mean Trophic level"
#      ,xaxt="n", ylim=c(thisYmin, thisYmax), xlim=c(0.5,5.5))
# abline(h=thisHlines, col=myGrey_trans, lwd=2)
# axis(at=seq(1,nh), labels = as.character(steepnessSens), side=1)
# for(s in 1:nTestScen){
#   thisXshift<-(s/nTestScen)*0.4 - 0.25
#   thisY<-TLarray[lookup_df$SSR!="bin" & lookup_df$scenario==testScenarios[s],thisTime]
#   points(x=seq(1,nh)+thisXshift, y=thisY, col=thisColors[s], type="h", lwd=5)
# }
# mtext("Cap", side=3, adj=0)
# }
# 
# 
# 
# 
