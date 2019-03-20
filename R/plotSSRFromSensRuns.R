source(paste(DIR$'General functions',"calcQ.r", sep=""))

Version<-"D"; versionText<-Version
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep=""); basePath<-thisPath
dataOutPath<-paste(thisPath,"SSRsens\\SSRsensTracers\\V",Version,sep="")

plotPath<-paste(DIR$'Base',"ATLANTISmodels\\Figures\\Recruitment\\V",Version,sep="")
# plotPath<-paste(DIR$'Figures',"Recruitment\\postReview\\",sep="") ##overwrites version in paper

groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

SSRdescriptions<-c("No cap", "Cap") ## check these match the SSR called out of lookup_df later. Used in plots

#first one created in populateTracersArrayFromSSRhModels.R - but could use alt scrip for different runs
load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracersArray", "storeNumbersArray", "storeWeightsArray", "storeYOYarray", "lookup_df
## . dim = dim=c(nruns, ntracers, nyears), c(nruns, ntracers,10, nyears), c(nruns, ntracers,10, nyears), c(nruns, ntracers, nyears) resp.
## lookup_df has definitions of runs including run colours

# # replaced by this version - holds catch scenarios too. For this, just take base catch scenario
# load(paste(thisPath,"SSRsens\\SSRsensTracers\\CatchScenarios_VDmodelTracers", sep=""))
# catchIndex <- lookup_df$scenario=="All100catch"
# lookup_df<-lookup_df[catchIndex,]; storeTracersArray<-storeTracersArray[catchIndex,,]

nruns<-dim(lookup_df)[1]
burnin<-110

steepnessSens<-seq(0.5,0.9,by=0.1)

thisBin<-lookup_df$SSR[1]; this_hindex<-lookup_df$h[1]
this_h<-steepnessSens[this_hindex]
thisScenario<-""; 
thisOut<-paste("outputSSR",thisBin,thisScenario,"h",this_hindex,sep="")
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\SSRsens",versionText,"\\",thisOut,"\\",sep="")

BaseFolder <- paste(DIR$'Base', "ATLANTISmodels\\ArchivedModels\\Base\\outputFish\\", sep="")
BaseNC.nc<-nc_open(paste(BaseFolder,"\\output.nc", sep=""))
thisVol<-ncvar_get(BaseNC.nc,"volume")
nts<-dim(thisVol)[3]
nyears<-nts; nlayers<-dim(thisVol)[1]

allTracers<-unique(names(BaseNC.nc$var))
xx<-grep("_N",allTracers); yy<-grep("_Nums",allTracers)
Ntracers<-allTracers[xx[!(xx %in% yy)]]; ntracers<-length(Ntracers)

timeAxisLabels<-seq(1975,2016,by=10); timeAxisAt<-seq(1,(nyears-burnin+1), by=10)

##steepness legend
pdf(paste(plotPath, "PFShLEGEND.pdf", sep=""), height=2, width=2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=seq(0.5,0.9,by=0.1), col=unique(lookup_df$hcolor), lwd=3, seg.len=3, x="center", bty="n", title="Steepness")
dev.off()

##steepness legend
pdf(paste(plotPath, "PFShLEGENDhoriz.pdf", sep=""), height=1, width=5)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), lend=1)
makeBlankPlot()
legend(legend=seq(0.5,0.9,by=0.1), col=unique(lookup_df$hcolor), lwd=5, seg.len=3, x="center", bty="n", title="Steepness", horiz=TRUE)
dev.off()

#########################################
# first plot the ssr curve for PFS for the runs - check they are actually different. If not, can't expect a change in the ecosystem
## also plot biomass tracer
get_KWRR<-function(Code){
  thisPar<-paste("KWRR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
get_KWSR<-function(Code){
  thisPar<-paste("KWSR_",Code,sep="")
  thisOut<-get_first_number(biolLines[grep(thisPar,biolLines)])
  return(thisOut)
}
baseBiolFile<-paste(basePath,"CRAM_BH_hybrid_biol.prm",sep=""); biolLines<-readLines(baseBiolFile)
base_h<-0.7
thisVar<-"BHalpha_PFS"; xx<-biolLines[grep(thisVar, biolLines)]; baseAlpha<-get_first_number(xx)
thisVar<-"BHbeta_PFS"; xx<-biolLines[grep(thisVar, biolLines)]; baseBeta<-get_first_number(xx)
ux<-5*base_h-1; thisB0<-(baseBeta*ux)/(1-base_h); thisR0<-(baseAlpha * ux)/ (4*base_h)

thisCode<-"PFS"; thisName<-str_trim(groupsDF$Name[groupsDF$Code==thisCode])
tracerIndex<-grep(thisName, Ntracers)    
timeIndex<-seq(burnin,nyears) #takes from 1970

test<-storeTracersArray[,tracerIndex,]

SSRs<-unique(lookup_df$SSR); nSSR<-length(SSRs); 

# par(mar=c(4.5,5,1,1), mfrow=c(2,1))
# for(s in 1:nSSR){
#   runIndex<-lookup_df$SSR==SSRs[s]
#   thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; thisYOYarray<-storeYOYarray[runIndex,tracerIndex,timeIndex]
#   this_nruns<-dim(thisSSBArray)[1]
#   thisYmax<-max(thisYOYarray, na.rm=TRUE); thisXmax<-max(thisSSBArray, na.rm=TRUE)
#   testSSB<-seq(0,thisXmax*1.1,length.out=1000)
#   # pdf(paste(plotPath,"SSRProjdemo.pdf",sep=""), height=4,width=5.5)
#   
#   plot(x=thisSSBArray[1,], y=thisYOYarray[1,], type="n", xlim=c(0,thisXmax),  ylim=c(0,thisYmax), xlab="SSB (tonnes)", ylab="Recruitment (numbers)");
#   mtext(SSRdescriptions[s], side=3, adj=0, font=2)
#   
#   # plot(x=thisSSBArray[1,], y=thisYOYarray[1,])
#   # abline(v=thisB0tonnes,lty=1,lwd=3,col=myGrey_trans); abline(h=thisR0, lty=1, lwd=3, col=myGrey_trans)
#   for(r in 1:this_nruns){
#     thisCol<-lookup_df$hcolor[r]; thisLty<-1
#     if(r<7){thisLty=2}
#     points(x=thisSSBArray[r,], y=thisYOYarray[r,],col=thisCol, pch=20)
#   }
# }
# 
# par(mar=c(4.5,5,1,1), mfrow=c(2,1))
# for(s in 1:nSSR){
#   runIndex<-lookup_df$SSR==SSRs[s]
#   thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; 
#   this_nruns<-dim(thisSSBArray)[1]
#   thisYmax<-max(thisYOYarray, na.rm=TRUE); thisXmax<-max(thisSSBArray, na.rm=TRUE)
#   testSSB<-seq(0,thisXmax*1.1,length.out=1000)
#   # pdf(paste(plotPath,"SSRProjdemo.pdf",sep=""), height=4,width=5.5)
#   
#   plot(thisSSBArray[1,], type="n", ylim=c(0,thisXmax), xaxt="n", ylab="SSB (tonnes)", xlab="Timestep (years)");
#   mtext(SSRdescriptions[s], side=3, adj=0, font=2)
#   axis(at=timeAxisAt, labels=timeAxisLabels, side=1)
#   
#   abline(v=thisB0,lty=1,lwd=3,col=myGrey_trans); abline(h=thisR0, lty=1, lwd=3, col=myGrey_trans)
#   for(r in 1:this_nruns){
#     thisCol<-lookup_df$hcolor[r]; thisLty<-1
#     if(r<7){thisLty=2}
#     points(thisSSBArray[r,],col=thisCol, lwd=2, type="l")
#   }
# }

thisYmax<-max(storeTracersArray[,tracerIndex,], na.rm=TRUE)*1.1

thisB0tonnes<-storeTracersArray[1,tracerIndex,1]/0.9

#############################
## in one plot - version in paper
###################################
pdf(paste(plotPath,"PFS_biomassHindcast.pdf", sep=""), height=4.5, width=12)
par(mar=c(4.5,5,1,3), fig=c(0,0.5,0.15,1), oma=c(0,0,0,0))
s=1
runIndex<-lookup_df$SSR==SSRs[s]
thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; 
this_nruns<-dim(thisSSBArray)[1]
plot(thisSSBArray[1,], type="n", xaxt="n", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="Timestep (years)");
axis(at=timeAxisAt, labels=timeAxisLabels, side=1)
mtext(paste("A: ", SSRdescriptions[s], sep=""), side=3, adj=0, font=2)
abline(h=thisB0tonnes,lty=3,lwd=3,col="red"); 
par(las=1)
axis(at=thisB0tonnes, labels = expression(B[0]), side=4, col="red", col.axis="red")
for(r in 1:this_nruns){
  thisCol<-lookup_df$hcolor[r]; thisLty<-1
  if(r<7){thisLty=2}
  points(thisSSBArray[r,],col=thisCol, lwd=2, type="l")
}
par(mar=c(4.5,5,1,3), fig=c(0.5,1,0.15,1), las=0, oma=c(0,0,0,0))
par(new=TRUE)
s=2
  runIndex<-lookup_df$SSR==SSRs[s]
  thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; 
  this_nruns<-dim(thisSSBArray)[1]

  plot(thisSSBArray[1,], type="n", xaxt="n", ylim=c(0,thisYmax), ylab="", xlab="Timestep (years)");
  mtext("SSB (tonnes)", side=2, line=3)
  axis(at=timeAxisAt, labels=timeAxisLabels, side=1)
  mtext(paste("B: ", SSRdescriptions[s], sep=""), side=3, adj=0, font=2)
  abline(h=thisB0tonnes,lty=3,lwd=3,col="red"); 
  par(las=1)
  axis(at=thisB0tonnes, labels = expression(B[0]), side=4, col="red", col.axis="red")
  for(r in 1:this_nruns){
    thisCol<-lookup_df$hcolor[r]; thisLty<-1
    if(r<7){thisLty=2}
    points(thisSSBArray[r,],col=thisCol, lwd=2, type="l")
  }
  # mtext("B", font=2, side=3, adj=0.5, outer=TRUE)
par(mar=c(0,0,0,0), fig=c(0,1,0,0.2), lend=1)
par(new=TRUE)
makeBlankPlot()
legend(legend=seq(0.5,0.9,by=0.1), col=unique(lookup_df$hcolor), lwd=5, seg.len=3, x="bottom", bty="n", title="Steepness", horiz=TRUE, inset=-0.1)

dev.off()

for(s in 1:nSSR){
  capDesc<-c("noCap","Cap")[s]
  pdf(paste(plotPath,"PFS_biomassHindcast_",capDesc,".pdf", sep=""), height=3, width=5)
  par(mar=c(4.5,5,1,3), mfrow=c(1,1))
  runIndex<-lookup_df$SSR==SSRs[s]
  thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; 
  this_nruns<-dim(thisSSBArray)[1]
  # thisYmax<-max(thisYOYarray, na.rm=TRUE); thisXmax<-max(thisSSBArray)
  testSSB<-seq(0,thisXmax*1.1,length.out=1000)

  plot(thisSSBArray[1,], type="n", xaxt="n", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="Timestep (years)");
  axis(at=timeAxisAt, labels=timeAxisLabels, side=1)
  mtext(SSRdescriptions[s], side=3, adj=0, font=2)
  # plot(x=thisSSBArray[1,], y=thisYOYarray[1,])
  abline(h=thisB0tonnes,lty=3,lwd=3,col="red"); 
  par(las=1)
  axis(at=thisB0tonnes, labels = expression(B[0]), side=4, col="red", col.axis="red")
  for(r in 1:this_nruns){
    thisCol<-lookup_df$hcolor[r]; thisLty<-1
    if(r<7){thisLty=2}
    points(thisSSBArray[r,],col=thisCol, lwd=2, type="l")
  }
  dev.off()
}
par(mfrow=c(2,2))
for(t in 1:ntracers){
  thisTracer <- Ntracers[t]; tracerIndex <- t
  thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; 
  this_nruns<-dim(thisSSBArray)[1]
  thisYmax<-max(thisSSBArray, na.rm=TRUE); 

  plot(thisSSBArray[1,], type="n", xaxt="n", ylim=c(0,thisYmax), ylab="SSB (tonnes)", xlab="Timestep (years)");
  axis(at=timeAxisAt, labels=timeAxisLabels, side=1)
  mtext(SSRdescriptions[s], side=3, adj=0, font=2)
  # plot(x=thisSSBArray[1,], y=thisYOYarray[1,])
  abline(h=thisB0tonnes,lty=3,lwd=3,col="red"); 
  par(las=1)
  axis(at=thisB0tonnes, labels = expression(B[0]), side=4, col="red", col.axis="red")
  for(r in 1:this_nruns){
    thisCol<-lookup_df$hcolor[r]; thisLty<-1
    if(r<7){thisLty=2}
    points(thisSSBArray[r,],col=thisCol, lwd=2, type="l")
  }
  mtext(thisTracer, side=3, adj=1)
}


# for(s in 1:nSSR){
runIndex<-lookup_df$h==1 #& lookup_df$SSR==SSRs[s]
thisSSBArray<-storeTracersArray[runIndex,tracerIndex,timeIndex]; 
this_nruns<-dim(thisSSBArray)[1]
thisYmax<-max(thisYOYarray, na.rm=TRUE); thisXmax<-max(thisSSBArray)
testSSB<-seq(0,thisXmax*1.1,length.out=1000)
# pdf(paste(plotPath,"SSRProjdemo.pdf",sep=""), height=4,width=5.5)

plot(thisSSBArray[1,], type="n", ylim=c(0,thisXmax), ylab="SSB (tonnes)", xlab="Timestep (years)");
# plot(x=thisSSBArray[1,], y=thisYOYarray[1,])
abline(v=thisB0,lty=1,lwd=3,col=myGrey_trans); abline(h=thisR0, lty=1, lwd=3, col=myGrey_trans)
for(r in 1:this_nruns){
  thisCol<-lookup_df$SSRcolor[runIndex][r]; thisLty<-1
  if(r<7){thisLty=2}
  points(thisSSBArray[r,],col=thisCol, lwd=2, type="l")
}
legend(legend=SSRdescriptions, col = lookup_df$SSRcolor[match(SSRs, lookup_df$SSR)], x="bottomright", lty=1,lwd=3, bty="n")
# dev.off()
# }

test_df <- data.frame(matrix(NA, ncol=5, nrow=10))
colnames(test_df)<-c("A", "B", "C", "D", "E")
test_df$A <- c(1,1,2,2,3,rep(4,5))
test_df$B <- rep(1,10)
