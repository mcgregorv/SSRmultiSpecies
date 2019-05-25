##tracks abundance, also converts to biomass, and calculates recruitment which is added to abundance with natural mortality then removed
basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(DIR$'Base',"\\Figures\\Recruitment\\postReview\\",sep="")

# recModelText<-"baseBiomass" #this goes into the plot file name
recModelText<-"simpleSwitchBiomass" #this goes into the plot file name
thisCex<-1.4; thisAxisCex<-1

##can replace this with new recruitment model and rerun this code
calc_rec<-function(B,a,b,B0,new_a){
  if(B<B0){
    R<-(a*B)/(b+B)
  }else{
    R<-(new_a*B)/(b+B)
  }
  return(R)
}

calc_a<-function(h,R0){
  a<-(4*h*R0)/(5*h-1)
  return(a)
}
calc_b<-function(h,B0){
  b<-(B0*(1-h))/(5*h-1)
  return(b)
}

R0<-800; B0<-600

aveWeight<-5e-5 #this is average individual weight in tonnes (5e-5 tonnes = 0.5 kg = 500 grams)

N0<-(B0/aveWeight)/1000 #numbers are in 1000's 

baseM<-(-1)*log(N0/(N0+R0),base=exp(1))

sinFn <- function(x, baseM){
  y <- baseM*0.8 * sin((x * pi)/100) + baseM
  return(y)
}
test <- seq(0,200, length.out=1000); testy <- unlist(lapply(test, sinFn, baseM=baseM))

plot(x=test, y=testy, type="l")


Mscalar=testy+rnorm(length(testy),0,0.005)

pdf(paste(plotPath, "CompareNewSimM.pdf", sep=""), height=4, width=7)
par(las=1, mar=c(4,4,1,4))
plot(x=test, y=Mscalar, type="l", col=myGrey, lwd=1.2, ylim=c(min(Mscalar), max(Mscalar)), ylab="M", xlab="Timestep (years)", xlim=c(0, 235))
points(x=test, y=testy, col="black", lty=2, lwd=2, type="l")
abline(h=baseM, col="red", lty=2, lwd=2)
axis(at=baseM, labels = "Base M", col.axis="red", side=4, col.ticks = "red")

#write M out so can be replicated exactly if needed
# write.csv(MbyTime, paste(plotPath,"Mbytime.csv", sep=""), row.names = FALSE)
# read in the old one to compare
MbyTime<- (read.csv(paste(plotPath,"Mbytime.csv", sep="")))[,1]
# par(new=TRUE)
points(x=0:(length(MbyTime)-1),y=MbyTime, type="l", col=myBlue, xlab="", ylab="")
dev.off()

#write M out so can be replicated exactly if needed
MbyTime_df <- data.frame(cbind("Timestep"=test, "M"=Mscalar, "SmoothM"=testy))
# write.csv(MbyTime, paste(plotPath,"MbytimeSimple.csv", sep=""), row.names = FALSE)

MbyTime_df<- (read.csv(paste(plotPath,"MbytimeSimple.csv", sep="")))

pdf(paste(plotPath,"MbyTime.pdf",sep=""),height=4, width=7)
par(mfrow=c(1,1),mar=c(4,4.5,2,0.5),oma=c(0,0,2,0))
plot(x=MbyTime_df$Timestep, y=MbyTime_df$M, type="l", col=myGrey,ylab="M", xlab="Timestep (years)", cex.axis=thisCex, cex.lab=thisCex)
points(x=MbyTime_df$Timestep, y=MbyTime_df$SmoothM, type="l", lty=2, lwd=2)
abline(h=baseM,col="red",lty=2, lwd=1.5)
dev.off()

this_h<-0.8; htext<-gsub("\\.","",as.character(this_h))

##do the cap population, which is the result of the lowest M, and if we are using the original BH, the lowest h
startPop<-N0; basePopulation<-startPop; cap_h<-min(this_h)

MbyTime <- MbyTime_df$M
nts<-length(MbyTime)

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
capN<-N0*10
test_h<-seq(0.35,0.95,by=0.05); nhs<-length(test_h)

colbyh<-colorRampPalette(colors=c(myAqua, myBlue,"midnightblue"))(nhs)
hTextIndex<-c(1,round(nhs/2),nhs); test_h_text<-rep("",nhs); test_h_text[hTextIndex]<-test_h[hTextIndex]

pdf(paste(plotPath,"h_LEGEND.pdf",sep=""),height=6.5,width=1.5)
par(mar=c(0,0,0,0))
makeBlankPlot()
legend(legend=test_h_text,col=colbyh,lty=NA,lwd=NA, pch=15,bty="n",x="center",title="h",cex=1.8, pt.cex=4, ncol=1)
dev.off()


recModelText<-"baseBiomass"
# recModelText<-"capR"

SSBarray<-array(NA, dim=c(nhs,nts)); recArray<-0*SSBarray
for(h in 1:nhs){
  this_h<-test_h[h]
  htext<-gsub("\\.","",as.character(this_h))
  
  a<-calc_a(h=this_h, R0=R0); b<-calc_b(h=this_h,B0=B0)
  # if(recModelText=="baseBiomass" ){
  new_a<-a
  # }else{
  #   new_a<-(1+this_xi*this_h)*R0
  # }
  
  basePopulation<-startPop
  recRecruitment<-calc_rec(B=B0,a=a,b=b,B0=B0,new_a=new_a)
  nts<-length(Mscalar)
  for(t in 2:nts){
    thisSSB<-basePopulation[t-1]*aveWeight*1000
    thisRec<-recRecruitment[t-1]
    if(recModelText=="capR" ){
      thisRec<-min(R0, thisRec)
    }
    # thisRec<-calc_rec(B=basePopulation[t-1],a=a,b=b,B0=B0,new_a=new_a)
    newPop<-(basePopulation[t-1]+thisRec)*exp(-MbyTime[t])
    basePopulation<-c(basePopulation,newPop)
    newSSB<-basePopulation[t]*aveWeight*1000
    newRec<-calc_rec(B=newSSB,a=a,b=b,B0=B0,new_a=new_a)
    if(recModelText=="capR" ){
      newRec<-min(R0, newRec)
    }
    recRecruitment<-c(recRecruitment,newRec)
  }
  SSBarray[h,]<-basePopulation*aveWeight*1000; recArray[h,]<-recRecruitment
}

## all in one figure version
pdf(paste(plotPath,"SimplePopulationSSRandBiomass.pdf",sep=""),height=5,width=15)
par(mfrow=c(1,2), oma=c(1,1,2,10), mar=c(4,4,1,1))

yMax<-max(recArray, na.rm=TRUE); xMax<-max(SSBarray, na.rm=TRUE)
yMax<-1200; xMax<-1800
plot(x=SSBarray[1,], y=recArray[1,], ylim=c(0, yMax), xlim=c(0,xMax), type="n", ylab="Recruitment", xlab="SSB", cex.axis=thisCex, cex.lab=thisCex)
for(h in 1:nhs){
  points(x=SSBarray[h,], y=recArray[h,], type="l", lwd=2, col=colbyh[h])
}
abline(h=R0, col="red", lty=2)
abline(v=B0, col="red", lty=2)
mtext("A", font=2, side=3, adj=0, line=0, outer=TRUE, cex=thisCex)
#
yMax<-max(SSBarray, na.rm=TRUE)
yMax<-4050
plot(SSBarray[1,], ylim=c(0, yMax),  type="n", ylab="SSB (tonnes)", xlab="Timestep (years)", cex.axis=thisCex, cex.lab=thisCex)
for(h in 1:nhs){
  points(SSBarray[h,], type="l", col=colbyh[h], lwd=2)
}
abline(h=B0, col="red", lty=2)
mtext("B", font=2, side=3, adj=0.5, line=0, outer=TRUE, cex=thisCex)
par(xpd=NA)
legend(legend=test_h_text,col=colbyh,lty=NA,lwd=NA, pch=15,bty="n",x="right",title="h",cex=1.5, pt.cex=3.5, ncol=1,inset=-0.3)
dev.off()



#plot SSR relationship
pdf(paste(plotPath,"SRRSimplePopulation",recModelText,".pdf",sep=""),height=4, width=5)
par(mfrow=c(1,1),mar=c(4,4.5,2,0.5),oma=c(0,0,2,0))
yMax<-max(recArray, na.rm=TRUE); xMax<-max(SSBarray, na.rm=TRUE)
yMax<-1200; xMax<-1800
plot(x=SSBarray[1,], y=recArray[1,], ylim=c(0, yMax), xlim=c(0,xMax), type="n", ylab="Recruitment", xlab="SSB", cex.axis=thisCex, cex.lab=thisCex)
for(h in 1:nhs){
  points(x=SSBarray[h,], y=recArray[h,], type="l", lwd=2, col=colbyh[h])
}
abline(h=R0, col="red", lty=2)
abline(v=B0, col="red", lty=2)
dev.off()

#plot SSR relationship
pdf(paste(plotPath,"SSBSimplePopulation",recModelText,".pdf",sep=""),height=4, width=5)
par(mfrow=c(1,1),mar=c(4,4.5,2,0.5),oma=c(0,0,2,0))
yMax<-max(SSBarray, na.rm=TRUE)
yMax<-1800
plot(SSBarray[1,], ylim=c(0, yMax),  type="n", ylab="SSB (tonnes)", xlab="Timestep (years)", cex.axis=thisCex, cex.lab=thisCex)
for(h in 1:nhs){
  points(SSBarray[h,], type="l", col=colbyh[h], lwd=2)
}
abline(h=B0, col="red", lty=2)
dev.off()
# 
# #plot SSR relationship
# pdf(paste(plotPath,"SSBSimplePopulation",recModelText,"plain.pdf",sep=""),height=4, width=5)
# par(mfrow=c(1,1),mar=c(4,4.5,2,0.5),oma=c(0,0,2,0))
# yMax<-max(SSBarray, na.rm=TRUE)
# plot(SSBarray[1,], ylim=c(0, yMax),  type="n", ylab="SSB (tonnes)", xlab="Timestep (years)", cex.axis=thisCex, cex.lab=thisCex)
# for(h in 1:nhs){
#   points(SSBarray[h,], type="l", col=colbyh[h], lwd=2)
# }
# abline(h=B0, col="red", lty=2)
# dev.off()
# 

###############################################################################################
# CAPPED VERSION ###

recModelText<-"capR"

SSBarray<-array(NA, dim=c(nhs,nts)); recArray<-0*SSBarray
for(h in 1:nhs){
  this_h<-test_h[h]
  htext<-gsub("\\.","",as.character(this_h))
  
  a<-calc_a(h=this_h, R0=R0); b<-calc_b(h=this_h,B0=B0)
  # if(recModelText=="baseBiomass" ){
  new_a<-a
  # }else{
  #   new_a<-(1+this_xi*this_h)*R0
  # }
  
  basePopulation<-startPop
  recRecruitment<-calc_rec(B=B0,a=a,b=b,B0=B0,new_a=new_a)
  nts<-length(Mscalar)
  for(t in 2:nts){
    thisSSB<-basePopulation[t-1]*aveWeight*1000
    thisRec<-recRecruitment[t-1]
    if(recModelText=="capR" ){
      thisRec<-min(R0, thisRec)
    }
    # thisRec<-calc_rec(B=basePopulation[t-1],a=a,b=b,B0=B0,new_a=new_a)
    newPop<-(basePopulation[t-1]+thisRec)*exp(-MbyTime[t])
    basePopulation<-c(basePopulation,newPop)
    newSSB<-basePopulation[t]*aveWeight*1000
    newRec<-calc_rec(B=newSSB,a=a,b=b,B0=B0,new_a=new_a)
    if(recModelText=="capR" ){
      newRec<-min(R0, newRec)
    }
    recRecruitment<-c(recRecruitment,newRec)
  }
  SSBarray[h,]<-basePopulation*aveWeight*1000; recArray[h,]<-recRecruitment
}

## all in one figure version
pdf(paste(plotPath,"SimplePopulationCappedSSRandBiomass.pdf",sep=""),height=5,width=15)
par(mfrow=c(1,2), oma=c(1,1,2,10), mar=c(4,4,1,1))

yMax<-max(recArray, na.rm=TRUE); xMax<-max(SSBarray, na.rm=TRUE)
yMax<-1200; xMax<-1800
plot(x=SSBarray[1,], y=recArray[1,], ylim=c(0, yMax), xlim=c(0,xMax), type="n", ylab="Recruitment", xlab="SSB", cex.axis=thisCex, cex.lab=thisCex)
for(h in 1:nhs){
  points(x=SSBarray[h,], y=recArray[h,], type="l", lwd=2, col=colbyh[h])
}
abline(h=R0, col="red", lty=2)
abline(v=B0, col="red", lty=2)
mtext("A", font=2, side=3, adj=0, line=0, outer=TRUE, cex=thisCex)
#
yMax<-max(SSBarray, na.rm=TRUE)
yMax<-4050
plot(SSBarray[1,], ylim=c(0, yMax),  type="n", ylab="SSB (tonnes)", xlab="Timestep (years)", cex.axis=thisCex, cex.lab=thisCex)
for(h in 1:nhs){
  points(SSBarray[h,], type="l", col=colbyh[h], lwd=2)
}
abline(h=B0, col="red", lty=2)
mtext("B", font=2, side=3, adj=0.5, line=0, outer=TRUE, cex=thisCex)
par(xpd=NA)
legend(legend=test_h_text,col=colbyh,lty=NA,lwd=NA, pch=15,bty="n",x="right",title="h",cex=1.5, pt.cex=3.5, ncol=1,inset=-0.3)
dev.off()



