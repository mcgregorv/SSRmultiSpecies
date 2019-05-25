basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
plotPath<-paste(basePath,"\\Figures\\Recruitment\\SSB\\",sep="")

thisCex<-1.8; thisAxisCex<-1.6


source(paste(DIR$'General functions',"makeBlankPlot.R",sep=""))

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

testS<-seq(0,10000,length.out = 1000)

R0<-800; B0<-600
h=0.25
baseBeta<-calc_b(B0=B0,h=h)
baseAlpha<-calc_a(h=h,R0=R0)

baseR<-unlist(lapply(testS,calc_rec,a=baseAlpha,b=baseBeta,B0=B0,new_a=baseAlpha))


thisCex<-1.3; thisAxisCex<-1.2
color2<-myRed
##steepness at 0.4
test_h<-seq(0.25,1,by=0.05); nh<-length(test_h)
col_h<-colorRampPalette(colors=c(myGold,myDarkGreen,"midnightblue"))(nh)
pdf(paste(plotPath,"BHExample_varySteepness.pdf",sep=""),height=5)
par(mar=c(4.5,5,5,6))
plot(x=testS,y=baseR,type="n",lwd=2,xlab="Spawning stock biomass",ylab="Recruitment",cex.axis=thisCex,cex.lab=thisCex)
for(i in 1:nh){
  h<-test_h[i]
  thisBeta<-calc_b(B0=B0,h=h)
  thisAlpha<-calc_a(h=h,R0=R0)
  thisR<-unlist(lapply(testS,calc_rec,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha))
  points(x=testS,y=thisR,type="l",lwd=2,col=col_h[i])
}
axis(at=B0,labels = expression('B'[0]),side=3,col.axis=color2,col=color2,cex.axis=thisAxisCex)
axis(at=R0,labels = expression('R'[0]),side=4,col.axis=color2,col=color2,cex.axis=thisAxisCex,las=2)
abline(v=B0,col=color2,lty=4,lwd=3.5)
abline(h=R0,col=color2,lty=4,lwd=3.5)
dev.off()

##lower part of the curve

testS<-seq(0,1000,length.out = 1000)
baseR<-unlist(lapply(testS,calc_rec,a=baseAlpha,b=baseBeta,B0=B0,new_a=baseAlpha))


pdf(paste(plotPath,"BHExample_varySteepnessLOWERcurve.pdf",sep=""),height=5)
par(mar=c(4.5,5,5,6))
plot(x=testS,y=baseR,type="n",lwd=2,xlab="Spawning stock biomass",ylab="Recruitment",cex.axis=thisCex,cex.lab=thisCex)
for(i in 1:nh){
  h<-test_h[i]
  thisBeta<-calc_b(B0=B0,h=h)
  thisAlpha<-calc_a(h=h,R0=R0)
  thisR<-unlist(lapply(testS,calc_rec,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha))
  points(x=testS,y=thisR,type="l",lwd=2,col=col_h[i])
}
axis(at=B0,labels = expression('B'[0]),side=3,col.axis=color2,col=color2,cex.axis=thisAxisCex)
axis(at=R0,labels = expression('R'[0]),side=4,col.axis=color2,col=color2,cex.axis=thisAxisCex,las=2)
abline(v=B0,col=color2,lty=4,lwd=3.5)
abline(h=R0,col=color2,lty=4,lwd=3.5)
dev.off()

###legend

keepLegend_h<-as.double(c(0.25,0.65,1)); legendIndex<-as.double(test_h)%in%keepLegend_h
legend_h_text<-test_h; legend_h_text[!legendIndex]<-""
pdf(paste(plotPath,"BHExample_varySteepnessLEGEND.pdf",sep=""),height=6,width=2)
par(mar=c(0,0,0,0),lend=1)
makeBlankPlot()
legend(legend=legend_h_text,col=col_h,x="center",lwd=4,seg.len=4,bty="n",title="Steepness",cex=thisCex)
dev.off()


##seepness example plots - plot for paper
thisAxisCex<-1.6; thisCex<-1.8; thisYlim<-950
pdf(paste(plotPath,"BHExample_h.pdf",sep=""),height=5, width=14)
par(mar=c(4.5,5,5,6), mfrow=c(1,2))
h<-0.4; hcol<-myBlue
h_text<-gsub("\\.","",h)
thisBeta<-calc_b(B0=B0,h=h)
thisAlpha<-calc_a(h=h,R0=R0)
thisR<-unlist(lapply(testS,calc_rec,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha))
R02<-calc_rec(0.2*B0,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha)

plot(x=testS,y=thisR,type="n",lwd=2,xlab="Spawning stock biomass",ylab="Recruitment",cex.axis=thisCex,cex.lab=thisCex, ylim=c(0, thisYlim))
axis(at=B0,labels = expression('B'[0]),side=3,col.axis="black",cex.axis=thisAxisCex,line=-0.2)
axis(at=R0,labels = expression('R'[0]),side=4,col.axis="black",cex.axis=thisAxisCex,las=2)
abline(v=B0,col="black",lty=4,lwd=3.5)
abline(h=R0,col="black",lty=4,lwd=3.5)
abline(v=0.2*B0,col=myGrey,lty=4,lwd=3.5)
abline(h=R02,col=myGrey,lty=4,lwd=3.5)
axis(at=0.2*B0,labels = expression(paste(0.2, 'B'[0])),side=3,col.axis=myGrey,cex.axis=thisAxisCex,line=1,col=myGrey,tck=0.1)
axis(at=R02,labels = expression('R'[0.2]),side=4,col.axis=myGrey,col=myGrey,cex.axis=thisAxisCex,las=2)
points(x=testS,y=thisR,type="l",lwd=3,col=hcol)
mtext(paste("h=",h,sep=""),side=3,adj=1,cex=thisCex)
mtext ("A", side=3, adj=0, line=-2, font=2, outer=TRUE, cex=thisCex)
## 
h<-0.8; hcol<-myDarkGreen
h_text<-gsub("\\.","",h)
thisBeta<-calc_b(B0=B0,h=h)
thisAlpha<-calc_a(h=h,R0=R0)
thisR<-unlist(lapply(testS,calc_rec,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha))
R02<-calc_rec(0.2*B0,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha)
plot(x=testS,y=thisR,type="n",lwd=2,xlab="Spawning stock biomass",ylab="Recruitment",cex.axis=thisCex,cex.lab=thisCex, ylim=c(0, thisYlim))
axis(at=B0,labels = expression('B'[0]),side=3,col.axis="black",cex.axis=thisAxisCex,line=-0.2)
axis(at=R0,labels = expression('R'[0]),side=4,col.axis="black",cex.axis=thisAxisCex,las=2)
abline(v=B0,col="black",lty=4,lwd=3.5)
abline(h=R0,col="black",lty=4,lwd=3.5)
abline(v=0.2*B0,col=myGrey,lty=4,lwd=3.5)
abline(h=R02,col=myGrey,lty=4,lwd=3.5)
axis(at=0.2*B0,labels = expression(paste(0.2, 'B'[0])),side=3,col.axis=myGrey,cex.axis=thisAxisCex,line=1,col=myGrey,tck=0.1)
axis(at=R02,labels = expression('R'[0.2]),side=4,col.axis=myGrey,col=myGrey,cex.axis=thisAxisCex,las=2)
points(x=testS,y=thisR,type="l",lwd=3,col=hcol)
mtext(paste("h=",h,sep=""),side=3,adj=1,cex=thisCex)
mtext ("B", side=3, adj=0.5, line=-2, font=2, outer=TRUE, cex=thisCex)
dev.off()

##########################
## zoomin
# h<-0.4; hcol<-myBlue
h<-0.8; hcol<-myDarkGreen
h_text<-gsub("\\.","",h)
thisBeta<-calc_b(B0=B0,h=h)
thisAlpha<-calc_a(h=h,R0=R0)
thisR<-unlist(lapply(testS,calc_rec,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha))
R02<-calc_rec(0.2*B0,a=thisAlpha,b=thisBeta,B0=B0,new_a=thisAlpha)

pdf(paste(plotPath,"BHExample_h",h_text,"zoomin.pdf",sep=""),height=5)
par(mar=c(4.5,5,5,6))
plot(x=testS,y=thisR,type="n",lwd=2,xlab="Spawning stock biomass",ylab="Recruitment",cex.axis=thisCex,cex.lab=thisCex,ylim=c(0,1.1*R0), xlim=c(0,1.1*B0))

axis(at=B0,labels = expression('B'[0]),side=3,col.axis="black",cex.axis=thisAxisCex,line=-0.2)
axis(at=R0,labels = expression('R'[0]),side=4,col.axis="black",cex.axis=thisAxisCex,las=2)
abline(v=B0,col="black",lty=4,lwd=3.5)
abline(h=R0,col="black",lty=4,lwd=3.5)

abline(v=0.2*B0,col=myGrey,lty=4,lwd=3.5)
abline(h=R02,col=myGrey,lty=4,lwd=3.5)
axis(at=0.2*B0,labels = expression(paste(0.2, 'B'[0])),side=3,col.axis=myGrey,cex.axis=thisAxisCex,line=1,col=myGrey,tck=0.1)
axis(at=R02,labels = expression('R'[0.2]),side=4,col.axis=myGrey,col=myGrey,cex.axis=thisAxisCex,las=2)

points(x=testS,y=thisR,type="l",lwd=3,col=hcol)
mtext(paste("h=",h,sep=""),side=3,adj=0,cex=thisCex)
# points(x=testS,y=thisR,type="l",lwd=3,col=myRed)
dev.off()

