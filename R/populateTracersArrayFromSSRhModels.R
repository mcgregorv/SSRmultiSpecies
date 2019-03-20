## reads in model outputs from sens runs on R cap and steepness
## populates arrays with N tracers, numbers and weights that are stored so can be used for calculating diversity indices

Version<-"B"; versionText<-""
Version<-"C"; versionText<-"C"
Version<-"D"; versionText<-"D"
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\",sep=""); basePath<-thisPath
dataOutPath<-paste(thisPath,"SSRsens\\SSRsensTracers\\V",Version,sep="")

baseBiolFile<-paste(basePath,"CRAM_BH_hybrid_biol.prm",sep=""); biolLines<-readLines(baseBiolFile)
thisCode="PFS"
steepnessSens<-seq(0.5,0.9,by=0.1); nh<-length(steepnessSens)
fishingLevels<-c(0,50,80,100,120,150); nscenarios<-length(fishingLevels)
binFolders<-c("bin", "binTESTrecR0"); nbins<-length(binFolders)

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

colbyScenario<-c(colorRampPalette(colors=c(myBlue,myLightBlue))(3), "black", colorRampPalette(colors=c(myGold, myOrange))(2))
colBySteepness<-colorRampPalette(colors=c(myGold, myGreen,myDarkGreen, "black"))(nh)  
colBySSR<-c(myOrange, myBlue)

##set up lookup array for steepness, and SSR
nruns<-nbins *  nh
lookup_df<-data.frame(matrix(NA, nrow=nruns, ncol=3)); colnames(lookup_df)<-c("SSR", "h", "run")
lookup_df$SSR<-rep(sort(rep(binFolders,nh))); lookup_df$SSRcolor<-colBySSR[match(lookup_df$SSR,binFolders)]
lookup_df$h<-rep(seq(1,nh),nbins ); lookup_df$hcolor<-colBySteepness[match(lookup_df$h,seq(1,nh))]
lookup_df$run<-seq(1,nruns)

## include catch scenarios
# nruns<-nbins * nscenarios * nh
# lookup_df<-data.frame(matrix(NA, nrow=nruns, ncol=4)); colnames(lookup_df)<-c("scenario","SSR", "h", "run")
# lookup_df$scenario<-sort(rep(fishingLevels, nh * nbins)); lookup_df$scenarioColor<-colbyScenario[match(lookup_df$scenario,fishingLevels)]
# lookup_df$SSR<-rep(sort(rep(binFolders,nh)),nscenarios); lookup_df$SSRcolor<-colBySSR[match(lookup_df$SSR,binFolders)]
# lookup_df$h<-rep(seq(1,nh),nbins * nscenarios); lookup_df$hcolor<-colBySteepness[match(lookup_df$h,seq(1,nh))]
# lookup_df$run<-seq(1,nruns)

base_h<-0.7
thisVar<-"BHalpha_PFS"; xx<-biolLines[grep(thisVar, biolLines)]; baseAlpha<-get_first_number(xx)
thisVar<-"BHbeta_PFS"; xx<-biolLines[grep(thisVar, biolLines)]; baseBeta<-get_first_number(xx)
ux<-5*base_h-1; thisB0<-(baseBeta*ux)/(1-base_h); thisR0<-(baseAlpha * ux)/ (4*base_h)

thisBin<-lookup_df$SSR[1]; this_hindex<-lookup_df$h[1]
this_h<-steepnessSens[this_hindex]
thisScenario<-""; 
thisOut<-paste("outputSSR",thisBin,thisScenario,"h",this_hindex,sep="")
thisPath<-paste(DIR$'Base',"ATLANTISmodels\\SSRsens",versionText,"\\",thisOut,"\\",sep="")

BaseNC.nc<-nc_open(paste(thisPath,"\\output.nc", sep=""))
thisVol<-ncvar_get(BaseNC.nc,"volume")
nts<-dim(thisVol)[3]
nyears<-nts; nlayers<-dim(thisVol)[1]

allTracers<-unique(names(BaseNC.nc$var))
xx<-grep("_N",allTracers); yy<-grep("_Nums",allTracers)
Ntracers<-allTracers[xx[!(xx %in% yy)]]; ntracers<-length(Ntracers)

groupsDF<-read.csv(paste(basePath,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]

storeTracersArray<-array(NA,dim=c(nruns, ntracers, nyears)) #tracers is the biomass
storeNumbersArray<-array(NA,dim=c(nruns, ntracers,10, nyears)) ; storeWeightsArray<-storeNumbersArray
storeYOYarray<-0*storeTracersArray

for(r in 1:nruns){
  thisBin<-lookup_df$SSR[r]; this_hindex<-lookup_df$h[r]
  this_h<-steepnessSens[this_hindex]
  thisScenario<-""; 
  # thisScenario<-lookup_df$scenario[r]; ## catch scenarios
  
  thisOut<-paste("outputSSR",thisBin,thisScenario,"h",this_hindex,sep="")
  
  thisPath<-paste(DIR$'Base',"ATLANTISmodels\\SSRsens",versionText,"\\",thisOut,"\\",sep="")
  ThisNC.nc<-nc_open(paste(thisPath,"output.nc", sep=""))
  thisVol<-ncvar_get(ThisNC.nc,"volume"); thisData<-ncvar_get(ThisNC.nc,"Pelagic_fish_sml_N")
  
  for(t in 1:ntracers){
    thisTracer<-Ntracers[t]; thisName<-gsub("_N","",thisTracer); xx<-grep(thisName,groupsDF$Name)
    thisCode<-groupsDF$Code[xx]; thisNumCohorts<-groupsDF$NumCohorts[xx]
    #do biomass for all
    thisData<-ncvar_get(ThisNC.nc,thisTracer)
    if(length(dim(thisData))==3){
      xx<-apply(thisData*thisVol,3,sum)*mg_2_tonne*X_CN
    } else{
      xx<-apply(thisData * thisVol[nlayers,,], 2, sum) * mg_2_tonne * X_CN
    }
    storeTracersArray[r,t,]<-xx[1:nyears]
    if(thisNumCohorts>1){
      for(c in 1:thisNumCohorts){
        thisTracer<-paste(thisName,c,"_Nums",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
        xx<-apply(thisData,3,sum, na.rm=TRUE)
        storeNumbersArray[r,t,c,]<-xx[1:nyears]
        thisTracer<-paste(thisName,c,"_ResN",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
        xx<-apply(thisData,3,nonZeroMean)
        thisTracer<-paste(thisName,c,"_StructN",sep=""); thisData<-ncvar_get(ThisNC.nc, thisTracer)
        yy<-apply(thisData,3,nonZeroMean)
        storeWeightsArray[r,t,c,]<-((xx+yy)*mg_2_tonne * X_CN)[1:nyears]
      }
      
      YOYfile<-paste(thisPath,"outputYOY.txt",sep="")
      if(file.exists(YOYfile)){
        testYOY<-readLines(YOYfile)
        if(length(testYOY)>0){
          YOY_df<-read.csv(YOYfile,sep=" "); 
          YOY<-YOY_df[,grep(thisCode,colnames(YOY_df))]
          thisKWRR<-get_KWRR(thisCode); thisKWSR<-get_KWSR(thisCode); thisRecruitWeight<-thisKWRR+thisKWSR
          thisRecruits<-YOY[-1]/(thisRecruitWeight*X_CN * mg_2_tonne)
          storeYOYarray[r,t,]<-thisRecruits[1:nyears]
        }
      }
    }
  }
}

#write them out - commented out so only write out if intending to
# save(list=c("storeTracersArray", "storeNumbersArray", "storeWeightsArray", "storeYOYarray", "lookup_df"),file=paste(dataOutPath,"modelTracers",sep=""))

# 
# load(paste(dataOutPath,"modelTracers",sep="")); ##to bring "storeTracersArray", "storeNumbersArray", "storeWeightsArray", "storeYOYarray"
# ## . dim = dim=c(nruns, ntracers, nyears), c(nruns, ntracers,10, nyears), c(nruns, ntracers,10, nyears), c(nruns, ntracers, nyears) resp.

# Version<-"C"; versionText<-"C"
# dataOutPath<-paste(thisPath,"SSRsens\\SSRsensTracers\\V",Version,sep="")

## catches too

g=grep("Pelagic_fish_sml_N", Ntracers)
testPFS<-storeTracersArray[,g,]

colByRun<-colorRampPalette(colors=c(myGold, myOrange,"red",myPurple,myBlue,myAqua))(10)

par(mfrow=c(2,1))
plot(testPFS[1,5:dim(testPFS)[2]], type="n", ylim=c(0, max(testPFS, na.rm=TRUE)))
for(i in 1:5){
  points(testPFS[i,5:dim(testPFS)[2]], type="l", col=lookup_df$hcolor[i])
}
abline(h=1.5e+6,col=myGrey_trans,lty=1, lwd=3)
abline(h=thisB0*mg_2_tonne*X_CN,col="red",lty=2)
plot(testPFS[1,], type="n", ylim=c(0, max(testPFS, na.rm=TRUE)))
for(i in 6:10){
  points(testPFS[i,], type="l", col=lookup_df$hcolor[i])
}
abline(h=thisB0*mg_2_tonne*X_CN,col="red",lty=2)
abline(h=1.5e+6,col=myGrey_trans,lty=1, lwd=3)
