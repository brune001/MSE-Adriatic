#-------------------------------------------------------------------------------
# Stock assessment script for 2019 Black Sea Turbot benchmark by Niels Hintzen (niels.hintzen@wur.nl)
#
# Last modified: 26/06/2020. vanja
#
#-------------------------------------------------------------------------------
library(FLCore)
library(FLSAM) 
library(FLEDA) 

#- Set your working directory paths
dataPath  <- "./Data/"
outPath   <- "./Results/"
codePath  <- "./"
#setwd("C:/Users/Vanja/Desktop/PROCJENE 2020/SARDINE 2019/data/run9")
#data.source <- file.path("C:/Users/Vanja/Desktop/PROCJENE 2020/SARDINE 2019/data/run9")  
#
#-------------------------------------------------------------------------------
#- Read in the stock data
#-------------------------------------------------------------------------------

Sar                 <- readFLStock(file.path(dataPath, "Sar17_18.ndx"),no.discards=TRUE)
Sar@catch           <- Sar@landings
units(Sar)[1:17]    <- as.list(c(rep(c("tonnes","thousands","kg"),4),
                                 rep("NA",2),"f",rep("NA",2)))
Sar@discards.n[]    <- 0; Sar@discards.wt[] <- 0; Sar@discards <- computeDiscards(Sar)


#Set object details
Sar@name                              <- "Sardine GSA 17-18"
range(Sar)[c("minfbar","maxfbar")]    <- c(1,3)
Sar                                   <- setPlusGroup(Sar,Sar@range["max"])

save(Sar,file="PILstk.RData")

mat(Sar)[1,]<-1      # change maturity so that SSB and TSB are the same 




#- Read in the tuning data
Sar.tun             <- readFLIndices(file.path(dataPath, "Sar17_18_fl.dat")) 
Sar.tun             <- lapply(Sar.tun,function(x) {x@type <- "number"; return(x)})
## Set the Index Type
Sar.tun[[1]]@type               <- "number"
Sar.tun[[2]]@type               <- "number"
Sar.tun[[3]]@type               <- "biomass"


# Creates the FLIndex for the Echo West Biomass
dmns     <-  dimnames(Sar.tun[[3]]@index)
dmns$age <- "all"
dmns     <- FLQuant(NA,dimnames=dmns)
#dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2016,unit="unique",season="all",area="unique",iter=1))
biomSurf   <- FLIndex(index=dmns)
biomSurf@index[]  <- Sar.tun[[3]]@index

range(biomSurf)["startf"] <- range(Sar.tun[[3]])["startf"]
range(biomSurf)["endf"]   <- range(Sar.tun[[3]])["endf"]
type(biomSurf)  <- "biomass"
EchoEastBiomass <- FLIndices("Echo East Biomass" = biomSurf)


## Trim indices to have the required ages
Sar.tun                         <- FLIndices("Echo West"=Sar.tun[[1]],
                                             "Echo East"=Sar.tun[[2]],
                                             "Echo East Biomass"= EchoEastBiomass[[1]])

#
## run data analyses ####
save(Sar.tun,file="PILtun.RData")





### Niels's version
Sar.tun[[1]] <- trim(Sar.tun[[1]],age=0:3) #age 4 has a lot of missing data from 2014 onwards (and 2013 is very low already)
Sar.tun[[2]]  <- trim(Sar.tun[[2]],age=0:2) #age 3 has only 1 datapoint so logical to drop


Sar.ctrl                    <- FLSAM.control(Sar,Sar.tun)
Sar.ctrl@f.vars[1,] <- c(0,1,1,2,2)
Sar.ctrl@residuals                                <-  F
Sar.ctrl                                          <- update(Sar.ctrl)

fit       <- FLSAM(Sar,Sar.tun,Sar.ctrl,return.fit=T) #model converges, residuals all estimates


Sar.ctrl2     <-Sar.ctrl
Sar.ctrl2@residuals                                <-  T
Sar.ctrl2                                          <- update(Sar.ctrl2)
SAR.sam       <- FLSAM(Sar,Sar.tun,Sar.ctrl2) #model converges, residuals all estimates




 #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
 # try fixing the process error
   init.sam          <- SAR.sam
   init.sam@params$value[which(init.sam@params$name=="logSdLogN")[2]] <- log(0.2)   
   fit0.2          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)

   SAR.sam0.2          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

runname <- "all surveys fixed PE0.2"
assess<- Sar
mat(assess)[1,]  <- 0
assess.sam <-SAR.sam0.2
assess.tun <- Sar.tun
assess.ctrl<-  Sar.ctrl2
source("createAssessmentPlots.r")


### something is wrong about the fit to echo east : all residuals have the same signe for age 1 , 
# this indicates that catchabilities need to be different for the ages 1 and 2

Sar.ctrl                    <- FLSAM.control(Sar,Sar.tun)
Sar.ctrl@f.vars[1,] <- c(0,1,1,2,2)
Sar.ctrl@catchabilities['Echo East' , ac(0:2)]  <- c(0,1,2) +101
Sar.ctrl@residuals                                <-  F
Sar.ctrl                              <- update(Sar.ctrl)
Sar.ctrl2     <-Sar.ctrl
Sar.ctrl2@residuals                                <-  T
Sar.ctrl2                                          <- update(Sar.ctrl2)

    fit0.2QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)

   SAR.sam0.2QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))


runname <- "all surveys fixed PE0.2 decoupled Q echo East"
assess<- Sar
mat(assess)[1,]  <- 0
assess.sam <-SAR.sam0.2QechoE
assess.tun <- Sar.tun
assess.ctrl<-  Sar.ctrl2
source("createAssessmentPlots.r")              ### seems much better



#######################3
#test sensitivity to the fixed PE values : compare 0.15,0.2,0.22
# I tried 0.1 and 0.25 0.3 but the model did converge... 
#  0.1 : Rvar = 0
#  0.25 : Rvar badly defined
#  0.3 :  Convergence failed



init.sam2<-init.sam
init.sam2@params$value[which(init.sam2@params$name=="logSdLogN")[2]] <- log(0.15)   
fit0.15QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam2,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)
SAR.sam0.15QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam2,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

init.sam3<-init.sam
init.sam3@params$value[which(init.sam3@params$name=="logSdLogN")[2]] <- log(0.22)
fit0.22QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam3,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)
SAR.sam0.22QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam3,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

restoplot<-FLSAMs(SAR.sam0.15QechoE,SAR.sam0.2QechoE,SAR.sam0.22QechoE,SAR.sam) 
names(restoplot) <- c("PE=0.15","PE=0.20","PE=0.22","PE free")                            
plot(restoplot)











#
#
# #### check with decoupled obsvar for the catches
#Sar.ctrl@obs.vars['catch unique' , ]  <- c(0,0,1,1,2) +101
#Sar.ctrl                              <- update(Sar.ctrl)
#Sar.ctrl2     <-Sar.ctrl
#Sar.ctrl2@residuals                                <-  T
#Sar.ctrl2                                          <- update(Sar.ctrl2)
#
#    fit0.2obsC          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
#                              starting.values=init.sam,
#                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
#                              return.fit=T)
#
#   SAR.sam0.2obsC          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
#                              starting.values=init.sam,
#                              map=list(logSdLogN=as.factor(c(-1.5,NA))))
#
#
#
#









#-------------------------------------------------------------------------------
#- Setup the control file
#-------------------------------------------------------------------------------
mat(Sar)[1,]<-1

Sar.ctrl                    <- FLSAM.control(Sar,Sar.tun2)

Sar.ctrl@f.vars[1,] <- c(0,1,1,2,2)
#- Set the number of catchabilities to estimate by index and age
#Sar.ctrl@catchabilities["Echo West",ac(0:3)]      <- c(1,2,3,4)          + 101
Sar.ctrl@residuals                                <-  F
Sar.ctrl                                          <- update(Sar.ctrl)


#-------------------------------------------------------------------------------
# Run the model (2 step approach)
#-------------------------------------------------------------------------------

rm(fit) ; fit            <- FLSAM(Sar,Sar.tun2,Sar.ctrl,return.fit=T)

Sar.ctrl@residuals                                <-  T
Sar.ctrl                                          <- update(Sar.ctrl)
Sar.sam            <- FLSAM(Sar,Sar.tun2,Sar.ctrl)

runname <- "all surveys only conf1"
source("createAssessmentPlots.r")
Sar.ctrl1<-Sar.ctrl


#-------------------------------------------------------------------------------
#- Setup the control file
#-------------------------------------------------------------------------------

Sar.ctrl                    <- FLSAM.control(Sar,Sar.tun2)
#- Set the number of catchabilities to estimate by index and age
Sar.ctrl@catchabilities["Echo West",ac(0:3)]      <- c(1,2,3,4)          + 101


#- Update the control object which automatically scales the parameter numbering
Sar.ctrl@residuals                                <-  F
Sar.ctrl                                          <- update(Sar.ctrl)


#-------------------------------------------------------------------------------
# Run the model (2 step approach)
#-------------------------------------------------------------------------------

rm(fit) ; fit            <- FLSAM(Sar,Sar.tun2,Sar.ctrl,return.fit=T)

Sar.sam            <- FLSAM(Sar,Sar.tun2,Sar.ctrl)

Sar.ctrl@residuals                                <-  T
Sar.ctrl                                          <- update(Sar.ctrl)
Sar.sam            <- FLSAM(Sar,Sar.tun2,Sar.ctrl)

runname <- "echoW only conf2"
source("createAssessmentPlots.r")




























#-------------------------------------------------------------------------------
# Additional runs
#-------------------------------------------------------------------------------
Sar.ctrl@residuals  <- FALSE

Sar.LOO             <- list()
for(i in 1:length(Sar.tun))
  Sar.LOO[[names(Sar.tun)[i]]]             <- FLSAM(Sar,Sar.tun[-i],drop.from.control(Sar.ctrl,fleets=names(Sar.tun)[i]))
Sar.LOO <- as(Sar.LOO,"FLSAMs")

#- Leave one in
Sar.LOI             <- list()
for(i in 1:length(Sar.tun))
 Sar.LOI[[names(Sar.tun)[i]]]             <- FLSAM(Sar,Sar.tun[i],drop.from.control(Sar.ctrl,fleets=names(Sar.tun)[-i]),
                                                    starting.values=init.sam,map=list(logSdLogN=as.factor(c(-1.5,NA))))
Sar.LOI <- as(Sar.LOI,"FLSAMs")
Sar.LOI             <- list()
for(i in 1:length(Sar.tun))
  Sar.LOI[[names(Sar.tun)[i]]]             <- FLSAM(Sar,Sar.tun[i],drop.from.control(Sar.ctrl,fleets=names(Sar.tun)[-i]))
Sar.LOI <- as(Sar.LOI,"FLSAMs")

Sar.retro <- list()
Sar.retro[[ac(2018)]] <- Sar.sam
for(iRetroYr in 2017:2015){
  rt.stck <- window(Sar,end=iRetroYr)
  rt.tun  <- Sar.tun
  for(iTun in names(rt.tun)){
    if(range(rt.tun[[iTun]])["maxyear"] >= iRetroYr)
      rt.tun[[iTun]] <- rt.tun[[iTun]][,-which(dimnames(rt.tun[[iTun]]@index)$year %in% ((iRetroYr+1):range(rt.tun[[iTun]])["maxyear"]))]
  }
  rt.ctrl <- Sar.ctrl
  rt.ctrl@range["maxyear"] <- iRetroYr
  rt.sam <- FLSAM(rt.stck,rt.tun,rt.ctrl)
  Sar.retro[[ac(iRetroYr)]] <- rt.sam
}


Sar.retro           <- as(Sar.retro,"FLSAMs")
mohns.rho(Sar.retro,ref.year=2018,type="rec",span=3)
mohns.rho(Sar.retro,ref.year=2018,type="ssb",span=3)
mohns.rho(Sar.retro,ref.year=2018,type="fbar",span=3)

storeMohnsRhos      <- matrix(NA,nrow=1,ncol=3,dimnames=list("modelRun",type=c("ssb","fbar","rec")))
storeMohnsRhos[1,1] <- mean(mohns.rho(Sar.retro,ref.year=2018,type="ssb",span=3)[1:3,1])
storeMohnsRhos[1,2] <- mean(mohns.rho(Sar.retro,ref.year=2018,type="fbar",span=3)[1:3,1])
storeMohnsRhos[1,3] <- mean(mohns.rho(Sar.retro,ref.year=2018,type="rec",span=3)[1:3,1])
print(storeMohnsRhos)


save(Sar.LOO,Sar.LOI,Sar.retro,file=file.path(outPath,"Saradditional.RData"))


AIC(Sar.sam)

