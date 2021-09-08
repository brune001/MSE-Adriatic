  #-------------------------------------------------------------------------------
# Stock assessment script for 2019 Black Sea Turbot benchmark by Niels Hintzen (niels.hintzen@wur.nl)
#
# Last modified: 26/06/2020. vanja
#
#-------------------------------------------------------------------------------
library(FLCore)
library(FLSAM) 
library(FLEDA) 
library(stockassessment)
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

Sar                 <- readFLStock(file.path(dataPath, "Sar17_18new.ndx"),no.discards=TRUE)
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
Sar.tun             <- readFLIndices(file.path(dataPath, "Sar17_18_fl_separated_newALK.dat")) 
Sar.tun             <- lapply(Sar.tun,function(x) {x@type <- "number"; return(x)})
## Set the Index Type
Sar.tun[[1]]@type               <- "number"
Sar.tun[[2]]@type               <- "number"
Sar.tun[[3]]@type               <- "number"
Sar.tun[[4]]@type               <- "number"
Sar.tun[[5]]@type               <- "biomass"


# Creates the FLIndex for the Echo West Biomass
dmns     <-  dimnames(Sar.tun[[5]]@index)
dmns$age <- "all"
dmns     <- FLQuant(NA,dimnames=dmns)
#dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2016,unit="unique",season="all",area="unique",iter=1))
biomSurf   <- FLIndex(index=dmns)
biomSurf@index[]  <- Sar.tun[[5]]@index

range(biomSurf)["startf"] <- range(Sar.tun[[5]])["startf"]
range(biomSurf)["endf"]   <- range(Sar.tun[[5]])["endf"]
type(biomSurf)  <- "biomass"
name(biomSurf) <- "Echo East Biomass"
EchoEastBiomass <- FLIndices("Echo East Biomass" = biomSurf)

# trim the surveys
Sar.tun[[1]]  <- window(Sar.tun[[1]] , start = 2008)
Sar.tun[[2]]  <- window(Sar.tun[[2]] , start = 2008)


## Trim indices to have the required ages
Sar.tun                         <- FLIndices("Echo W17"=Sar.tun[[1]],
                                             "Echo W18"=Sar.tun[[2]],
                                             "Echo L"=Sar.tun[[3]],
                                             "Echo East"=Sar.tun[[4]] ,
                                             "Echo East Biomass"= EchoEastBiomass[[1]])



#
## run data analyses ####
save(Sar.tun,file="PILtun_trimmed.RData")





### Niels's version
Sar.tun[[2]] <- trim(Sar.tun[[2]],age=0:3) #age 4 has a lot of missing data from 2014 onwards (and 2013 is very low already)
Sar.tun[[3]]  <- trim(Sar.tun[[3]],age=0:3) #age 4 has only 1 datapoint so logical to drop
Sar.tun[[4]]  <- trim(Sar.tun[[4]],age=1:4) #remove age0

Sar.ctrl                    <- FLSAM.control(Sar,Sar.tun)
Sar.ctrl@states [1,]    <- c(0,1,2,3,3)
Sar.ctrl@f.vars[1,] <- c(1,1,1,2,2)
Sar.ctrl@cor.F <-2
Sar.ctrl@catchabilities["Echo W17",ac(0:4)]<- c(0,1,2,2,2)
Sar.ctrl@catchabilities["Echo W18",ac(0:3)]<- c(0,1,2,2) +101
Sar.ctrl@catchabilities["Echo L"  ,ac(0:3)]<- c(0,1,1,1) +201
Sar.ctrl@catchabilities["Echo East"  ,ac(1:4)]<- c(1,2,2,2) +301
Sar.ctrl@catchabilities["Echo East Biomass"  ,ac(0)] <- 501

Sar.ctrl@obs.vars["Echo W17",ac(0:4)]<- c(1,1,1,2,2)
Sar.ctrl@obs.vars["Echo W18",ac(0:3)]<- c(0,0,1,1) +101
Sar.ctrl@obs.vars["Echo L"  ,ac(0:3)]<- c(0,0,1,1) +201



Sar.ctrl@residuals                                <-  F
Sar.ctrl                                          <- update(Sar.ctrl)

fit       <- FLSAM(Sar,Sar.tun,Sar.ctrl,return.fit=T) #model converges, residuals all estimates


Sar.ctrl2     <-Sar.ctrl
Sar.ctrl2@residuals                                <-  F
Sar.ctrl2                                          <- update(Sar.ctrl2)
SAR.sam       <- FLSAM(Sar,Sar.tun,Sar.ctrl2) #model converges, residuals all estimates

runname <- "all surveys conf 1 trimmed"
assess<- Sar
mat(assess)[1,]  <- 0
assess.sam <-SAR.sam
assess.tun <- Sar.tun
assess.ctrl<-  Sar.ctrl2
library(stockassessment)
source("createAssessmentPlots.r")







 #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
 # try fixing the process error
   init.sam          <- SAR.sam   
   init.sam@params$value[which(init.sam@params$name=="logSdLogN")[2]] <- log(0.2)   
   fit0.2          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)

Sar.ctrl3     <-Sar.ctrl2
Sar.ctrl3@residuals                                <-  T
Sar.ctrl3                                          <- update(Sar.ctrl3)

   SAR.sam0.2          <- FLSAM(Sar,Sar.tun,Sar.ctrl3,
                              starting.values=init.sam,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

runname <- "all surveys fixed PE0.2 trimmed"
assess<- Sar
mat(assess)[1,]  <- 0
assess.sam <-SAR.sam0.2
assess.tun <- Sar.tun
assess.ctrl<-  Sar.ctrl3
library(stockassessment)
source("createAssessmentPlots.r")



######################3
#test sensitivity to the fixed PE values : compare 0.15,0.2,0.22
# I tried 0.1 and 0.25 0.3 but the model did converge... 
#  0.1 : Rvar = 0
#  0.25 : Rvar badly defined
#  0.3 :  Convergence failed



init.sam2<-init.sam
init.sam2@params$value[which(init.sam2@params$name=="logSdLogN")[2]] <- log(0.10)   
fit0.1QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam2,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)
SAR.sam0.1QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl3,
                              starting.values=init.sam2,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

init.sam3<-init.sam
init.sam3@params$value[which(init.sam3@params$name=="logSdLogN")[2]] <- log(0.3)
fit0.3QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam3,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)
SAR.sam0.3QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam3,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

restoplot<-FLSAMs(SAR.sam0.1QechoE,SAR.sam0.2,SAR.sam0.3QechoE,SAR.sam) 
names(restoplot) <- c("PE=0.1","PE=0.2","PE=0.3","PE free")                            
plot(restoplot)




init.sam05<-init.sam
init.sam05@params$value[which(init.sam05@params$name=="logSdLogN")[2]] <- log(0.05)
fit0.05QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam05,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)


init.sam025<-init.sam
init.sam025@params$value[which(init.sam025@params$name=="logSdLogN")[2]] <- log(0.025)
fit0.025QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam025,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)


init.sam0125<-init.sam
init.sam0125@params$value[which(init.sam0125@params$name=="logSdLogN")[2]] <- log(0.0125)
fit0.0125QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam0125,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)


init.sam0063<-init.sam
init.sam0063@params$value[which(init.sam0063@params$name=="logSdLogN")[2]] <- log(0.0063)
fit0.0063QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam0063,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)

init.sam00125<-init.sam
init.sam00125@params$value[which(init.sam00125@params$name=="logSdLogN")[2]] <- log(0.00125)
fit0.00125QechoE          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam00125,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)


res<-list(fit0.0125QechoE,fit0.025QechoE,fit0.05QechoE,fit0.1QechoE,fit0.2,fit0.3QechoE)
names(res)   <- c("0.0125" ,"0.025","0.05","0.1","0.2","0.3")
class(res)<- "samset" 
library(stockassessment)
par(mfrow=c(2,1))
ssbplot(res ) 
fbarplot(res)
#
#
# #### check with decoupled obsvar for the catches
Sar.ctrl@obs.vars['catch unique' , ]  <- c(0,1,1,1,2) +101
Sar.ctrl                              <- update(Sar.ctrl)
Sar.ctrl2     <-Sar.ctrl
Sar.ctrl2@residuals                                <-  T
Sar.ctrl2                                          <- update(Sar.ctrl2)

    fit0.2obsC          <- FLSAM(Sar,Sar.tun,Sar.ctrl,
                              starting.values=init.sam05,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))),
                              return.fit=T)

   SAR.sam0.2obsC          <- FLSAM(Sar,Sar.tun,Sar.ctrl2,
                              starting.values=init.sam05,
                              map=list(logSdLogN=as.factor(c(-1.5,NA))))

#
obv <- obs.var(SAR.sam0.2obsC)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv <- obv[order(obv$value),]
bp <- barplot(obv$value,ylab="Observation Variance",
              main="Observation variances by data source",col=factor(obv$fleet))
axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
#







