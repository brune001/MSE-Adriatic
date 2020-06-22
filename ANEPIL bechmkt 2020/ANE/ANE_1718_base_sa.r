######################################################################################################
# Malin FLSAM Assessment
####################################################################################################

rm(list=ls()); gc(); graphics.off(); start.time <- proc.time()[3]
options(stringsAsFactors=FALSE)

library(FLSAM)
library(FLCore)
library(ggplotFL)

### ======================================================================================================
### Define parameters for use in the assessment code here
### ======================================================================================================
#path                <- "C:/Users/Lusseaus/Documents/ICES Working Groups/HAWG/2016/Repository/wg_HAWG/VIa/"
#path                <- "C:/Users/s.angelini/CNR/GFCM 2017/SMALL PELAGICS/FINAL RUN/ANE_results_mat0.5age0"
#try(setwd(path),silent=TRUE)

data.source         <- file.path(getwd(),"data")      #Data source, not code or package source!!!
#data.source <- file.path("C:/Users/s.angelini/CNR/GFCM MEETING/GFCM_Small pelagics benchmark_2019_2020/Adriamed_2020/ANE/files ANE/")    #Data source, not code or package source!!!
# indices.data        <- file.path("./data")
#indices.data        <- file.path("C:/Users/s.angelini/CNR/GFCM MEETING/GFCM_Small pelagics benchmark_2019_2020/Adriamed_2020/ANE/files ANE/")
output.dir          <- file.path(".","results")       #Output directory
output.base         <- file.path(output.dir)          #Output base filename, including directory. Other output filenames are built by appending onto this one

### ======================================================================================================
### Prepare stock object for assessment
### ======================================================================================================
ANE                         <- readFLStock(file.path(data.source, "Ane17_18.ndx"),no.discards=TRUE)
ANE@catch.n                 <- ANE@landings.n
ANE@catch                   <- ANE@landings
ANE@catch.wt                <- ANE@discards.wt   <- ANE@landings.wt
ANE@stock.wt                <- ANE@catch.wt
m.spwn(ANE)                 <-0.5
harvest.spwn(ANE)           <-0.5 
# ANE@stock.wt[1,ac(2013)]    <- yearMeans(ANE@stock.wt[1,ac(2010:2012)])
# ANE@stock.wt[1,ac(2015)]    <- yearMeans(ANE@stock.wt[1,ac(2010:2014)])
# ANE@stock.wt[1,ac(2016)]    <- yearMeans(ANE@stock.wt[1,ac(2013:2015)])
# ANE@stock.wt[1,ac(2017)]    <- yearMeans(ANE@stock.wt[1,ac(2013:2015)])
units(ANE)[1:17]            <- as.list(c(rep(c("tonnes","thousands","kg"),4), rep("NA",5)))

#Set fbar
range(ANE)[c("minfbar","maxfbar")] <- c(1,2)

#Set plus group
#ANE <- setPlusGroup(ANE,3)   # no anchovy older than age 3

# remove this empty column for 2019
#ANE   <- trim(ANE,year=2000:2016)
units(harvest(ANE))  <- "f"

# sop(ANE,"catch") # check for the SOP correction

#Set stock object name - this is propagated through into the figure titles
ANE@name    <- "Anchovy - Adriatic Sea - GSA 17 and 18"

save(ANE , file = file.path(data.source,"ANEstk.RData"))

### ======================================================================================================
### Prepare index object for assessment
### ======================================================================================================
#Load and modify all index data
ANE.tun   <- readFLIndices(file.path(data.source, "Ane17_18_fl.dat"))

## Set the Index Type
ANE.tun[[1]]@type               <- "number"
ANE.tun[[2]]@type               <- "number"
ANE.tun[[3]]@type               <- "biomass"

#Set names
names(ANE.tun)                  <- lapply(ANE.tun,name)


## Set the plus group or specify there isnt one
ANE.tun[[1]]@range["plusgroup"] <- NA
ANE.tun[[2]]@range["plusgroup"] <- NA
ANE.tun[[3]]@range["plusgroup"] <- NA

# Creates the FLIndex for the Echo West Biomass
dmns     <-  dimnames(ANE.tun[[3]]@index)
dmns$age <- "all"
dmns     <- FLQuant(NA,dimnames=dmns)
#dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2016,unit="unique",season="all",area="unique",iter=1))
biomSurf   <- FLIndex(index=dmns)
biomSurf@index[]  <- ANE.tun[[3]]@index

range(biomSurf)["startf"] <- range(ANE.tun[[3]])["startf"]
range(biomSurf)["endf"]   <- range(ANE.tun[[3]])["endf"]
type(biomSurf)  <- "biomass"
EchoEastBiomass <- FLIndices("Echo East Biomass" = biomSurf)

## Trim indices to have the required ages for echo west and remove echo E
ANE.tun                         <- FLIndices("Echo West"=ANE.tun[[1]],
                                             "Echo East"=ANE.tun[[2]],
                                             "Echo East Biomass"= EchoEastBiomass[[1]])

save(ANE.tun , file = file.path(data.source,"ANEtun.RData"))

### ======================================================================================================
### Prepare FLSAM object
### ======================================================================================================
ANE.ctrl <- FLSAM.control(ANE,ANE.tun)

#All fishing mortality states are free except
#oldest ages to ensure stablity
ANE.ctrl@states["catch unique",]            <- c(1,2,3,3)

#Correlated Random walks for fishing mortalities - Default = FALSE = independent)
ANE.ctrl@cor.F                              <- 2

# Catchabilities
ANE.ctrl@catchabilities["Echo West",1:2]        <- c(1,2)
# ANE.ctrl@catchabilities["Echo East",1:3]       <- c(4,5,6)
ANE.ctrl@catchabilities["Echo East Biomass",1]  <- c(7)


#Fishing mortality RWs are set from an analysis of ICA VPA results
ANE.ctrl@f.vars["catch unique",]            <- c(1,2,2,3)
# ANE.ctrl@f.vars["catch unique",]            <- c(1,2,3,3) #c(1,rep(2,8))
# ANE.ctrl@f.vars["catch unique",]            <- c(1,1,2,2) #c(1,rep(2,8))
# ANE.ctrl@f.vars["catch unique",]            <- c(1,2,1,1)

#Set the variances. Separate variance for recruitment and plus group
ANE.ctrl@logN.vars[]                        <- c(0,1,1,1)

#Bind the observation variances
ANE.ctrl@obs.vars["catch unique",]          <- c(1,2,3,4)
#ANE.ctrl@obs.vars["Echo West",1:3]          <- c(101, 101, 101) #rep(101,3) #c(101, rep(102,2)) #c(101,102,103)
#ANE.ctrl@obs.vars["Echo East",1:3]          <- c(201,201,201)
#ANE.ctrl@obs.vars["Echo East Biomass",1]    <- 301
#
ANE.ctrl                                    <- update(ANE.ctrl)

### ======================================================================================================
### Perform the assessment
### ======================================================================================================
ANE.sam             <- FLSAM(ANE,ANE.tun,ANE.ctrl)
ANE.ctrl@residuals  <- FALSE
ANE.retro           <- retro(ANE,ANE.tun,ANE.ctrl,4)
save(ANE,ANE.sam,ANE.tun,ANE.ctrl,ANE.retro,file=file.path(output.dir,"IBP_VIaHerring_baserun.Rdata"))

run_name            <- "baserun"
source("./createAssessmentPlots.r")
#- Additional diagnostics and runs
xyplot(value ~ an(name) | fleet,group=age,data=catchabilities(ANE.sam),type="b",pch=19,scales=list(y="free"))

MSH.loi <- looi(MSH,MSH.tun,MSH.ctrl,type="loo",MSH.sam)

MSH.noMSHERAS_retro  <- retro(MSH,MSH.tun[-1],drop.from.control(MSH.ctrl,fleet="MS HERAS"),7)
MSH.noWoSHERAS_retro <- retro(MSH,MSH.tun[-2],drop.from.control(MSH.ctrl,fleet="WoS HERAS"),7)
MSH.noIBTS_Q1_retro  <- retro(MSH,MSH.tun[-3],drop.from.control(MSH.ctrl,fleet="IBTS_Q1"),7)
MSH.noIBTS_Q4_retro  <- retro(MSH,MSH.tun[-4],drop.from.control(MSH.ctrl,fleet="IBTS_Q4"),7)