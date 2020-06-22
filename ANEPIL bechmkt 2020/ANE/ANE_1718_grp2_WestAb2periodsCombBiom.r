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
ANE <- setPlusGroup(ANE,2)   # extremely low age 3 catches

# remove this empty column for 2019
#ANE   <- trim(ANE,year=2000:2016)
units(harvest(ANE))  <- "f"

# sop(ANE,"catch") # check for the SOP correction

#Set stock object name - this is propagated through into the figure titles
ANE@name    <- "Anchovy - Adriatic Sea - GSA 17 and 18"

save(ANE , file = file.path(data.source,"ANEstk2plus.RData"))

### ======================================================================================================
### Prepare index object for assessment
### ======================================================================================================
#Load and modify all index data
ANE.tun   <- readFLIndices(file.path(data.source, "Ane17_18_flWestAb2PeriodsCombBiom.dat"))

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
CbBiomass <- FLIndices("Combined Biomass" = biomSurf)

## Trim indices to have the required ages for echo west and remove echo E
ANE.tun                         <- FLIndices("Echo West1"=trim(ANE.tun[[1]],age=0:1),
                                             "Echo West2"=trim(ANE.tun[[2]],age=1:2),
                                             "CbBiomass"= CbBiomass[[1]])

save(ANE.tun , file = file.path(data.source,"ANEAne17_18_flWestAbCombBiom.RData"))

### ======================================================================================================
### Prepare FLSAM object
### ======================================================================================================
ANE.ctrl <- FLSAM.control(ANE,ANE.tun)

#All fishing mortality states are free except
#oldest ages to ensure stablity
  ANE.ctrl@states["catch unique",]            <- c(1,2,3)

#Correlated Random walks for fishing mortalities - Default = FALSE = independent)
  ANE.ctrl@cor.F                              <- 0

# Catchabilities
ANE.ctrl@catchabilities["Echo West1",1:2]        <- c(0,1)
ANE.ctrl@catchabilities["Echo West2",2:3]        <- c(1,3)
ANE.ctrl@catchabilities["CbBiomass",1]  <- c(5)
##

#Fishing mortality RWs are set from an analysis of ICA VPA results
# ANE.ctrl@f.vars["catch unique",]            <- c(1,2,2)

#Set the variances. Separate variance for recruitment and plus group
ANE.ctrl@logN.vars[]                        <- c(0,1,2)

#Bind the observation variances
#ANE.ctrl@obs.vars["catch unique",]          <- c(201,202,203)
#ANE.ctrl@obs.vars["Echo West1",1:3]          <- c( 101, 102,102) #rep(101,3) #c(101, rep(102,2)) #c(101,102,103)
#ANE.ctrl@obs.vars["Echo West2",1:3]          <- c( 101, 102,102) #rep(101,3) #c(101, rep(102,2)) #c(101,102,103)
#ANE.ctrl@obs.vars["Echo East Biomass",1]    <- 301
##
ANE.ctrl@biomassTreat <- 1          # total biomass index

ANE.ctrl@residuals <- T
#
ANE.ctrl                                    <- update(ANE.ctrl)

### ======================================================================================================
### Perform the assessment
### ======================================================================================================
ANE.sam             <- FLSAM(ANE,ANE.tun,ANE.ctrl,return.fit=T,newtonsteps=0)

stockassessment::partable(ANE.sam)

ANE.sam2             <- FLSAM(ANE,ANE.tun,ANE.ctrl)


ANE.ctrl@residuals  <- FALSE
ANE.retro           <- retro(ANE,ANE.tun,ANE.ctrl,4)


 pdf("./results/assessment resultsWestAb2perdiosCombBiom.pdf")
 source('createAssessmentPlots.r')

