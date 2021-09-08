###### STECF EWG 17-09
# Split, 23 - 29 September 2017

# source("http://flr-project.org/R/instFLR.R")

library(FLCore)
library(FLEDA)
library(FLXSA)
library(FLAssess)
library(FLash)
library(FLBRP) 
library(FLSAM)
# library(SQLiteFL)
library(doBy)
library(reshape)

### ============================================================================
### Misc
### ============================================================================
data.source <- file.path("C:/Users/s.angelini/CNR/STECF_small pelagics 2017/LAST STOCK OBJ LIGAS/ANE_sa/")    #Data source, not code or package source!!!

### ============================================================================
### Prepare stock object for assessment
### ============================================================================
#Load object
ANCHOVY           <- readFLStock(file.path(data.source, "Ane17_18.ndx"),no.discards=TRUE)

#Catch is calculated from: catch.wt * catch.n, however, the reported landings are
#normally different (due to SoP corrections). Hence we overwrite the calculate landings
ANCHOVY@catch           <- ANCHOVY@landings
units(ANCHOVY)[1:17]    <- as.list(c(rep(c("tonnes","thousands","tonnes"),4), 
                                     rep("NA",2),"f",rep("NA",2)))

#Set object details
ANCHOVY@name                              <- "Anchovy GSA 17_18"
range(ANCHOVY)[c("minfbar","maxfbar")]    <- c(1,2)
ANCHOVY                                   <- setPlusGroup(ANCHOVY,ANCHOVY@range["max"])

save(ANCHOVY, file='ANCHOVY.RData') # save the stock object

### ============================================================================
### Prepare index object for assessment
### ============================================================================
#Load and modify all numbers at age data
# load("C:/Users/PC09/Documents/Work/Adriatic/GFCM_2016/Anchovy/PC/Run_2idx_same2015/ANCHOVY2idx.RData)
ANCHOVY.tun   <- readFLIndices(file.path(data.source,"Ane17_18_fl.dat"))
ANCHOVY.tun   <- lapply(ANCHOVY.tun,function(x) {x@type <- "number"; return(x)})
names(ANCHOVY.tun) <- c("Echo West", "Echo East", 'Echo East Biomass')
ANCHOVY.tun[["Echo West"]]@range["plusgroup"] <- 4 # for newpg run 2
ANCHOVY.tun[["Echo East"]]@range["plusgroup"] <- 4 # for newpg run 1
ANCHOVY.tun[["Echo East Biomass"]]@range["plusgroup"] <- NA

# 2 iNDEXES
ANCHOVY.tun   <- readFLIndices(file.path(data.source,"Ane17_18_fl.dat"))
ANCHOVY.tun   <- lapply(ANCHOVY.tun,function(x) {x@type <- "number"; return(x)})
names(ANCHOVY.tun) <- c("Echo West", "Echo East")
ANCHOVY.tun[["Echo West"]]@range["plusgroup"] <- NA
ANCHOVY.tun[["Echo East"]]@range["plusgroup"] <- NA


# Creates the FLIndex for the Echo East Biomass
dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2012,unit="unique",season="all",area="unique",iter=1))
#dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2016,unit="unique",season="all",area="unique",iter=1))
biomSurf   <- FLIndex(index=dmns)
biomSurf@index[]  <- c(6223, 81866, 132340, 142089, 56488, 110290, 122170, 166325, 46472, 11639)
range(biomSurf)["startf"] <- 0.75
range(biomSurf)["endf"] <- 0.83
type(biomSurf)  <- "biomass"
EchoEastBiomass <- FLIndices(EchoEastBiomass=biomSurf)

# ANCHOVY.tun <- FLIndices(c(ANCHO.tun[["Echo West"]], ANCHO.tun[["Echo East"]],EchoEastBiomass[["EchoEastBiomass"]])) 

ANCHOVY.tun <- FLIndices(list(ANCHOVY.tun[["Echo West"]], ANCHOVY.tun[["Echo East"]],EchoEastBiomass[["EchoEastBiomass"]]))
names(ANCHOVY.tun)<-c("Echo West","Echo East","Echo East Biomass")



### ============================================================================
### Apply plusgroup to all data sets
### ============================================================================
pg <- 4

#- This function already changes the stock and landings.wts correctly
ANCHOVY <- setPlusGroup(ANCHOVY,pg)


######## ONLY ONE SURVEY
#ANCHOVY1idx.tun   <- readFLIndices(file.path(data.source,"Ane17_18_fl_1idx.dat"))
#ANCHOVY1idx.tun   <- lapply(ANCHOVY1idx.tun,function(x) {x@type <- "number"; return(x)})
#names(ANCHOVY1idx.tun) <- c("Echo WestEast")
#ANCHOVY1idx.tun[["Echo WestEast"]]@range["plusgroup"] <- NA
#ANCHOVY.tun[["Echo West18"]]@range["plusgroup"] <- NA



#################################
#### set up the control object
#################################

# library(FLSAM)

ANCHOVY.ctrl <- FLSAM.control(ANCHOVY,ANCHOVY.tun)

#Set the variances. Separate variance for recruitment and plus group
ANCHOVY.ctrl@logN.vars[] <- c(1, 2, rep(3, 2), 4)
ANCHOVY.ctrl@f.vars["catch",] <- c(rep(1,2), rep(2,3))

#All fishing mortality states are free except
#oldest ages to ensure stablity
ANCHOVY.ctrl@states["catch",] <- seq(dims(ANCHOVY)$age)
# ANCHOVY.ctrl@states["catch",ac(0:4)] <- 101

#Group observation variances of catches to ensure stability
ANCHOVY.ctrl@obs.vars["catch",ac(0)]  <- 201
ANCHOVY.ctrl@obs.vars["catch",ac(1)]  <- 202
ANCHOVY.ctrl@obs.vars["catch",ac(2)]  <- 203
ANCHOVY.ctrl@obs.vars["catch",ac(3:4)]  <- 204

ANCHOVY.ctrl@obs.vars["Echo West",ac(0)]    <- 205
ANCHOVY.ctrl@obs.vars["Echo West",ac(1)]  <- 206
#ANCHOVY.ctrl@obs.vars["Echo West",ac(2)]  <- 206 # for echosurvey age 1-4
ANCHOVY.ctrl@obs.vars["Echo West",ac(2:3)]  <- 207
ANCHOVY.ctrl@obs.vars["Echo West",ac(4)]  <- 208
#ANCHOVY.ctrl@obs.vars["Echo West",ac(3:4)]  <- 207 # for echosurvey age 1-1

#ANCHOVY.ctrl@obs.vars["Echo East",ac(0:2)]    <- 209
ANCHOVY.ctrl@obs.vars["Echo East",ac(0)]    <- 209
ANCHOVY.ctrl@obs.vars["Echo East",ac(1)]    <- 210
ANCHOVY.ctrl@obs.vars["Echo East",ac(2)]    <- 211
#ANCHOVY.ctrl@obs.vars["Echo East",ac(1)]    <- 208 # for echosurvey age 1-4
#ANCHOVY.ctrl@obs.vars["Echo East",ac(2)]    <- 209 # for echosurvey age 1-4
#ANCHOVY.ctrl@obs.vars["Echo East",ac(3)]    <- 210 # for echosurvey age 1-4
# ANCHOVY.ctrl@obs.vars["Echo East",ac(3:4)]  <- 209 # 210

# ANCHOVY.ctrl@obs.vars["Echo East Biomass",ac(1)]  <- 211


#Group catchability parametesr
ANCHOVY.ctrl@catchabilities["Echo West",ac(0:4)]  <- c(rep(101,1),rep(102,1),rep(103,1),rep(104,2))
ANCHOVY.ctrl@catchabilities["Echo East",ac(0:2)]  <- c(rep(105,1),rep(106,1),rep(107,1))
#ANCHOVY.ctrl@catchabilities["Echo West",ac(1:4)]  <- c(rep(101,1),rep(102,1),rep(103,1),rep(104,1)) # echosurvey age 1-4
#ANCHOVY.ctrl@catchabilities["Echo East",ac(1:3)]  <- c(rep(104,1), rep(105,1), rep(106,1)) # echosurvey age 1-4


#Finalise
ANCHOVY.ctrl@name <- "Final Assessment"
ANCHOVY.ctrl <- update(ANCHOVY.ctrl)




### ============================================================================
### Run the assessment
### ============================================================================

#Perform assessment
ANCHOVY.sam <- FLSAM(ANCHOVY,ANCHOVY.tun,ANCHOVY.ctrl, run.dir="c:\\samtemp\\")
name(ANCHOVY.sam) <- "Anchovy GSA1718 2017"

#Update stock object
ANCHOVY      <- ANCHOVY + ANCHOVY.sam
ANCHOVY@stock <- computeStock(ANCHOVY)

# Save results
output.dir          <-  file.path("C:/Users/s.angelini/CNR/STECF_small pelagics 2017/LAST STOCK OBJ LIGAS/ANE_sa/results_newpg/")
save(ANCHOVY,ANCHOVY.tun,ANCHOVY.ctrl,ANCHOVY.sam,file=file.path(output.dir,paste(name(ANCHOVY),".RData",sep="")))


#ane15 <- load("C:/Users/PC09/Documents/Work/Adriatic/GFCM_SAF/Results/Anchovy_SAMassessmentGFCM2014/ANCHOVY1718_Run6.RData")
#ane15sam <- get(ane15[1])

#comparison <- FLStocks(ane15=ane15sam, ane16=ANCHOVY)
#tiff("PC/comparison15-16.tiff", width = 1800, height = 2000, res=200)
#plot(comparison) + theme_light(base_size = 16)
#dev.off()

###### Comparison between the assessment
# ane15 <- load("C:/Users/PC09/Documents/Work/Adriatic/GFCM_SAF/Results/Anchovy_SAMassessmentGFCM2014/ANCHOVY1718_Run6.RData")
# ane15sam <- get(ane15[1])

# ane16_2idx <- ANCHOVY
# ane16_2idx_tun <- ANCHOVY.tun
# ane16_2idx_sam <- ANCHOVY.sam
# ane16_2idx_ctrl <- ANCHOVY.ctrl

# comparison <- FLStocks(ane15=ane15sam, ane16=ANCHOVY, ane16_1idx=ANCHOVY1idx)
# tiff("PC/comparison15-16.tiff", width = 1800, height = 2000, res=200)
# plot(comparison) + theme_light(base_size = 16)
# dev.off()


######GRAFOVI-FUNCTIONS:###################
### Data exploration plots
### ======================================================================================================

#Ratio of mature and immature biomass              
a

### ======================================================================================================
### Reference points
### ======================================================================================================
library(FLBRP)
ref. <- brp(FLBRP(ANCHOVY,fbar=seq(0,1,length.out=101),nyears=3))
print(refpts(ref.))

ANCHOVY.SRR <- FLSR(
  rec = rec(ANCHOVY)[,ac((range(ANCHOVY)["minyear"]+1): range(ANCHOVY)["maxyear"])],
  ssb = ssb(ANCHOVY)[,ac((range(ANCHOVY)["minyear"])  :(range(ANCHOVY)["maxyear"]-1))],
  model='segreg')
ANCHOVY.SRR <- fmle(ANCHOVY.SRR)
plot(ANCHOVY.SRR)

newData <- predict(ANCHOVY.SRR,ssb=FLQuant(seq(0,max(ssb(ANCHOVY)),length.out=200)))
yrange  <- range(pretty(c(0,range(rec(ANCHOVY)))))/1e6; xrange <- range(pretty(c(0,range(ssb(ANCHOVY)))))/1e6
plot(y=newData/1e6,x=seq(0,max(ssb(ANCHOVY)),length.out=200)/1e6,type="l",lwd=2,
     xlab="SSB (million tonnes)",ylab="Recruitment (billions)",xlim=xrange,ylim=yrange,
     las=1,cex.lab=1.3,cex.axis=1.1,xaxs="i",yaxs="i")
points(y=rec(ANCHOVY)/1e6,x=ssb(ANCHOVY)/1e6)

dev.off()

### ======================================================================================================
### Document Assessment
### ======================================================================================================
##Create a standard table
stockSummaryTable <- cbind(rec(ANCHOVY.sam)$year,
                           rec(ANCHOVY.sam)$value,      rec(ANCHOVY.sam)$lbnd,    rec(ANCHOVY.sam)$ubnd,
                           tsb(ANCHOVY.sam)$value,      tsb(ANCHOVY.sam)$lbnd,    tsb(ANCHOVY.sam)$ubnd,
                           ssb(ANCHOVY.sam)$value,      ssb(ANCHOVY.sam)$lbnd,    ssb(ANCHOVY.sam)$ubnd,
                           catch(ANCHOVY.sam)$value,    catch(ANCHOVY.sam)$lbnd,  catch(ANCHOVY.sam)$ubnd,
                           catch(ANCHOVY.sam)$value / ssb(ANCHOVY.sam)$value, catch(ANCHOVY.sam)$lbnd / ssb(ANCHOVY.sam)$lbnd, catch(ANCHOVY.sam)$ubnd / ssb(ANCHOVY.sam)$ubnd,
                           fbar(ANCHOVY.sam)$value,     fbar(ANCHOVY.sam)$lbnd,   fbar(ANCHOVY.sam)$ubnd,
                           c(quantMeans(harvest(ANCHOVY.sam)[ac(0:1),])),
                           c(sop(ANCHOVY),NA))

colnames(stockSummaryTable) <-
  c("Year",paste(rep(c("Recruits Age 0 (Thousands)","Total biomass (tonnes)","Spawing biomass (tonnes)",
                       "Landings (tonnes)","Yield / SSB (ratio)","Mean F ages 2-6"),each=3),c("Mean","Low","High")),"Mean F ages 0-1","SoP (%)")
stockSummaryTable[nrow(stockSummaryTable),] <- NA
stockSummaryTable[nrow(stockSummaryTable),"Spawing biomass (tonnes) Mean"] <- 2271364
stockSummaryTable[nrow(stockSummaryTable),2:4] <- c(rec(ANCHOVY.sam)$value[nrow(rec(ANCHOVY.sam))],rec(ANCHOVY.sam)$lbnd[nrow(rec(ANCHOVY.sam))],rec(ANCHOVY.sam)$ubnd[nrow(rec(ANCHOVY.sam))])
write.csv(stockSummaryTable,file=file.path(output.dir,"stockSummaryTable.csv"))

### ============================================================================
### Finish
### ============================================================================


#RETROSPECTIVE
# dir.create("c:\\samtemp\\")
retro_ane <- retro(ANCHOVY,ANCHOVY.tun,ANCHOVY.ctrl, retro=1)
plot(retro_ane)
retro_sar <- retro_ane
save(retro_sar, file='retro_sar.Rdata')
getwd()
