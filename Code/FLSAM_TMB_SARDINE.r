
#==============================================================================
# libraries and constants
#==============================================================================
rm(list=ls())

library(FLa4a)
library(FLash)
library(FLAssess)
library(ggplotFL)
library(FLBRP)
library(FLSAM)
library(FLCore)
library(FLEDA)
library(doBy)
library(reshape)
library(devtools)
#install_github("ices-tools-prod/msy")
library(msy)


#
# path to local github repository
setwd("C:/Users/brune001/my git files/MSE-Adriatic/")

# source needed functions
# source('C:/Users/MORELLO/Desktop/SP_MSE/PIL/PIL_feb_2017/MSE_funs_LAST.R')
source('./Code/MSE_funs_LAST.R')

#==============================================================================
# Read data
# Inputs and outputs of the SARDINE SAM accepted at 2017 GFCM WGSAD
#==============================================================================

load("./Data/SARDINE/SARDINE GSA 17_18 (1).RData")


stk                           <- SARDINE
sam                           <- SARDINE.sam
ids                           <- SARDINE.tun
sam.ctrl                      <- SARDINE.ctrl

# need to reformat the control object to macht with the new version of FLSAM
       sam.ctrl.new <- FLSAM.control(stk,ids)
#         for (slt in slotNames(sam.ctrl.new))  try( slot(sam.ctrl.new,slt) <- slot(sam.ctrl,slt))


sam.ctrl.new@states["catch unique",]                      <- sam.ctrl@states["catch",]
sam.ctrl.new@logN.vars[]                                  <- sam.ctrl@logN.vars[]
sam.ctrl.new@catchabilities["Echo West",ac(0:4)]          <- sam.ctrl@catchabilities["Echo West",ac(0:4)]
sam.ctrl.new@catchabilities["Echo East",ac(0:2)]          <- sam.ctrl@catchabilities["Echo East",ac(0:2)]
sam.ctrl.new@catchabilities["Echo East Biomass",ac(0)]    <- 101

sam.ctrl.new@f.vars["catch unique",]                      <- sam.ctrl@f.vars["catch",]
sam.ctrl.new@obs.vars["catch unique",ac(0:4)]             <- sam.ctrl@obs.vars["catch",ac(0:4)]
sam.ctrl.new@obs.vars["Echo West",ac(0:4)]                <- sam.ctrl@obs.vars["Echo West",ac(0:4)]
sam.ctrl.new@obs.vars["Echo East",ac(0:2)]                <- sam.ctrl@obs.vars["Echo East",ac(0:2)]
sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]          <- 101
sam.ctrl.new@cor.F <- 0
sam.ctrl.new@residuals <- F
sam.ctrl.new <- update(sam.ctrl.new)

# just a quick check that we can reproduce the assessment

  SARDINE2.sam <- FLSAM(stk,ids,sam.ctrl.new)
  SARDINE2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))
  
  SARDINE2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=FLSAM2par(SARDINE2.sam))
   matplot(t(rbind(ssb(SARDINE2.sam)[,2],ssb(SARDINE.sam)[,2])),type="l")


#- Fix logN.vars first
sam.ctrl.new@logN.vars[]            <- 0
SARDINE2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[-5,][,drop=T]),logN=log(stk@stock.n[,drop=T])))#FAIL

#- Try to bind some obs.var paramters together + logN.vars
sam.ctrl.new@obs.vars["catch unique",ac(0:4)]             <- 0
sam.ctrl.new@obs.vars["Echo West",ac(0:4)]                <- c(0,0,1,1,2) + 101
sam.ctrl.new@obs.vars["Echo East",ac(0:2)]                <- c(0,1,2) + 201
sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]          <- 301
sam.ctrl.new                                              <- update(sam.ctrl.new)
SARDINE2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[-5,][,drop=T]),logN=log(stk@stock.n[,drop=T])))

#- put cor.F to 2
sam.ctrl.new@cor.F <- 2
system.time(SARDINE2.sam <- FLSAM(stk,ids,sam.ctrl.new))

#- Final model
sam.ctrl.new@logN.vars[]                                  <- 0
sam.ctrl.new@obs.vars["catch unique",ac(0:4)]             <- 0
sam.ctrl.new@obs.vars["Echo West",ac(0:4)]                <- c(0,0,1,1,2) + 101
sam.ctrl.new@obs.vars["Echo East",ac(0:2)]                <- c(0,1,2) + 201
sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]          <- 301
sam.ctrl.new                                              <- update(sam.ctrl.new)
SARDINE2.sam                                              <- FLSAM(stk,ids,sam.ctrl.new)


#- Impact matrix
- logN -> failure
- obs.var -> failure
- logN + obs.var -> stable but with covar issues
- cor.F=2 -> failure
- logN + cor.F=2 -> stable no covar issues, fast,
- logN + obs.var + cor.F=2 -> stable, fast, lower AIC, no covar issues




par(mfrow=c(3,1),mar=c(1.5,5,1.5,1.5),oma=c(3,3,3,1))
matplot(y=t(rbind(ssb(SARDINE2.sam)[,2],ssb(SARDINE.sam)[,2])),x=1975:2016,type="l",main="SSB",ylab="",las=1,xlab="")
mtext(side=2,text="Tonnes",line=5)
title("SARDINE",line=2,cex=2,outer=T)
grid()
matplot(t(rbind(fbar(SARDINE2.sam)[,2],fbar(SARDINE.sam)[,2])),x=1975:2016,type="l",main="Fbar",ylab="",las=1,xlab="")
mtext(side=2,text="Year^-1",line=5)
grid()
matplot(t(rbind(rec(SARDINE2.sam)[,2],c(SARDINE.sam@stock.n[1,]))),x=1975:2016,type="l",main="Recruitment",ylab="",las=1,xlab="")
mtext(side=2,text="Thousands",line=5)
grid()

