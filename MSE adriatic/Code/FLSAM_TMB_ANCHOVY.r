
#==============================================================================
# libraries and constants
#==============================================================================
rm(list=ls())


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

#===============================================================================
# Read data                                                                    |
# Inputs and outputs of the Anchovy SAM accepted at 2017 GFCM WGSAD            |
#===============================================================================

load("./Data/ANCHOVY/Anchovy GSA 17-18 (1).RData")


stk                           <- ANCHOVY
sam                           <- ANCHOVY.sam
ids                           <- ANCHOVY.tun
sam.ctrl                      <- ANCHOVY.ctrl

# need to reformat the control object to macht with the new version of FLSAM
       sam.ctrl.new <- FLSAM.control(stk,ids)
#         for (slt in slotNames(sam.ctrl.new))  try( slot(sam.ctrl.new,slt) <- slot(sam.ctrl,slt))


sam.ctrl.new@states["catch unique",]                      <- sam.ctrl@states["catch",]
sam.ctrl.new@logN.vars[]                                  <- sam.ctrl@logN.vars[]
sam.ctrl.new@catchabilities["Echo West",ac(0:4)]          <- sam.ctrl@catchabilities["Echo West",ac(0:4)]
sam.ctrl.new@catchabilities["Echo East",ac(0:2)]          <- sam.ctrl@catchabilities["Echo East",ac(0:2)]
sam.ctrl.new@catchabilities["Echo East Biomass",ac(0)]    <- 8

sam.ctrl.new@f.vars["catch unique",]                      <- sam.ctrl@f.vars["catch",]
sam.ctrl.new@obs.vars["catch unique",ac(0:4)]             <- sam.ctrl@obs.vars["catch",ac(0:4)]
sam.ctrl.new@obs.vars["Echo West",ac(0:4)]                <- sam.ctrl@obs.vars["Echo West",ac(0:4)]
sam.ctrl.new@obs.vars["Echo East",ac(0:2)]                <- sam.ctrl@obs.vars["Echo East",ac(0:2)]
sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]          <- 12
sam.ctrl.new@cor.F <- 0
sam.ctrl.new@residuals <- F
sam.ctrl.new <- update(sam.ctrl.new)

# just a quick check that we can reproduce the assessment

#  ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new)
#  ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))
#  
#  ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=FLSAM2par(ANCHOVY2.sam))
#   matplot(t(rbind(ssb(ANCHOVY2.sam)[,2],ssb(ANCHOVY.sam)[,2])),type="l")
#
#
##- Fix logN.vars first
#sam.ctrl.new@logN.vars[]            <- 0
#ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))#FAIL
#
##- Try to bind some obs.var paramters together + logN.vars
#sam.ctrl.new@obs.vars["catch unique",ac(0:4)]             <- c(0,1,2,2,2)
#sam.ctrl.new@obs.vars["Echo West",ac(0:4)]                <- c(0,0,1,1,1) + 101
#sam.ctrl.new@obs.vars["Echo East",ac(0:2)]                <- c(0,1,2) + 201
#sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]          <- 301
#sam.ctrl.new                                              <- update(sam.ctrl.new)
#ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))
#
##- Drop echo east age 2
#ids.new                                                   <- ids
#ids.new[["Echo East"]]                                    <- trim(ids.new[["Echo East"]],age=0:1)
#sam.ctrl.new@catchabilities["Echo East",ac(2)]            <- -1
#sam.ctrl.new@obs.vars["Echo East",ac(2)]                  <- -1
#sam.ctrl.new                                              <- update(sam.ctrl.new)
#ANCHOVY2.sam <- FLSAM(stk,ids.new,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))
#
##- put cor.F to 2
#sam.ctrl.new@cor.F <- 2
#ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))
#
##- Set selection at oldest age same as oldest age -1
#sam.ctrl.new@states["catch unique",ac(4)] <- 3
#ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=list(logF=log(stk@harvest[,drop=T]),logN=log(stk@stock.n[,drop=T])))
#
#
#
##- Fix logN.vars first
#sam.ctrl.new@logN.vars[]                                  <- 0
#sam.ctrl.new@obs.vars["catch unique",ac(0:4)]             <- c(0,1,2,2,2)
#sam.ctrl.new@obs.vars["Echo West",ac(0:4)]                <- c(0,0,1,1,1) + 101
#sam.ctrl.new@obs.vars["Echo East",ac(0:2)]                <- c(0,1,2) + 201
#sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]          <- 301
#sam.ctrl.new                                              <- update(sam.ctrl.new)
#ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new)
#
#- Final model
sam.ctrl.new@logN.vars[]                                  <- 0
sam.ctrl.new@states["catch unique",ac(4)]                 <- 3
sam.ctrl.new                                              <- update(sam.ctrl.new)
ANCHOVY2.sam                                              <- FLSAM(stk,ids,sam.ctrl.new)

ANCHOVY2     <- ANCHOVY + ANCHOVY2.sam
ANCHOVY2.tun <- ids
ANCHOVY2.ctrl<- sam.ctrl.new

save(ANCHOVY2,ANCHOVY2.sam,ANCHOVY2.tun,ANCHOVY2.ctrl ,file ="./Data/ANCHOVY/Anchovy GSA 17-18_tbmSAM.RData" )

sstk <- monteCarloStockTMB ( ANCHOVY2 , ANCHOVY2.tun , ANCHOVY2.sam , 2 )


#- Impact matrix
#- obs.var -> failure
#- logN -> failure
#- logN + obs.var -> stable but long runtime
#- cor.F=2 -> extreme values, no proper convergence
#- logN + cor.F=2 -> extreme values, no proper convergence
#- logN + obs.var + cor.F=2 -> markedly different numbers
#- statesp -> failure
#- statesp + logN -> good conversion, fast runtime
#- statesp + logN + obs.var -> good conversion, fast runtime
#
par(mfrow=c(3,1),mar=c(1.5,5,1.5,1.5),oma=c(3,3,3,1))
matplot(y=t(rbind(ssb(ANCHOVY2.sam)[,2],ssb(ANCHOVY.sam)[,2])),x=1975:2016,type="l",main="SSB",ylab="",las=1,xlab="")
mtext(side=2,text="Tonnes",line=5)
title("Anchovy",line=2,cex=2,outer=T)
grid()
matplot(t(rbind(fbar(ANCHOVY2.sam)[,2],fbar(ANCHOVY.sam)[,2])),x=1975:2016,type="l",main="Fbar",ylab="",las=1,xlab="")
mtext(side=2,text="Year^-1",line=5)
grid()
matplot(t(rbind(rec(ANCHOVY2.sam)[,2],c(ANCHOVY.sam@stock.n[1,]))),x=1975:2016,type="l",main="Recruitment",ylab="",las=1,xlab="")
mtext(side=2,text="Thousands",line=5)
grid()

