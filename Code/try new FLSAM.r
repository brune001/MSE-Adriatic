
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
# Inputs and outputs of the Anchovy SAM accepted at 2017 GFCM WGSAD
#==============================================================================

load("./Data/ANCHOVY/Anchovy GSA 17-18 (1).RData")


stk                           <- ANCHOVY
sam                           <- ANCHOVY.sam
ids                           <- ANCHOVY.tun
sam.ctrl                      <- ANCHOVY.ctrl

# need to reformat the control object to macht with the new version of FLSAM
       sam.ctrl.new <- FLSAM.control(stk,ids)
#         for (slt in slotNames(sam.ctrl.new))  try( slot(sam.ctrl.new,slt) <- slot(sam.ctrl,slt))


sam.ctrl.new@states["catch unique",] <- sam.ctrl@states["catch",]
sam.ctrl.new@states["catch unique",ac(4)] <- 3
sam.ctrl.new@logN.vars[]             <- sam.ctrl@logN.vars[]
sam.ctrl.new@logN.vars[5]            <- 3

sam.ctrl.new@catchabilities["Echo West",ac(0:4)]          <- sam.ctrl@catchabilities["Echo West",ac(0:4)]
sam.ctrl.new@catchabilities["Echo East",ac(0:2)]          <- sam.ctrl@catchabilities["Echo East",ac(0:2)]
sam.ctrl.new@catchabilities["Echo East Biomass",ac(0)]    <- 101

sam.ctrl.new@f.vars["catch unique",]  <- sam.ctrl@f.vars["catch",]

sam.ctrl.new@obs.vars["catch unique",ac(0:4)]       <- sam.ctrl@obs.vars["catch",ac(0:4)]
sam.ctrl.new@obs.vars["Echo West",ac(0:4)]          <- sam.ctrl@obs.vars["Echo West",ac(0:4)]
sam.ctrl.new@obs.vars["Echo West",ac(4)]            <- 6
sam.ctrl.new@obs.vars["Echo East",ac(0:2)]          <- sam.ctrl@obs.vars["Echo East",ac(0:2)]
sam.ctrl.new@obs.vars["Echo East Biomass",ac(0)]    <- 101

sam.ctrl.new@cor.F <- 0
sam.ctrl.new <- update(sam.ctrl.new)
# just a quick check that we can reproduce the assessment

  ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new)
  ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new,starting.values=FLSAM2par(ANCHOVY2.sam))
   rbind(ssb(ANCHOVY2.sam)[,2],ssb(ANCHOVY.sam)[,2])

