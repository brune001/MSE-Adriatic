
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
         for (slt in slotNames(sam.ctrl.new))  try( slot(sam.ctrl.new,slt) <- slot(sam.ctrl,slt))

# just a quick check that we can reproduce the assessment
  ANCHOVY2.sam <- FLSAM(stk,ids,sam.ctrl.new)
   rbind(ssb(ANCHOVY2.sam)[,2],ssb(ANCHOVY.sam)[,2])

