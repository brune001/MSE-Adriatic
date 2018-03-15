###############################################################################
# 
#EJ(20150119)
# Tests on the evaluation of the NS-MAP using COD
# NOTE1: The final analysis is in the report file.
#  One may tangle the Snw file to get the R script trimmed of comments.
#
# Modified from Piera Carpi, AgurtzPIL Urtizberea Ijurco & Miguel Bernal
# 20 January 2016
# Changed for ANCHOVY MSE for GFCM 
#
# Modified by Betulla Morello
# February/MARCH 2017
# Updated for 2017 GFCM WKMSE, using WGSAD2016 assessment data (last year=dy=2015)
# FINAL VERSION UPLOADED ONTO GFCM SERVER
#
# Scenario numbers in brackets are equivalent to scenario numbers in the report
#
#
#   NEW VERSION  MARCH 2018 : Thomas BRUNEL (WMR)
#   SAM is fully integrated (assessment in the loop, uncertainty on starting conditions and model parameters incorporated in the OM)
#   now presented in different modules to facilitate implementation
#
###############################################################################

#==============================================================================
# libraries 
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


#==============================================================================
# running repertory and paths 
#==============================================================================

# path to local github repository
setwd("C:/Users/brune001/my git files/MSE-Adriatic/")

# source needed functions 
source('./Code/MSE_funs_LAST.R')


#==============================================================================
# Basic configuration of the simulations
#==============================================================================

species <- "ANCHOVY"  # "ANCHOVY" or "SARDINE"
assess.name <- "Anchovy GSA 17-18 (1)"

number.years.simulated  <- 5
number.replicates.stock <- 2
management.scenario     <- "mgt1"




blim <- 45936
bpa <- 2*blim

# Load assessment data
source('./Code/by modules/01_ADRIAMED_MSE_Load assessment data.r')
# Set up objects and configuration of the MSE
source('./Code/by modules/02_ADRIAMED_MSE_Set up objects and configuration.r')
# Define BRPs and Management Scenarios



# NEED AN ASSUMPTION   TO FORECAST THE OM IN THE FIRST YEAR
TAC[,(ac(iy))] <- 15000 #  WHAT IS A REALISTIC VALUE ?
# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- 45936
bpa <- blim*2
Btrig <- blim+((bpa-blim)/2)
dt <- date()
fsq <- mean(c(fbar(stk)[,ac(dy)]))                                                   



# run the simulation
source('./Code/by modules/04_ADRIAMED_MSE_Simulations.r')



# breakpoint mean(SSB)
ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0 <- pstk
ANE_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0 <- stk0
plot(ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
save(list = c("ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0", "ANE_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0"), file = "Results/ANE/ANE_Opt1a_mse_SegregmeanSSB_statusQuo_250it.RData")
png("ANE_Opt1a_mse_SegregmeanSSB_statusQuo_250it.png", width=700, height=700)
plot(ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
dev.off()




