###############################################################################
#
#   NEW VERSION  MARCH/APRIL 2018 : Thomas BRUNEL and Niels HINTZEN (WMR)
#
#   SAM is fully integrated (assessment in the loop, uncertainty on starting conditions and model parameters incorporated in the OM)
#   now presented in different modules to facilitate implementation
#
###############################################################################

# developped from

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



#==============================================================================
# libraries 
#==============================================================================
rm(list=ls())
library(FLash)
#library(FLasher)
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
source('./Code/by_modulesTMB/MSE_functions.R')


#==============================================================================
# Basic configuration of the simulations
#==============================================================================




species <- "ANCHOVY"  # "ANCHOVY" or "SARDINE"
assess.name <- "Anchovy GSA 17-18_tbmSAM"
###
#species <- "SARDINE"  # "ANCHOVY" or "SARDINE"
#assess.name <- "Sardine GSA 17-18_tbmSAM"
#####


### to create new starting conditions
update.objects <- F

if (update.objects)
{
# for full
number.years.simulated  <- 20
number.replicates.stock <- 250
##for short#
#number.years.simulated  <- 12
#number.replicates.stock <- 2
##

# Load assessment data
check.assess <- F
source('./Code/by_modulesTMB/01_ADRIAMED_MSE_Load assessment data.r')
# Set up objects and configuration of the MSE
source('./Code/by_modulesTMB/02_ADRIAMED_MSE_Set up objects and configuration.r')
save.image(file=paste0("./Results/",species,"/",assess.name,"_",it,"iters_",ny,"yrs_blank_objects_MSE.RData"))
}

# to read in existing starting conditions (need to us this for the final run of simulatoins
# to make sure they use the same conditionning.

# choose for full or short MSE
run <- "short"

if(run == "full")  fname <-  paste0("./Results/",species,"/",assess.name,"_250iters_20yrs_blank_objects_MSE.RData")
if(run == "short") fname <-  paste0("./Results/",species,"/",assess.name,"_2iters_12yrs_blank_objects_MSE.RData")

load(fname)

# Define BRPs and Management Scenarios
source('./Code/by_modulesTMB/03_ADRIAMED_MSE_BRPs and Scenarios.r')
save.image(file=fname)
# Save the environment at the start of the simulations




# run the simulation
#scenario <- c("F.low")

strt <- proc.time()
for (sc in scenario)  source('./Code/by_modulesTMB/04_ADRIAMED_MSE_Simulations.r')       # run the simulation for each management scenario
proc.time() - strt



