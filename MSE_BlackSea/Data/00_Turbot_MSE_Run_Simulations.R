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
###This needs to be run only for the batch mode, not if you run one scenario at the time in R
args=(commandArgs(TRUE))
print(args)
argsOrig <- args
scenario <- strsplit(args,"=")[[1]][2]
###########

rm(list=ls())
library(FLash)
#library(FLasher)
library(FLAssess)
library(ggplotFL)
library(FLBRP)
library(FLSAM)  #devtools::install_github("flr/FLSAM", ref="develop_V2")
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
#save(file=Turbot_250iters_20yrs_blank_objects_MSE.RData)

# path to local github repository
#setwd("C:/Users/brune001/my git files/MSE-Adriatic/")
setwd("D:/Max files for backup/Documents/Commitees/GFCM/Turbot MSE/MSE SAM/")
# source needed functions 
source('./MSE_functions.R')

#==============================================================================
# Basic configuration of the simulations
#==============================================================================

species <- "Turbot"  
assess.name <- "Black Sea turbot TMB"

#Trial with one scenario
#scenario <- "CTACNOIUU"
#sc <- "CTACNOIUU"
 
### to create new starting conditions; first time you need to create the empty worspace and run it as F but run the setting in the loop. If you increase the number of iterations, then run it as T
update.objects <- F

if (update.objects)
{
  # for full
  number.years.simulated  <- 20
  number.replicates.stock <- 250
  ##for short#
  number.years.simulated  <- 8
  number.replicates.stock <- 2
  ##
  
  # Load assessment data
  check.assess <- F
  source('./01_Turbot_MSE_Load_assessment_data.r')
  # Set up objects and configuration of the MSE
  source('./02_Turbot_MSE_Set_up_objects_and_configuration.r')
  save.image(file=paste0("./Results/",species,"/",assess.name,"_",it,"iters_",ny,"yrs_blank_objects_MSE.RData"))
}

# to read in existing starting conditions (need to us this for the final run of simulatoins
# to make sure they use the same conditionning.

# choose for full or short MSE
#run <- "full"
run <- "short"

if(run == "full")  fname <-  paste0("./Results/",species,"/",assess.name,"_250iters_20yrs_blank_objects_MSE.RData")
if(run == "short") fname <-  paste0("./Results/",species,"/",assess.name,"_2iters_8yrs_blank_objects_MSE.RData")

load(fname)

if(run == "full")  sname <-  paste0("./Results/",species,"/",assess.name,"_",scenario,"_250iters_20yrs_blank_objects_MSE.RData")
if(run == "short") sname <-  paste0("./Results/",species,"/",assess.name,"_2iters_8yrs_blank_objects_MSE.RData")

# Define BRPs and Management Scenarios
source('./03_Turbot_MSE_BRPs_and_Scenarios.r')
save.image(file=sname)
# Save the environment at the start of the simulations

# run the simulation (only if running the batch file)
scenario <- strsplit(argsOrig,"=")[[1]][2]
strt <- proc.time()

for (sc in 1:length(management.scenarios))  source('./04_Turbot_MSE_Simulations.r')       # run the simulation for each management scenario
proc.time() - strt

