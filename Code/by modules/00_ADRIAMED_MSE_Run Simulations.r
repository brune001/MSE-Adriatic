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




update.objects <- T

if (update.objects)
{
# Load assessment data
source('./Code/by modules/01_ADRIAMED_MSE_Load assessment data.r')
# Set up objects and configuration of the MSE
source('./Code/by modules/02_ADRIAMED_MSE_Set up objects and configuration.r')
save.image(file=paste0("./Results/",species,"/MSE_",assess.name,"_blank_objects_MSE_",".RData"))
}


# Define BRPs and Management Scenarios
load(paste0("./Results/",species,"/MSE_",assess.name,"_blank_objects_MSE_",".RData"))
source('./Code/by modules/03_ADRIAMED_MSE_BRPs and Scenarios.r')
# Save the environment at the start of the simulations




# run the simulation
source('./Code/by modules/04_ADRIAMED_MSE_Simulations workin version.r')      # load the function go_fish

lapply(management.scenarios , function(scen) go_fish(scen$name))    # run the simulation for each management scenario



