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
#
###############################################################################



#==============================================================================
# Read data
# Inputs and outputs of the Anchovy SAM accepted at 2017 GFCM WGSAD 
#==============================================================================

load(paste0("./Data/",species,"/",assess.name,".RData"))

# find out file name for the assessment
if (species == "ANCHOVY")
{
stk                           <- ANCHOVY
sam                           <- ANCHOVY.sam
ids                           <- ANCHOVY.tun
sam.ctrl                      <- ANCHOVY.ctrl   }

if (species == "SARDINE")
{
stk                           <- SARDINE
sam                           <- SARDINE.sam
ids                           <- SARDINE.tun
sam.ctrl                      <- SARDINE.ctrl   }


# need to reformat the control object to macht with the new version of FLSAM 
#       sam.ctrl.new <- FLSAM.control(stk,ids)
#       for (slt in slotNames(sam.ctrl.new))  try( slot(sam.ctrl.new,slt) <- slot(sam.ctrl,slt))

# just a quick check that we can reproduce the assessment
 
if (check.assess)
{ 
  sam2 <- FLSAM(stk,ids,sam.ctrl)
  plot(ssb(sam2)[,2],ssb(sam)[,2], main="comparison official assessmnet V.S local run")
}
