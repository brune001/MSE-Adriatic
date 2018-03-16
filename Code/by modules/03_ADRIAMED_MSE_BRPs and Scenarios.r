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



##########################################
# define scenarios 1 by 1 

#-----scenario mgt1------------------------
blim <- 45936
bpa <- blim*2
Btrig <- blim+((bpa-blim)/2)
fsq <- mean(c(fbar(stk)[,ac(dy)]))                                                   

mgt1 <- list(     name = "mgt1" ,
                  BRP = c(
                      Blim  = blim,
                      Bpa   = blim*2 ,
                      Btrig = blim+((bpa-blim)/2) ,
                      Ftarget = mean(c(fbar(stk)[,ac(dy)]))
                              ) ,
                  HCR = c(Btrig = NA , Ftarget = Ftarget , TAC_IAV = NA )
                  )    




#-------scenario mgt2----------------------
blim <- 45936
bpa <- blim*2
Btrig <- bpa
fsq <- mean(c(fbar(stk)[,ac(dy)]))                                                   

mgt2 <- list(     name = "mgt2" ,
                  BRP = c(
                      Blim  = blim,
                      Bpa   = blim*2 ,
                      Btrig = blim+((bpa-blim)/2) ,
                      Ftarget = mean(c(fbar(stk)[,ac(dy)]))
                              ) ,
                  HCR = c(Btrig = NA , Ftarget = Ftarget , TAC_IAV = NA )
                  )    



########################################
# combine them in a list

scenarios <- list(mgt1,mgt2)
names(scenarios) <- lapply (scenarios , function(x) x[[1]])

