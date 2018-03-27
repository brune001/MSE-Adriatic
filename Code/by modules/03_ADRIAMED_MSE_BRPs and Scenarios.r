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

#  define BRPs
if (species == "ANCHOVY")
{
blim  <- 45936
bpa   <- blim*2
Btrig <- blim+((bpa-blim)/2)
Fmsy  <- 0.55                                                   
}

if (species == "SARDINE")
{
blim  <- 125318
bpa   <- blim*2
Btrig <- blim+((bpa-blim)/2)
Fmsy  <- 0.72                                                   
}


# define the function for the HCRs

# to apply a constant F
 HCR.cstF <- function(stoc, target , year)
                  {
                  targ<-list()
                  targ$quant = "f"
                  targ$val = c(fbar(stoc)[,ac(range(stoc)["maxyear"])] , mgt.target["Ftarget"])
                  targ$rel = c(NA,NA)
                  return(targ)
                  }
                  
# to apply a reduction in F towards a Ftarget to reach at a given year 
# each year, the annual percentage of reduction needed to reach the target
# is updated based on the latest assessment
 HCR.time.target.F <- function(stoc, target)
                  {
                  ylast <- range(stoc)["maxyear"]    #  last assessment year
                  ytarg <- target[["Yr.targ"]]
                  flast <- fbar(stoc)[,ac(ylast)]   
                  perc.red  <- (flast - Fmsy ) / (ytarg - ylast) # percentage reduction
                  mult      <-  1-perc.red                       # multiplier to apply
                  
                  targ<-list()
                  targ$quant = "f"
                  targ$val = mult
                  targ$rel = ac(ylast:(ylast+1))
                  return(targ)
                  }                  

##########################################
# define scenarios 1 by 1 

#-----scenario 1 : F status quo ------------------------

# basis Ftarget = Fsq used to give advice starting in 2018
# Fsq = F2014-2016


Fsq <- list( name = "Fsq" ,
             target = list(Ftarget = mean(c(fbar(stk)[,ac(dy)]))) ,
             HCR =  HCR.cstF
           )    



#-----scenario 2 : 	Fmsy2020 ------------------------
# basis  : Linear reduction of F towards FMSY by 2020
# 


Fmsy2020 <- list( name = "Fmsy2020" ,
                  target = list(Yr.targ = 2020 , Ftarget = Fmsy) ,
                  HCR =  HCR.time.target.F
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

scenarios <- list(Fsq,mgt2)
names(scenarios) <- lapply (scenarios , function(x) x[[1]])

