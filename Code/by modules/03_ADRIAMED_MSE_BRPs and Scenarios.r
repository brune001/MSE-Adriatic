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


################################################################################
#  define BRPs
################################################################################


if (species == "ANCHOVY")
{
blim  <- 45936
bpa   <- blim*2
Fmsy  <- 0.55                                                   
}

if (species == "SARDINE")
{
blim  <- 125318
bpa   <- blim*2
Fmsy  <- 0.72                                                   
}



################################################################################
# define the function for the HCRs
################################################################################

# to apply a constant F
 HCR.cstF <- function(stoc, target)
                  {
                  targ<-list()
                  targ$quant = "f"
                  targ$val = list(y1 = c(fbar(stoc)[,ac(range(stoc)["maxyear"])]) , y2 = target$Ftarget)
                  targ$rel = c(NA,NA)
                  return(targ)
                  }
                  
# to apply a reduction in F towards a Ftarget to reach at a given year 
# each year, the annual percentage of reduction needed to reach the target
# is updated based on the latest assessment
 HCR.time.target.F <- function(stoc, target)
                  {
                  # define time and F 
                  ylast <- range(stoc)["maxyear"]    #  last assessment year
                  ytarg <- target[["Yr.targ"]]
                  flast <- fbar(stoc)[,ac(ylast)] 
                  ftarget <- target$Ftarget 
                  # reduction to apply
                  perc.red  <- (flast - ftarget ) / (ytarg - ylast) # percentage reduction
                  mult1 <-   1 - perc.red                       # multiplier to apply in the first year
                  mult2 <-   1 - 2*perc.red                     # multiplier to apply in the 2nd year, compared to the current F
                  if(ytarg == (ylast+1)) mult2 <- ftarget/flast         # in case is reached during the intermediate year
                  # build the target object                  
                  targ<-list()
                  targ$quant = "f"
                  targ$val = list(y1 = c(flast*mult1) , y2 = c(flast*mult2))
                  targ$rel = c(NA,NA)
                  
                  # if we are past the target year to reach Fmsy, just apply a flat Fmsy
                  if(ylast>=ytarg) targ <- HCR.cstF (stoc, target) 
                  
                  
                  return(targ)
                  }                  

# to apply a constant catch
 HCR.cstC <- function(stoc, target)
                  {     
                  targ<-list()
                  targ$quant = "catch"
                  targ$val = list(y1 = c(catch(stoc)[,ac(range(stoc)["maxyear"])]) , y2 = target$Ctarget)
                  targ$rel = c(NA,NA)
                  return(targ)
                  }
                  
 
# to apply the HCR defined by GFCM 2013 recommandation
# F=Fmax if B>Btrig
# F=0    if B<Blim
# F= Fmax * (B-Bpa)/(Bpa-Blim)  if blim <= B <= Btri

HCR.gfcm <- function(stoc, target)
                  {
                  # management points
                  Fmax      <- target$Fmax
                  Btrig     <- target$Btrig
                  ssb.now <- ssb(stoc)[,ac(range(stoc)["maxyear"])]
                  
                  #  compute the Ftarget
                  Ftarget   <- c(fbar(stoc)[,ac(range(stoc)["maxyear"])])
                  Ftarget[] <- Fmax  #empty object to store the result
                  
                  sliding.slope <-  c((ssb.now - blim) / (bpa-blim))
                  sliding.slope[sliding.slope<0]   <- 0
                  sliding.slope[sliding.slope>0.5] <- 1                  
                  
                  Ftarget <- Ftarget * sliding.slope
                  
                  targ<-list()
                  targ$quant = "f"
                  targ$val = list(y1 = c(fbar(stoc)[,ac(range(stoc)["maxyear"])]) , y2 = Ftarget)
                  targ$rel = c(NA,NA)
                  return(targ)
                  }
   
 
 
                  
                  
################################################################################
# define scenarios 1 by 1 
################################################################################



################################################################################
#----- Constant F scenarios ------------------------

# scenario F status quo : 
      # basis Ftarget = Fsq used to give advice starting in 2018
      # Fsq = F2014-2016
F.sq <- list( name = "F.sq" ,
             target = list(Ftarget = mean(c(fbar(stk)[,ac(dy)]))) ,
             HCR =  HCR.cstF ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
           )    



# scenario Fmsy :
      # basis Ftarget = Fmsy used to give advice starting in 2018
      
F.msy <- list( name = "F.msy" ,
             target = list(Ftarget = Fmsy) ,
             HCR =  HCR.cstF ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
           )    




################################################################################
#----- scenarios to reach a target F in a given year ------------------------


# Linear reduction of F towards FMSY by 2020
Fmsy2020 <- list( name = "Fmsy2020" ,
                  target = list(Yr.targ = 2020 , Ftarget = Fmsy) ,
                  HCR =  HCR.time.target.F,
                  spatial.closure = F ,
                  additionnal.F.reduction = NA
                )    


# Linear reduction of F towards FMSY by 2020
Fmsy2025 <- list( name = "Fmsy2025" ,
                  target = list(Yr.targ = 2025 , Ftarget = Fmsy) ,
                  HCR =  HCR.time.target.F,
                  spatial.closure = F ,
                  additionnal.F.reduction = NA
                )    


                 

################################################################################
#----- Constant catch scenarios ------------------------

# scenario Catch at 2014 level
C2014 <- list( name = "C2014" ,
             target = list(Ctarget = c(catch(stk)[,ac(2014)])) ,
             HCR =  HCR.cstC ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
           )    

# scenario Catch at historic lowest level
Chistmin <- list( name = "Chistmin" ,
             target = list(Ctarget = c(min(catch(stk)))) ,
             HCR =  HCR.cstC ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
           )    


################################################################################
#----- scenarios X% annual decrease in the catches untill Bpa is reached ------------------------

# scenario Catch at 2014 level
C5red <- list( name = "C5red" ,
             target = list(prec.red = 0.05, Btarget = Bpa) ,
             HCR =  HCR.Cred ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
           )    



################################################################################
#----- 	GFCM recommendation (2013)  ------------------------

# scenario based on the GFCM proposed haverst control rule

GFCM.HCR <- list( name = "GFCM.HCR" ,
             target = list(Fmax = Fmsy , Btrig = mean(c(bpa,blim))) ,
             HCR =  HCR.gfcm ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
           )    









########################################
# combine them in a list

# management.scenarios <- list(Fsq,Fmsy,Fmsy2020,Fmsy2020,Fmsy2025,C2014 , Chistmin)
management.scenarios <- list(Fmsy2020)
names(management.scenarios) <- lapply (management.scenarios , function(x) x[[1]])

