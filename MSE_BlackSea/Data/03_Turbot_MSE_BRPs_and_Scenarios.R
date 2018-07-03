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
#Modified by Max Cardinale for GFCM
################################################################################
#  define BRPs
################################################################################

if (species == "Turbot")
{
  blim  <- 3535
  bpa   <- 4949
  Fmsy  <- 0.26
  CTAC <- 644
  CTAC50 <- 945
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
  step.change  <- (ftarget - flast ) / (ytarg - ylast) # annual reduction
  targ1 <-   flast +  step.change                      # Ftarget first year : here we make the assumption in the STF that we applied the same decrease last year. not exactly true, but better than a status quo assumption for the intermediate year
  targ2 <-   flast + 2*step.change                     # Ftarget for the advice year
  if(ytarg == (ylast+1)) targ2 <- targ1                # in case is reached during the intermediate year
  # build the target object                  
  targ<-list()
  targ$quant = "f"
  targ$val = list(y1 = c(targ1) , y2 = c(targ2))
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

# to apply an annual X% decrease in the catches until B>=Bpa & F<=Fmsy
HCR.Cred <- function(stoc, target)
{     
  red <- target$prec.red
  btarg<-target$Btarget
  cref <- catch(stoc)[,ac(2014)]
  yref    <- ay
  yad   <- range(stoc)["maxyear"]+2
  ny<-yad-yref
  #current SSB
  ssbnow <- ssb(stoc)[,ac(range(stoc)["maxyear"])]
  
  if (ssbnow <  btarg)
  {
    targ<-list()
    targ$quant = "catch"
    targ$val = list(y1 = c(cref)*(1-red)^(ny-2) , y2 = c(cref)*(1-red)^(ny-1) )
    targ$rel = c(NA,NA)
  } else {               # apply Fmsy
    
    targ<-HCR.gfcm.modified(stoc , list(Fmax = Fmsy , Btrig = mean(c(bpa,blim))) )
    # apply the HCR
  }
  

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

# same as above expect that reduction in F starts from 1  as soon as B<Bpa 
HCR.gfcm.modified <- function(stoc, target)
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
  sliding.slope[sliding.slope>1] <- 1                  
  
  Ftarget <- Ftarget * sliding.slope
  
  targ<-list()
  targ$quant = "f"
  targ$val = list(y1 = c(fbar(stoc)[,ac(range(stoc)["maxyear"])]) , y2 = Ftarget)
  targ$rel = c(NA,NA)
  return(targ)
} 


# to first go  to Bpa  in a givne target year and then apply the HCR GFCM
HCR.Bpa.Fmsy<- function(stoc, target)
{
  # define time and F 
  ystart  <- ay+1                      # first simulation year
  ylast   <- range(stoc)["maxyear"]    # last assessment year in the current loop
  blast   <- ssb(stoc)[,ac(ylast)] 
  ynow    <- ylast + 1                 # current simulation year
  ytarg   <- target[["Yr.targ"]]       # time target to be at Ftarget
  flast   <- fbar(stoc)[,ac(ylast)]    # F in the last assessment year in the current loop
  ftarget <- target$Ftarget 
  btarget <- target$Btarget 
  
  # first target is Btarget in the first year
  if(ynow  <ytarg) {
    
    # find the Fbar which brings SSB at Bpa
    targ<-list()
    targ$quant = c("f","ssb")
    B.annual.change <- (btarget - blast) / (ytarg - ynow)
    btargetnow <-  blast + 2*B.annual.change
    targ$val = list(y1 = c(fbar(stoc)[,ac(range(stoc)["maxyear"])]) , y2 = c(btargetnow))
    targ$rel = c(NA,NA)
  } else {               # apply Fmsy
    
    targ<-HCR.gfcm.modified(stoc , list(Fmax = Fmsy , Btrig = mean(c(bpa,blim))) )
    # apply the HCR
  }
  
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
              target = list(Ftarget = mean(c(fbar(stk)[,ac((dy-2):dy)]))) ,
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
#----- scenarios with a catch ban ------------------------

# Catch ban followed by TAC
S5a <- list( name = "S5a" ,
               target = list(Ctarget = ifelse(range(stk)["maxyear"] %in% 2017:2022, 0, CTAC)),
                HCR =  HCR.cstC,
                spatial.closure = F ,
                additionnal.F.reduction = NA
              ) 

# Catch ban followed by Fmsy
S5b <- list( name = "S5b" ,
             target = list(Ftarget = 0.01) ,
             HCR =  HCR.cstF ,
             spatial.closure = F ,
             additionnal.F.reduction = NA
)

################################################################################
#----- Constant catch scenarios ------------------------

# scenario Catch at TAC + 50%IUU
CTACIUU50 <- list( name = "CTACIUU50" ,
               target = list(Ctarget = CTAC50) ,
               HCR =  HCR.cstC ,
               spatial.closure = F ,
               additionnal.F.reduction = NA
)    

# scenario Catch at TAC
CTACNOIUU <- list( name = "CTACNOIUU" ,
                  target = list(Ctarget = CTAC) ,
                  HCR =  HCR.cstC ,
                  spatial.closure = F ,
                  additionnal.F.reduction = NA
)    

################################################################################
#----- 	Bpa objective for 2020 and Fmsy for 2020  ------------------------

Bpa.Fmsy2020 <- list( name= "Bpa.Fmsy2020",
                      target =  list(Yr.targ = 2020 , Ftarget = Fmsy , Btarget = bpa),
                      HCR = HCR.Bpa.Fmsy ,
                      spatial.closure = F ,
                      additionnal.F.reduction = NA
) 

########################################
# combine them in a list

management.scenarios <- list(F.sq, F.msy, CTACNOIUU, CTACIUU50,     
                              Bpa.Fmsy2020, S5a, S5b)

names(management.scenarios) <- lapply (management.scenarios , function(x) x[[1]])