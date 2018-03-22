###############################################################################
#
#  look at SE
#
##############################################################################

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
#       sam.ctrl.new <- FLSAM.control(stk,ids)
#       for (slt in slotNames(sam.ctrl.new))  try( slot(sam.ctrl.new,slt) <- slot(sam.ctrl,slt))

# try with a 0% mature at age 0 
 stk2<-stk
 stk2@mat[1] <- 0
  ANCHOVY2.sam <- FLSAM(stk2,ids,sam.ctrl)
   rbind(ssb(ANCHOVY2.sam)[,2],ssb(ANCHOVY.sam)[,2])
stk2<-stk2+ANCHOVY2.sam


plot(FLStocks(stk,stk2))
#==============================================================================
# Single species MSE
#==============================================================================

#==============================================================================
# Stochastic projections to show example of envelope analysis
#------------------------------------------------------------------------------
# Frange: 0.227-0.364
# Btrig: 194000
# Bpa: 194000 using 138500 
# Blim: 138500 * 0.5
# Fmsy: 0.3 
#==============================================================================


plot(ssb(stk2),rec(stk2),xlim=c(0,160000) , ylim= c(0,2e8))
#------------------------------------------------------------------------------

# fit hockey stick
# THIS PART AS BEEN MODIFIED SO THAT A MODEL IS FITTED FOR EACH ITERATION
# THERE ARE NOW AS MANY PARAMETERS AS ITERATIONS
sr <- fmle(as.FLSR(sstk, model="segreg"), fixed=list(b=mean(ssb(stk),na.rm=T)))
                             # this fits the model for each iteration, 
                             # but I don't know how to use the mean(ssb) specific to each iteration
                             # using yearMeans(ssb(sstk) does not work               
                             # so this line is just to create a FLPar object of the right dimension
# now do the actual parameter estimation per iteration
for (its in 1:it)  iter(params(sr),its) <- params(fmle(as.FLSR(iter(sstk,its), model="segreg"), fixed=list(b=mean(iter(ssb(sstk),its),na.rm=T))))


# calculate residuals of the fit and look
plot(sr)
sr.res <- residuals(sr)
plot(sr.res)


# THIS IS MODIFIED SO THAT AN ARIMA MODEL IS FITTED FOR THE RESIDUALS OF EACH ITERATION
# AND USE TO PRODUCE THE FUTURE DEVIATIONS FOR THE CORRESPONDING ITERATION
# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
### S/R residuals - with autocorrelation
rec.res <- residuals(sr)

# autoregressive model order 1
set.seed(108)
# a list with one model per iteration
arima.fit.lst <- lapply(as.list(dimnames(rec(sr))$iter) ,  function(its) {arima(an(iter(rec.res,its)), order = c(1, 0, 0))})

# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.res <- make.arma.resid.lst(arima.fit.lst, age = 0, years = iy:(iy + ny-1))
plot(sr.res)




# WITHOUT AUTOCORRELATION:
##sr.res1[] <- sample(c(residuals(sr)), ny*it, replace=TRUE)

# PLOTS

# Confidence interval plot
ssb <- apply(matrix(as.vector(stk@stock.n*stk@mat*stk@stock.wt), nrow=1+stk@range[[2]] ), 2, sum)       ### this does seem a little strange? why computing SSB this way (and not taking spawning time into account
ssb <- c(ssb(stk)@.Data)

segreg.meanssb  <- function(ab, ssb) log(ifelse(ssb >= mean(ssb), ab$a*mean(ssb), ab$a*ssb))
fit.meanssb <- eqsr_fit(stk,nsamp=2000, models = c("segreg.meanssb"))

png("Results/ANE/PLOTS/ANE_segreg_bkptmeanSSB_confidence.png", width=700, height=700)
eqsr_plot(fit.meanssb)
dev.off()

# GENERAL PLOTS
# Labeled plot (with years)
fitted <- as.data.frame(fitted(sr)); names(fitted)[7] <- "fitted"
Ssb    <- as.data.frame(ssb(sr))   ; names(Ssb)[7]    <- "ssb"
Obs    <- as.data.frame(rec(sr))   ; names(Obs)[7]    <- "obs"



df.sr.bktp.meanssb <- cbind(fitted[,-1] , data.frame(ssb = Ssb[,7] ,obs = Obs[,7]))
df.sr.bktp.meanssb.med <- aggregate(cbind(fitted,ssb,obs)~year , df.sr.bktp.meanssb , median)

png("Results/ANE/PLOTS/ANE_SR_labeled.png", width=700, height=700)
ggplot(df.sr.bktp.meanssb, aes(x=ssb, y=obs , colour = factor(year))) + 
  geom_point() +
  geom_point(data = df.sr.bktp.meanssb.med ,  size = 8) +
  geom_line(data = df.sr.bktp.meanssb.med , aes(colour="black"))  +
  
  theme_light(base_size=12) + theme(legend.position = "none") + 
  xlab("SSB") + ylab("R") + ggtitle("Anchovy GSA1718") +
  theme(axis.title.x = element_text(face="bold", size=25),
        axis.title.y = element_text(face="bold", size=25),
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size=25, face="bold"))
dev.off()

png("Results/ANE/PLOTS/ANE_SSB.png", width=700, height=700)
ggplot(df.sr.bktp.meanssb, aes(x=year, y=ssb , colour = iter ,group=iter)) + 
  geom_point(size=4) +     geom_line()  +
  theme_light(base_size=12) + theme(legend.position = "none") + 
  xlab("Year") + ylab("SSB (tonnes)") + ggtitle("Anchovy GSA1718 - SSB") +
  theme(axis.title.x = element_text(face="bold", size=25),
        axis.title.y = element_text(face="bold", size=25),
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size=25, face="bold"))
dev.off()

png("Results/ANE/PLOTS/ANE_REC.png", width=700, height=700)
ggplot(df.sr.bktp.meanssb, aes(x=year, y=obs , colour = iter ,group=iter)) + 
  geom_point(size=4) +     geom_line()    +
  theme_light(base_size=12) + theme(legend.position = "none") + 
  xlab("Year") + ylab("R") + ggtitle("Anchovy GSA1718 - Recruitment") +
  theme(axis.title.x = element_text(face="bold", size=25),
        axis.title.y = element_text(face="bold", size=25),
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size=25, face="bold"))
dev.off()









#######################################################################################################################################

# Fixed objects
TAC <- FLQuant(NA, dimnames=list(TAC="all", year=vy, iter=1:it))
BB  <- FLQuant(0, dimnames=list(TAC="all", year=vy, iter=1:it))

# 
                        
# short term forecast: start with a projection of F into the future to ny (16 yrs)
# this serves as a starting point for projecting the stock
pstk <- stf(sstk, ny, 3, 3)             # harvest is average last 3 years

landings.n(pstk) <- propagate(landings.n(pstk), it)           # NOT NEEDED already the right NB ITERS
discards.n(pstk) <- propagate(discards.n(pstk), it)

#### NOTE THAT BY DOING THIS ONLY, WE HAVE NOT VARIABILITY IN THE BIOLOGY OR FISHERIES SELECTION PATTERN IN THE FUTURE !!!!!



# Prepare index object for future surveys USING THE FORMALISE OF SAM :
# for each iteration, the corresponding catchbilities and observation SD are used
# future simulated observations  (index values) are based the SAM formula  : log(Iobs a,y)  = log (Imod a,y) + epsilon, 
#                                                                             with epsilon = N(0,SdObs)
#                                                                             and  Imod a,y = Q * Na,y exp ( - Za,y  * time of survey)
#
#  future observation will be computed as   Iobs a,y =  Imod a,y * dev,idx,
#   where dev.idx is the exp(rnorm (1 , mu = 0 , sd = SdObs)
#
# to do so we need to extract from the parameters of SAM (for a given iteration) the Q to  be used (going in the slot index.q) and extract the SdOBs, 
# make random draws, and store the exp of them, the actual dev.idx, in the slot index.var

idx <- ids
# loading the parameteres specific to each iteration
load("./RESULTS/ANE/random.param.RData")  #getting object random.param

logQ        <- random.param[,colnames(random.param)=="logFpar"]
logSdObs    <- random.param[,colnames(random.param)=="logSdLogObs"]
logQSsb     <- random.param[,colnames(random.param)=="logScaleSSB"]
logSdObsSsb <- random.param[,colnames(random.param)=="logSdSSB"]

# multiplying the number of iterations in the indices and expanding the time frame
for (i in 1:length(ids)) name(ids[[i]]) <- names(ids)[i]

idx <- lapply(ids,function(x) 
      {
      n. <- name(x)
      cat(n.)      
      # multiplying the number of iterations in the indices and expanding the time frame
      x <- propagate(x,it)
      x <- window(x,start=range(x)["minyear"],end=fy)
      y0idx <- range(x)["minyear"]
      
      fleet.type <- sam.ctrl@fleets[n.]
     
      # fill in the catchability slot and index deviations for age structured surveys
      if (fleet.type == 2 )
      {
      agesQ     <-  names(sam.ctrl@catchabilities[n.,]) [!is.na(sam.ctrl@catchabilities[n.,])]
      agesSdObs <-  names(sam.ctrl@obs.vars[n.,]) [!is.na(sam.ctrl@obs.vars[n.,])]
      
      for (its in 1:it) 
          {
          # each iteration the estimated catchability in the slot index.q for each survey :
          iter( index.q(x)[agesQ,]   , its )     <- FLQuant( matrix(rep(exp(logQ[its,sam.ctrl@catchabilities[n.,agesQ]]),dim(index.q(x)[agesQ,])[2]),nrow = dim(index.q(x)[agesQ,])[1]), dimnames = dimnames(iter( index.q(x)[agesQ,] , its )))
          #  the observation SD in the slot index.var
          idx.dev<- FLQuant( matrix(rep(exp(logSdObs[its,sam.ctrl@obs.vars[n.,agesSdObs]]),dim(index.var(x)[agesSdObs,])[2]),nrow = dim(index.var(x)[agesSdObs,])[1]), dimnames = dimnames(iter( index.var(x)[agesSdObs,] , its )))
          # we want to replace that by actual deviations that we can just multiply by the modeled index
          iter( index.var(x)[agesSdObs,] , its )  <-  exp(apply(idx.dev , c(1,2) , function (x) rnorm(1,0,x)))
          # we don't want to apply a deviation on the historical data, so impose a value of 1
          index.var(x)[agesSdObs,ac(y0idx:ay),,,,its]   <- 1
          }
      }
      
      
      # fill in the catchability slot and index deviations  for SSB surveys
      if (fleet.type == 3 )
      {
      for (its in 1:it) 
          {
# each iteration the estimated catchability in the slot index.q for each survey :
          iter( index.q(x)  , its )  <- FLQuant( matrix(rep(exp(logQSsb[its]),dim(index.q(x))[2]),nrow = dim(index.q(x))[1]), dimnames = dimnames(iter( index.q(x) , its )))
#  the observation SD in the slot index.var
          idx.dev <- FLQuant( matrix(rep(exp(logSdObs[its]),dim(index.var(x))[2]),nrow = dim(index.var(x))[1]), dimnames = dimnames(iter( index.var(x) , its )))
          # we want to replace that by actual deviations that we can just multiply by the modeled index          
          iter( index.var(x) , its ) <-   exp(apply(idx.dev , c(1,2) , function (x) rnorm(1,0,x)))
          # we don't want to apply a deviation on the historical data, so impose a value of 1
          index.var(x)[,ac(y0idx:ay),,,,its]   <- 1
          
############ this is how it should be done if the survey was to be continued in the future, but apparently, this SSB survey stopped after 2012,
# so replace all the above by NAs          
          index.var(x)[,ac(iy:fy),,,,its]   <- NA  # so that we generate NAs in the future ,
                                                   # we could just not include this survey in the part where new index values are generate (loops),
                                                   # but in case one day this survey constinue, maybe it is better to keep it in the framework.         
          }
      }    
       
    return(x)  
      })





## similarly generate a future errors for the catch matrix, deviations to be use as multipliers
      agesSdObs <-  names(sam.ctrl@obs.vars["catch",]) [!is.na(sam.ctrl@obs.vars["catch",])]
      catch.dev <- window(catch.n(sstk) , start=range(sstk)["minyear"],end=fy)

      for (its in 1:it) 
        {
        cdev <- matrix(rep(logSdObs[its,sam.ctrl@obs.vars["catch",agesSdObs]],dim(catch.dev)[2]),nrow = dim(catch.dev)[1])
        cdev <- exp(cdev)
        cdev <- FLQuant(cdev, dimnames = dimnames(iter( catch.dev , its )))
        cdev <- apply(cdev , c(1,2) , function (x) rnorm(1,0,x))
        cdev[,ac(y0:ay)]    <- 0
        cdev                <- exp(cdev) # just transform in the dimension in which it can be used as a multiplier on the catches
        iter(catch.dev,its) <- cdev
        }





#------------------------------------------------------------------------------
# 1a. (1a) STATUS QUO SCENARIO - SEGMENTED STOCK RECRUITMENT WITH BPT AT MEAN SSB
#------------------------------------------------------------------------------

# NEED AN ASSUMPTION   TO FORECAST THE OM IN THE FIRST YEAR
TAC[,(ac(iy))] <- 15000 #  WHAT IS A REALISTIC VALUE ?
# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- 45936
bpa <- blim*2
Btrig <- blim+((bpa-blim)/2)
dt <- date()

#########################################################
# go fish!

for(i in vy[-length(vy)]){   #a[-(15:16)]
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  iay <- an(i)   # an is equivalent to as.numeric
  cat(i, ">")
  vy0 <- 1:(iay-y0) # data years (positions vector)
  sqy <- (iay-y0-nsqy+1):(iay-y0) # status quo years (positions vector)
  #sqy <- (iay-y0-nsqy+1):(iay-y0)
  # define stock0 from pstk until the last populated year
  # pstk is at the beginning only populated into the future for F
  # the rest is only 1975-2016 but as the loop progresses through the projection
  # years the object is populated with projected numbers
  
  # CREATE PERCEIVED STOCK
  stk0 <- pstk[,vy0]
  # APPLY UNCERTAINTY ON CATCHES (THE WHOLE ARRAY IS UPDATED, BUT SINCE THE CATCH.DEV ARE PREDEFINED, THE OBSERVED CATCH MATRIX DOES NOT CHANGE AS THE SIMULATION PROGRESSES
  # NOTE THAT THE HISTORIC PART OF THE CATCH MATRIX IS NOT AFFECTED
  # add 1 to everything to avoid zeros
  catch.n(stk0) <- catch.n(stk0) * catch.dev[,vy0] # avoid zeros
  
  
  
  # CREATE PERCEIVED TUNNING INDICES
  idx0 <- window(idx , end = iay-1) 
  # COMPUTE SURVEY INDICES BASED ON THE ESTIMTAED CATCHABILITY AND ADD UNCERTAINTY
  # NOTE THAT HISTORICAL PART OF THE TIME SERIES ARE NOT MODIFIED
  if (iay > iy) # don't do that for the first year because it correspond to the actual year the current assessment was carried out
  {
      for (idx.nb in 1:length(idx))
      {
        n. <- name(idx[[idx.nb]])
        # define necessary dimensions
        idsty<-range(idx[[idx.nb]])["minyear"]                            # first year of the survey
        idx0[[idx.nb]] <- idx[[idx.nb]][,ac(idsty:(iay-1))]
        type.idx <- sam.ctrl@fleets[n.]                                   # type of index (abund or ssb)
        age.surv <- dimnames(index(idx0[[idx.nb]]))$age                   # age range
        time.surv <- mean(c(range(idx0[[idx.nb]])[c("startf","endf")]))   # time of the year when survey is carried out
        
        if (type.idx == 2) 
        {
        mod.n <- stock.n(stk0) * exp (- time.surv * ( harvest(stk0) + m(stk0)))
        index(idx0[[idx.nb]])[age.surv,ac(iy:(iay-1))] <- mod.n[age.surv,ac(iy:(iay-1))]*index.q(idx[[idx.nb]])[age.surv,ac(iy:(iay-1))]     # modelled index
        index(idx0[[idx.nb]])[age.surv,ac(iy:(iay-1))] <- index(idx0[[idx.nb]])[age.surv,ac(iy:(iay-1))] * index.var(idx[[idx.nb]])[age.surv,ac(iy:(iay-1))]   # adding uncertainty
        }
        
        if (type.idx == 3) 
        {
        index(idx0[[idx.nb]])[,ac(iy:(iay-1))] <- ssb(pstk)[,ac(iy:(iay-1))]*index.q(idx[[idx.nb]])[,ac(iy:(iay-1))]           # SSB index    
        index(idx0[[idx.nb]])[,ac(iy:(iay-1))] <- index(idx0[[idx.nb]])[,ac(iy:(iay-1))] *  index.var(idx[[idx.nb]])[,ac(iy:(iay-1))]    # adding uncertainty
        }             
      }
  }
  ##
 
#### DO THE ASSESSMENT FOR EACH ITERATION SUCCESIVELY  AND UPDATE THE PERCIEVED STOCK 
  sam0.ctrl <- sam.ctrl
  sam0.ctrl@nohess <- T
  sam0.ctrl@range["maxyear"]  <- iay-1

  for (its in 1:it)
    {
    idx0its <- FLIndices(lapply(idx0 , function(x) iter(x,its)))
    sam0its <- FLSAM(iter(stk0,its), idx0its , sam0.ctrl)
    iter(stk0,its) <- iter(stk0,its) + sam0its
    }
  
  
  
### STF ON THE PERCIEVED STOCK TO PRODUCE AND ADVICE
  # fwd control
  # what is F status quo? is it fixed or does it vary as you progress in the projection?
  fsq0 <- fsq # status quo 2013-2015 from SAM (deterministic)
  # dnms <- list(iter=1:it, year=c(iay, iay + 1), c("min", "val", "max"))   # I changed that because the order of the dimensions did not correspond to those in the crtl object
                                                                            # this was wrong and did not project the stock correctly
                                                                            # for instance when doing fbar(stkTmp), we did not get fsq0
  
  dnms <- list(year=c(iay, iay + 1), c("min", "val", "max"),iter=1:it) 
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(harvest = 1)
  ftrg.vec <- an(fsq0) # Ftarget = status quo
  #Bescape <- blim
  arr0[,"val",] <- c(fsq0, ftrg.vec)                                        # changed as above
  #arr0[,,"min"] <- c(rep(NA, 2 * it), rep(Bescape, it))
  #arr0 <- aperm(arr0, c(2,3,1))
  # in Control you define what you want to vary in iay and iay+1 (which is F)
  ctrl <- fwdControl(data.frame(year=c(iay, iay+1), quantity=c('f', 'f'), val=c(fsq0, ftrg.vec)))
  ctrl@trgtArray <- arr0
  ## Short term forecast of stk0
  stkTmp <- stf(stk0, 2)
  # project forward with the control you want and the SR rel you defined above, with residuals
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr) #, sr.residuals = exp(sr.res[,ac(iay:(iay+1))]), sr.residuals.mult = TRUE) #  !!!!!! There should not be any residuals here
  TAC[,ac(iay+1)] <- catch(stkTmp)[,ac(iay+1)]



### UPDATE THE OM BASED ON THE TAC ADVICE (with on year lag)
  # OM proj
 ctrl@target <- ctrl@target[2,]
 ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]   ### I don't  agree with this... we cannot apply the Ftarget to the OM because the realised F is different due to assessment errors

# instead we have to find the F to applied, based on the TAC, assuming catch = TAC  
 dnms <- list(year=c(iay), c("min", "val", "max"),iter=1:it)
 arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
 arr0[,"val",] <- TAC[,ac(iay)]
 ctrlOM <- fwdControl(data.frame(year=c(iay), quantity=c('catch'), val=c(iterMeans(TAC[,ac(iay)])@.Data)))
 ctrlOM@trgtArray <- arr0
 # update pstk with stkTmp
 pstk <- fwd(pstk, ctrl=ctrlOM, sr=sr, sr.residuals = exp(sr.res[,ac(iay)]), sr.residuals.mult = TRUE) #
}


return(val)
date()

# breakpoint mean(SSB)
ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0 <- pstk
ANE_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0 <- stk0
plot(ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
save(list = c("ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0", "ANE_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0"), file = "Results/ANE/ANE_Opt1a_mse_SegregmeanSSB_statusQuo_250it.RData")
png("ANE_Opt1a_mse_SegregmeanSSB_statusQuo_250it.png", width=700, height=700)
plot(ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
dev.off()




