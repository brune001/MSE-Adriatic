###############################################################################
# EJ(20150119)
# Tests on the evaluation of the NS-MAP using COD
# NOTE1: The final analysis is in the report file.
#  One may tangle the Snw file to get the R script trimmed of comments.
#
# Modified from Piera Carpi, AgurtzPIL Urtizberea Ijurco & Miguel Bernal
# 20 January 2016
# Changed for SARDINE MSE for GFCM 
#
# Modified by Betulla Morello
# February 2017
# Updated for 2017 GFCM WKMSE, using WGSAD2016 assessment data (last year=dy=2015)
###############################################################################

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

# source needed functions 
setwd("N:/Projecten/MSEadriatic/Work")
source('./Code/MSE_funs_LAST.R')

#==============================================================================
# Read data
# Inputs and outputs of the Sardine SAM accepted at 2016 GFCM WGSAD 
#==============================================================================

load("./Data/PIL stock and SAM.RData")

#==============================================================================
# Fit a4a model to replicate SAM as much as possible
#==============================================================================

# my stock
stk <- stk2 
# my tuning index
ids <- idx1
ids <- FLIndices(ids)

# Replacing zeros with small numbers in the tuning index
ids[[1]]@index[ids[[1]]@index == 0] <- 0.1
ids[[1]]@catch.n[ids[[1]]@catch.n == 0] <- 0.1
ids[[2]]@index[ids[[2]]@index == 0] <- 0.1
ids[[2]]@catch.n[ids[[2]]@catch.n == 0] <- 0.1
# Assign names to tuning indices
names(ids) <- c("EchoWest","EchoEast")

#m.spwn(stk) <- 0
#harvest.spwn(stk) <- 0

#####################################################
# Run FLa4a to replicate SAM accepted assessment
#####################################################

# set single GAMs for the things you think are important: q, F and R

#####
# Catchability model
qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4)) # final best model: depending on age (with 5 nodes) and year (with 4 nodes) for each survey

#####
#F model
mod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5))    # good! Final best model

#####
# Recruitment model 
rmodel4 <- ~s(year, k=20)                                               # Final best model


#### 
# Fit a4a
fit <- a4aSCA(stock = stk, srmodel=rmodel4, fmodel = fmod8, qmodel=qmod5, indices = ids, verbose = FALSE, fit = "assessment")

# Create new object with new results
fitA4A <- stk+fit


# We now have an operational model!

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

#------------------------------------------------------------------------------
# variables to fix and OM conditioning
#------------------------------------------------------------------------------

it <- 250                                # iterations - should be 250
y0 <- range(stk2)["minyear"]     # year zero (initial) = 1975
ny <- 16                                # number of years to project - Usually 20
# In order for this code to run iy = dy
dy <- 2015                              # data year
ay <- 2015                              # assessment year
iy <- 2015                              # initial projections year (also intermediate)
fy <- iy + ny -1                        # final year
vy <- ac(iy:fy)
nsqy <- 3                               # number of SQ years upon which to average results

mny <- 2020                             #2016 # min year to get to trg
mxy <- 2020                             # 2016 # max year to get to trg


# Management quantities
#flo <- 0.23
#fup <- 0.36
#fmsy <- 0.55
# 1. F status quo: maintain F from 2015
fsq <- mean(c(fbar(stk)[,ac(dy)]))
# 2. Blim from 2015 benchmark reference points (tonnes)
blim <- 125318
# 3. Bpa from 2015 benchmark reference points (tonnes)
bpa <- 2*blim

#------------------------------------------------------------------------------
# S/R
#------------------------------------------------------------------------------

#############################################
# 1. hockey stick with break point at Bpa
#############################################
# fit hockey stick
sr.bkpt.bpa <- fmle(as.FLSR(stk, model="segreg"), fixed=list(b=bpa))# method="L-BFGS-B"mean(ssb(stk))
# calculate residuals of the fit and look
sr.res.bkpt.bpa <- residuals(sr.bkpt.bpa)
plot(sr.res.bkpt.bpa)

# For the autocorrelation issue with segreg model
# get a and b from sr model
a.bkpt.bpa <- as.numeric(sr.bkpt.bpa@params["a"]) #  <- a <- 713.1137
b.bkpt.bpa <- as.numeric(sr.bkpt.bpa@params["b"]) # <- b <- bpa

## S/R residuals - with autocorrelation
rec.res.bkpt.bpa <- residuals(sr.bkpt.bpa)
plot(rec.res.bkpt.bpa)

# autoregressive model order 1
set.seed(108)
# mean
arima.fit.bkpt.bpa <- arima(an(rec.res.bkpt.bpa), order = c(1, 0, 0))
# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.res.bkpt.bpa <- make.arma.resid(arima.fit.bkpt.bpa, age = 0, years = iy:(iy + ny-1), nit = it)
plot(sr.res.bkpt.bpa)

#############################################
# 2. hockey stick with break point at meanSSB over the entire time series:note that the mean(ssb(stk) is 362979.7; very different to bpa
#############################################
# fit hockey stick
sr <- fmle(as.FLSR(stk, model="segreg"), fixed=list(b=mean(ssb(stk))))# method="L-BFGS-B"mean(ssb(stk))

# calculate residuals of the fit and look
sr.res <- residuals(sr)
plot(sr.res)

# For the autocorrelation issue with segreg model
# get a and b from sr model
a <- as.numeric(sr@params["a"]) #  <- a <- 713.1137
b <- as.numeric(sr@params["b"]) # <- b <- mean(ssb(stk))

## S/R residuals - with autocorrelation
rec.res <- residuals(sr)

# autoregressive model order 1
set.seed(108)
# mean
arima.fit <- arima(an(rec.res), order = c(1, 0, 0))
# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.res <- make.arma.resid(arima.fit, age = 0, years = iy:(iy + ny-1), nit = it)




# Fixed objects
TAC <- FLQuant(NA, dimnames=list(TAC="all", year=vy, iter=1:it))
BB <- FLQuant(0, dimnames=list(TAC="all", year=vy, iter=1:it))

# # Prepare stock objects we need, with iterations and propagate towards final year
stk <- iter(SARDINE2idx, 1)
# simulate new stock ased on a4a final fit
sstk <- stk + simulate(fit, it)
summary(sstk)
plot(sstk)
# short term forecast: start with a projection of F into the future to ny (16 yrs)
pstk <- stf(sstk, ny, 3, 3)             # harvest is average last 3 years

landings.n(pstk) <- propagate(landings.n(pstk), it)
discards.n(pstk) <- propagate(discards.n(pstk), it)


# Prepare index object 

idx <- ids

# create an index of abundance using the error associated to the real index
for (i in 1:length(idx)){
  lst <- mcf(list(idx[[i]]@index, stock.n(stk)))
  idx.lq <- log(lst[[1]]/lst[[2]])
  idx.lq[is.infinite(idx.lq)] <- NA # fix zeros
  idx.qmu <- idx.qsig <- stock.n(iter(pstk,1))
  idx.qmu[] <- yearMeans(idx.lq)
  idx.qsig[] <- log((sqrt(yearVars(idx.lq))/yearMeans(idx.lq))^2 + 1)
  idx.q <- idx_temp <- FLQuant(NA, dimnames=dimnames(stock.n(pstk)))
  idx.q[,ac(y0:dy)] <- propagate(exp(idx.lq[,ac(y0:dy)]), it)
  idx.q[!is.na(idx.qmu)] <- rlnorm(it, idx.qmu[!is.na(idx.qmu)], idx.qsig[!is.na(idx.qmu)])
  plot(idx.q)
  idx_temp <- idx.q * stock.n(pstk)
  idx[[i]] <- FLIndex(index=idx_temp, index.q=idx.q)
  range(idx[[i]])[c("startf", "endf")] <- c(0, 0)
  plot(index(idx[[i]]))
}

#------------------------------------------------------------------------------
# 1a. STATUS QUO SCENARIO - SEGMENTED STOCK RECRUITMENT WITH BPT AT MEAN SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- min(ssb(stk))
bpa <- blim*2
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
dt <- date()

#########################################################
# go fish!

for(i in vy[-length(vy)]){   #a[-(15:16)]
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  ay <- an(i)   # an is equivalent to as.numeric
  cat(i, ">")
  vy0 <- 1:(ay-y0) # data years (positions vector)
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
  #sqy <- (ay-y0-nsqy+1):(ay-y0)
  # define stock0 from pstk until the last populated year
  # pstk is at the beginning only populated into the future fro F
  # the rest is only 1975-2015 but as the loop progresses through the projection
  # years the object is populated with projected numbers
  stk0 <- pstk[,vy0]
  # add 1 to everything to avoid zeros
  catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
  ## note that vy0 is changing below so index is being updated
  for (index_counter in 1:length(idx)){
    idx0[[index_counter]] <- idx[[index_counter]][,vy0]
    index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
  }
  ##
  qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
  fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
  rmodel4 <- ~ s(year, k=20)  
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5)
  stk0 <- stk0 + fit
  # fwd control
  # what is F status quo? is it fixed or does it vary as you progress in the projection?
  fsq0 <- fsq # status quo 2013-2015 from SAM (deterministic)
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(harvest = 1)
  ftrg.vec <- an(fsq0) # Ftarget = status quo
  #Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  #arr0[,,"min"] <- c(rep(NA, 2 * it), rep(Bescape, it))
  #arr0 <- aperm(arr0, c(2,3,1))
  # in Control you define what you want to vary in ay and ay+1 (which is F)
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
  ## Short term forecast of stk0
  stkTmp <- stf(stk0, 2)
  # project forward with the control you want and the SR rel you defined above, with residuals
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  # update pstk with stkTmp
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
}


return(val)
date()

# breakpoint mean(SSB)
Pil_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0 <- pstk
Pil_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0 <- stk0
plot(Pil_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
save(list = c("Pil_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0", "Pil_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0"), file = "Pil_Opt1a_mse_SegregmeanSSB_statusQuo_250it.RData")
png("PIL_Opt1a_mse_SegregmeanSSB_statusQuo_250it.png", width=700, height=700)
plot(Pil_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
dev.off()



#------------------------------------------------------------------------------
# 4. SCENARIO BASED BOTH ON F AND BIO - CHOOSE THE MORE PRECAUTIONARY OPTION
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
Fmsy <- 0.715
Ftar <- Fmsy
dt <- date()
trgyF <- 2020
trgyB <- 2018
Fmin <- 0.1
#########################################################
# go fish
for(i in vy[-length(vy)]){
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  ay <- an(i)   # an is equivalent to as.numeric
  cat(i, "\n")
  vy0 <- 1:(ay-y0) # data years (positions vector)
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
  stk0 <- pstk[,vy0]
  catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
  ## note that vy0 is changing below so index is being updated
  for (index_counter in 1:length(idx)){
    idx0[[index_counter]] <- idx[[index_counter]][,vy0]
    index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
  }
  qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
  fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
  rmodel4 <- ~ s(year, k=20)  
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5) #, qmodel=qmod5)
  #rmodel <- ~geomean(CV = 0.4)
  #fmod <- ~ s(age, k = 5) + factor(year)
  #fmod <- ~ s(age, k=5) + s(year, k = 14) + s(age, year, k=30)  # Good (still high in the past, 
  #fmod <- ~s(age, k = 5) + s(year, k=8) + te(age, year, k = c(4,5))
  #fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel, qmodel=qmod)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1, ay + 1), c("min", "val", "max"))
  
  ftrg.q <- hcr.Ftar.Btar(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = data.frame(ssb = 1, harvest = 1), Btrig = Btrig, Fmin = 0.1, Blim = blim, trgyF = trgyF)
  ftrg.vec <- an(ftrg.q)
  
  SSBtrg0 <- ssb(stk0)[, ac(an(i) - 1)] - (ssb(stk0)[, ac(an(i) - 1)]-Btrig)/ifelse(trgyB - ay < 1, 1, trgyB - ay)
  SSBtrg.vec <- an(SSBtrg0)
  
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- c(fsq0, ftrg.vec, rep(NA, it))
  arr0[,,"min"] <- c(rep(NA, 2 * it), SSBtrg.vec)
  arr0 <- aperm(arr0, c(2,3,1))
  
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1, ay + 1), quantity=c('f', 'f', 'ssb'), val=NA))
  ctrl@trgtArray <- arr0
  
  stkTmp <- stf(stk0, 3)
  stkTmp10 <- fwd(stkTmp[,,,,,], ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #,exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2:3,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2:3,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}


##
return(val)
date()
Pil_Opt4_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it <- pstk
Pil_Opt4_mse.stk0.GFCM_segreg_BtarFtar_trgy_250it <- stk0
save(list = c("Pil_Opt4_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it", "Pil_Opt4_mse.stk0.GFCM_segreg_BtarFtar_trgy_250it"), file = "PIL_Opt4_mse_Segreg_FtarBtar_250it.RData")
png("Pil_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it.png", width=700, height=700)
plot(window(Pil_Opt4_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it, end=2030)) + ggtitle("Sardine - Target Fmsy & recover to Bpa") + theme_light(base_size=12)
dev.off()

iterMeans(ssb(Pil_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it) < blim)


#------------------------------------------------------------------------------
# 6a. ENFORCEMENT OF FISHING DAYS REDUCTION - REDUCTION OF 9% OF F IN 2015
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
dt <- date()

#########################################################
# go fish
for(i in vy[-length(vy)]){
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  ay <- an(i)   # an is equivalent to as.numeric
  cat(i, "\n")
  vy0 <- 1:(ay-y0) # data years (positions vector)
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
  stk0 <- pstk[,vy0]
  catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
  ## note that vy0 is changing below so index is being updated
  for (index_counter in 1:length(idx)){
    idx0[[index_counter]] <- idx[[index_counter]][,vy0]
    index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
  }
  ##
  qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
  fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
  rmodel4 <- ~ s(year, k=20)  
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5) #, qmodel=qmod5)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,c(sqy)])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  #refpt <- data.frame(ssb = 1, harvest = 1)
  #ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, refpt = refpt, Btrig = Btrig, Fmin = 0, Blim = blim, Bpa=bpa)
  #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
  #ftrg.vec <- an(ftrg.q)
  #Bescape <- blim
  ftrg.q <- fbar(stk)[,"2014",,,,] * 0.91 
  ftrg.vec <- rep(an(ftrg.q),it)
  arr0[,,"val"] <- c(fsq0, ftrg.vec) #rep(NA, it)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
  #future_catch <- c(catch(stk0)[,"2013"]) * 0.9
  #ctrl_catch <- fwdControl(data.frame(year=an(ay:(ay+1)), quantity = "catch", val=future_catch))
  #ctrl_target <- ctrl_target[order(ctrl_target$year),]
  #ctrl <- fwdControl(ctrl_catch)
  #ctrl <- fwdControl(data.frame(year=c(ay, ay+1, ay + 3), quantity=c('f', 'f', 'ssb'), val=NA))
  #ctrl@trgtArray <- arr0
  ## 
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}



##

date()
PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it <- pstk
PIL_Opt6a_mse.stk0.GFCM_segreg_9perc_red_250it <- stk0
#plot(PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it)
save(list = c("PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it", "PIL_Opt6a_mse.stk0.GFCM_segreg_9perc_red_250it"), file = "PIL_Opt6a_mse_segreg_9perc_red_250it.RData")
# png("PIL_Opt6a_mse_Segreg__9perc_red_250it_250it.png", width=700, height=700)
# plot(PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it)
# dev.off()



