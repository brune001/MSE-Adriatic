###############################################################################
# EJ(20150119)
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
library(FLCore)
library(FLEDA)
library(FLXSA)
library(SQLiteFL)
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

load("./Data/ANE/Anchovy GSA 17-18 (1).RData")


stk                           <- ANCHOVYnew
sam                           <- ANCHOVY.sam
ids                           <- ANCHOVY.tun
sam.ctrl                      <- ANCHOVY.ctrl

# just a quick check that we can reproduce the assessment
ANCHOVY2.sam <- FLSAM(ANCHOVYnew,ANCHOVY.tun,ANCHOVY.ctrl)
rbind(ssb(ANCHOVY2.sam)[,2],ssb(ANCHOVY.sam)[,2])


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

#-------------------------------------------------------------------------------
# 1 : setting dimensions
#-------------------------------------------------------------------------------

it <- 25     # 250                      # iterations - should be 250
y0 <- range(ANCHOVYnew)["minyear"]      # year zero (initial) = 1975
ny <- 20                                # number of years to project - Usually 20
# In order for this code to run iy = dy
dy <- range(ANCHOVYnew)["maxyear"]      # data year
ay <- dy                                # assessment year
iy <- ay+1                              # initial projections year (also intermediate)
fy <- iy + ny -1                        # final year
py <- y0:dy                             # past years
vy <- ac(iy:fy)                         # future years
nsqy <- 3                               # number of SQ years upon which to average results

mny <- 2020                             #2016 # min year to get to trg
mxy <- 2020                             # 2016 # max year to get to trg

#-------------------------------------------------------------------------------
# 2 : Create stock object & use vcov for new realisations
#-------------------------------------------------------------------------------

outdir <- paste0(getwd(),"/Results/ANE/")    # location where the model parameters for the nits replicates will be stored

sstk                            <- monteCarloStock2 ( stk , sam , it , run.dir=outdir,seed_number=floor(pi*10000))
#sstk                            <- window(sstk, start = y0 , end = fy)
#sstk@catch.n                    <- sstk@stock.n * sstk@harvest / (sstk@harvest + sstk@m) * (1 - exp(-sstk@harvest - sstk@m))
#sstk@landings.n                 <- sstk@catch.n





# Management quantities
#flo <- 0.23
#fup <- 0.36
#fmsy <- 0.55
# 1. F status quo: maintain F from 2015
fsq <- mean(c(fbar(stk)[,ac(dy)]))                                                   ###### why taking the mean across iterations and not the median??
# 2. Blim from 2015 benchmark reference points (tonnes)
blim <- 45936
# 3. Bpa from 2015 benchmark reference points (tonnes)
bpa <- 2*blim

#--------------------------------------------------------------------------------------------
# 4 :  S/R
#--------------------------------------------------------------------------------------------
# Tested different options (Segreg model, Ricker model, Random rec (??))
# the problem with the SR relationship is that no functional form fits: 
# there seems to be a 1:1 relationship
# demonstrated by the fact that estimated rec and ssb trends are very similar
# this is probably due to an autocorrelation of recruitment 
# (i.e. rec in one yr is correlated to rec in the previous year)
# or due to external factors (e.g. the interaction between the 2 spp and/or environment)
# so we tried 4 different things: 
# A. hockey stick with break point at meanSSB over the entire time series
# B. hockey stick with break point at maxSSB over the entire time series
# C. hockey stick with break point at meanSSB over a reduced time series (1985-2015)
# D. hockey stick with break point at maxSSB over a reduced time series (1985-2015)
# E. hockey stick with break point at Bpa
# F. geometric mean

############################################################################
#      4A. hockey stick with break point at meanSSB over the entire time series:
############################################################################
# fit hockey stick

sr <- fmle(as.FLSR(stk, model="segreg"), fixed=list(b=mean(ssb(stk),na.rm=T)))# method="L-BFGS-B"mean(ssb(stk))
#sr <- fmle(as.FLSR(stk, model="ricker")) # You could also try rickerAR1
#sr <- fmle(as.FLSR(stk, model="shepherd")) 


# NEW VERSION WITH A SR MODEL PER ITERATION
#          SR <- rep(list(as.FLSR(stk, model="segreg")),it)
#          SR <- lapply( SR ,  function(x) fmle(x, fixed=list(b=mean(ssb(stk))))
#
# calculate residuals of the fit and look
plot(sr)
sr.res <- residuals(sr)
plot(sr.res)

# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
# For the autocorrelation issue with segreg model
# get a and b from sr model
a <- as.numeric(sr@params["a"]) #  <- a <- 713.1137
b <- as.numeric(sr@params["b"]) # <- b <- mean(ssb(stk))

## S/R residuals - with autocorrelation
#ssb <- an(ssb(sr))
#pred.r <- ifelse(c(ssb) <= b, a * c(ssb), a * b)
#rec.res <- log(rec(sr)) - log(pred.r)
##plot(rec.res)
rec.res <- residuals(sr)

# autoregressive model order 1
set.seed(108)
# mean
arima.fit <- arima(an(rec.res), order = c(1, 0, 0))
# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.res <- make.arma.resid(arima.fit, age = 0, years = iy:(iy + ny-1), nit = it)
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
df.sr.bktp.meanssb <- as.data.frame(cbind(fit=as.vector(fitted(sr)), 
                                              ssb=as.vector(ssb(sr)), 
                                              obs=as.vector(rec(sr)), 
                                              year=as.vector(py)))

png("Results/ANE/PLOTS/ANE_SR_labeled.png", width=700, height=700)
ggplot(df.sr.bktp.meanssb, aes(x=ssb, y=obs)) + 
  geom_point() +
  geom_text(aes(label=year), size=8) +
  theme_light(base_size=12) + theme(legend.position = "none") + 
  xlab("SSB") + ylab("R") + ggtitle("Anchovy GSA1718") +
  theme(axis.title.x = element_text(face="bold", size=25),
        axis.title.y = element_text(face="bold", size=25),
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size=25, face="bold"))
dev.off()

png("Results/ANE/PLOTS/ANE_SSB.png", width=700, height=700)
ggplot(df.sr.bktp.meanssb, aes(x=year, y=ssb)) + 
  geom_point(size=4) +
  theme_light(base_size=12) + theme(legend.position = "none") + 
  xlab("Year") + ylab("SSB (tonnes)") + ggtitle("Anchovy GSA1718 - SSB") +
  theme(axis.title.x = element_text(face="bold", size=25),
        axis.title.y = element_text(face="bold", size=25),
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size=25, face="bold"))
dev.off()

png("Results/ANE/PLOTS/ANE_REC.png", width=700, height=700)
ggplot(df.sr.bktp.meanssb, aes(x=year, y=obs)) + 
  geom_point(size=4) +
  theme_light(base_size=12) + theme(legend.position = "none") + 
  xlab("Year") + ylab("R") + ggtitle("Anchovy GSA1718 - Recruitment") +
  theme(axis.title.x = element_text(face="bold", size=25),
        axis.title.y = element_text(face="bold", size=25),
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size=25, face="bold"))
dev.off()








skip.other.SR <- T

if (!skip.other.SR)
{
#############################################
# B. hockey stick with break point at maxSSB over the entire time series:note that the mean(ssb(stk) is 362979.7; very different to bpa
#############################################
# fit hockey stick
sr.all.max <- fmle(as.FLSR(stk, model="segreg"), fixed=list(b=max(ssb(stk))))# method="L-BFGS-B"mean(ssb(stk))
#sr <- fmle(as.FLSR(stk, model="ricker")) # You could also try rickerAR1
#sr <- fmle(as.FLSR(stk, model="shepherd")) 

# calculate residuals of the fit and look
sr.all.max.res <- residuals(sr.all.max)
plot(sr.all.max.res)

# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
# For the autocorrelation issue with segreg model
# get a and b from sr model
a <- as.numeric(sr.all.max@params["a"]) #  <- a <- 713.1137
b <- as.numeric(sr.all.max@params["b"]) # <- b <- mean(ssb(stk))

## S/R residuals - with autocorrelation
#ssb <- an(ssb(sr))
#pred.r <- ifelse(c(ssb) <= b, a * c(ssb), a * b)
#rec.res <- log(rec(sr)) - log(pred.r)
##plot(rec.res)
rec.all.max.res <- residuals(sr.all.max)

# autoregressive model order 1
set.seed(108)
# mean
arima.fit <- arima(an(rec.all.max.res), order = c(1, 0, 0))
# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.all.max.res <- make.arma.resid(arima.fit, age = 0, years = iy:(iy + ny-1), nit = it)
plot(sr.all.max.res)

# WITHOUT AUTOCORRELATION:
##sr.res1[] <- sample(c(residuals(sr)), ny*it, replace=TRUE)

##### Bootstrap residuals from SR model with geometric mean.
# OPTION 1 RANDOM RECRUITMENT INSIDE AN INTERVAL
# res.boot <- function(sr, ay) {
#   auxSR = residuals(sr)[,sample(dimnames(residuals(sr))$year, dims(sr)$year, replace=T),,,,]
#   dimnames(auxSR)$year <- ac(y0:(ay+1))
#   return(auxSR)
# }

# OPTION 2 STOCHASTIC UNCERTAINTY ON R
# multi_rec_residuals <- FLQuant(NA, dimnames = list(year=2014:2024, iter=1:it))
# residuals(sr)
# # We want to fill up our multi_rec_residuals FLQuant by randomly sampling from these log residuals
# sample_years <- sample(dimnames(residuals(sr))$year, it * 11, replace = TRUE)
# multi_rec_residuals[] <- exp(residuals(sr)[,sample_years])
# We can use this multi_rec_residuals in sr.residuals


# PLOTS

# Confidence interval plot
ssb <- apply(matrix(as.vector(stk@stock.n*stk@mat*stk@stock.wt), nrow=1+stk@range[[2]] ), 2, sum)

segreg.maxssb  <- function(ab, ssb) log(ifelse(ssb >= max(ssb), ab$a*max(ssb), ab$a*ssb))
fit.maxssb <- eqsr_fit(stk,nsamp=2000, models = c("segreg.maxssb"))

png("Results/ANE/PLOTS/ANE_segreg_bkptmaxSSB_confidence.png", width=700, height=700)
eqsr_plot(fit.maxssb)
dev.off()

##########################################################################################
# C. hockey stick with break point at meanSSB reduced (1985-2015) time series
##########################################################################################
# fit hockey stick
sr.short <- fmle(as.FLSR(trim(stk, year=1985:2015), model="segreg"), fixed=list(b=mean(ssb(stk))))# method="L-BFGS-B"mean(ssb(stk))
#sr <- fmle(as.FLSR(stk, model="ricker")) # You could also try rickerAR1
#sr <- fmle(as.FLSR(stk, model="shepherd")) 

# calculate residuals of the fit and look
sr.short.res <- residuals(sr.short)
plot(sr.short.res)

# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
# For the autocorrelation issue with segreg model
# get a and b from sr model
a <- as.numeric(sr.short@params["a"]) #  <- a <- 713.1137
b <- as.numeric(sr.short@params["b"]) # <- b <- mean(ssb(stk))

## S/R residuals - with autocorrelation
#ssb <- an(ssb(sr))
#pred.r <- ifelse(c(ssb) <= b, a * c(ssb), a * b)
#rec.res <- log(rec(sr)) - log(pred.r)
##plot(rec.res)
rec.short.res <- residuals(sr.short)

# autoregressive model order 1
set.seed(108)
# mean
arima.fit <- arima(an(rec.short.res), order = c(1, 0, 0))
# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.short.res <- make.arma.resid(arima.fit, age = 0, years = iy:(iy + ny-1), nit = it)
plot(sr.short.res)

# WITHOUT AUTOCORRELATION:
##sr.res1[] <- sample(c(residuals(sr)), ny*it, replace=TRUE)

# PLOTS

# Confidence interval plot
redstk <- trim(stk, year=1985:2015)

ssb <- apply(matrix(as.vector(redstk@stock.n*redstk@mat*redstk@stock.wt), nrow=1+redstk@range[[2]] ), 2, sum)

segreg.meanssb.red  <- function(ab, ssb) log(ifelse(ssb >= mean(ssb), ab$a*mean(ssb), ab$a*ssb))
fit.meanssb.red <- eqsr_fit(redstk,nsamp=2000, models = c("segreg.meanssb.red"))

png("Results/ANE/PLOTS/ANE_segreg_bkptmeanSSB_red_confidence.png", width=700, height=700)
eqsr_plot(fit.meanssb.red)
dev.off()

#############################################
# D. hockey stick with break point at maxSSB restricted time series
#############################################
# fit hockey stick
sr.max <- fmle(as.FLSR(trim(stk, year=1985:2015), model="segreg"), fixed=list(b=max(ssb(stk))))# method="L-BFGS-B"mean(ssb(stk))
#sr <- fmle(as.FLSR(stk, model="ricker")) # You could also try rickerAR1
#sr <- fmle(as.FLSR(stk, model="shepherd")) 

# calculate residuals of the fit and look
sr.max.res <- residuals(sr.max)
plot(sr.max.res)

# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
# For the autocorrelation issue with segreg model
# get a and b from sr model
a <- as.numeric(sr.max@params["a"]) #  <- a <- 713.1137
b <- as.numeric(sr.max@params["b"]) # <- b <- mean(ssb(stk))

## S/R residuals - with autocorrelation
#ssb <- an(ssb(sr))
#pred.r <- ifelse(c(ssb) <= b, a * c(ssb), a * b)
#rec.res <- log(rec(sr)) - log(pred.r)
##plot(rec.res)
rec.max.res <- residuals(sr.max)

# autoregressive model order 1
set.seed(108)
# mean
arima.fit <- arima(an(rec.max.res), order = c(1, 0, 0))
# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.max.res <- make.arma.resid(arima.fit, age = 0, years = iy:(iy + ny-1), nit = it)
plot(sr.max.res)

# WITHOUT AUTOCORRELATION:
##sr.res1[] <- sample(c(residuals(sr)), ny*it, replace=TRUE)

# PLOTS

# Confidence interval plot
redstk <- trim(stk, year=1985:2015)

ssb <- apply(matrix(as.vector(redstk@stock.n*redstk@mat*redstk@stock.wt), nrow=1+redstk@range[[2]] ), 2, sum)

segreg.maxssb.red  <- function(ab, ssb) log(ifelse(ssb >= max(ssb), ab$a*max(ssb), ab$a*ssb))
fit.maxssb.red <- eqsr_fit(redstk,nsamp=2000, models = c("segreg.maxssb.red"))

png("Results/ANE/PLOTS/ANE_segreg_bkptmaxSSB_red_confidence.png", width=700, height=700)
eqsr_plot(fit.maxssb.red)
dev.off()

# #############################################
# # E. hockey stick with break point at Bpa
# #############################################
# # fit hockey stick
# # sr.bkpt.bpa <- fmle(as.FLSR(trim(stk, year=1980:2015), model="segreg"), fixed=list(b=bpa))# method="L-BFGS-B"mean(ssb(stk))
# sr.bkpt.bpa <- fmle(as.FLSR(stk, model="segreg"), fixed=list(b=bpa))# method="L-BFGS-B"mean(ssb(stk))
# # calculate residuals of the fit and look
# sr.res.bkpt.bpa <- residuals(sr.bkpt.bpa)
# plot(sr.res.bkpt.bpa)
# 
# 
# ssb <- apply(matrix(as.vector(stk@stock.n*stk@mat*stk@stock.wt), nrow=1+stk@range[[2]] ), 2, sum)
# 
# segreg.bpa  <- function(ab, ssb) log(ifelse(ssb >= bpa, ab$a*bpa, ab$a*ssb))
# fit.bpa <- eqsr_fit(stk,nsamp=2000, models = c("segreg.bpa"))
# 
# png("Results/ANE/Pil_segreg_bkptBPA_confidence.png", width=700, height=700)
# eqsr_plot(fit.bpa)
# dev.off()
# 
# # I want to add something that takes autocorrelation in SR relationship into account
# # to do this I use an arima model
# # For the autocorrelation issue with segreg model
# # get a and b from sr model
# a.bkpt.bpa <- as.numeric(sr.bkpt.bpa@params["a"]) #  <- a <- 713.1137
# b.bkpt.bpa <- as.numeric(sr.bkpt.bpa@params["b"]) # <- b <- bpa
# 
# ## S/R residuals - with autocorrelation
# #ssb <- an(ssb(sr))
# #pred.r <- ifelse(c(ssb) <= b, a * c(ssb), a * b)
# #rec.res <- log(rec(sr)) - log(pred.r)
# ##plot(rec.res)
# rec.res.bkpt.bpa <- residuals(sr.bkpt.bpa)
# plot(rec.res.bkpt.bpa)
# 
# # autoregressive model order 1
# set.seed(108)
# # mean
# arima.fit.bkpt.bpa <- arima(an(rec.res.bkpt.bpa), order = c(1, 0, 0))
# # create autocorrelation in residuals and propagate throughout stock into the future
# # from initial year of projections (iy) to last of projections (ny-1)
# sr.res.bkpt.bpa <- make.arma.resid(arima.fit.bkpt.bpa, age = 0, years = iy:(iy + ny-1), nit = it)
# plot(sr.res.bkpt.bpa)
# 
# # WITHOUT AUTOCORRELATION:
# ##sr.res2[] <- sample(c(residuals(sr.bkpt.bpa.res)), ny*it, replace=TRUE)
# 
# #############################################
# # F. geometric mean of recruitment 
# #############################################
# 
# # geometric mean: CHECK fixed=list(a=1) with MIGUEL
# sr.geom.mean <- fmle(as.FLSR(stk, model="geomean"), fixed=list(a=1)) 
# plot(rec(ANCHOVY2idx)~ssb(ANCHOVY2idx))
# lines(ssb(ANCHOVY2idx), fitted(sr.geom.mean), col="red")
# 
# multi_rec_residuals_geomean <- FLQuant(NA, dimnames = list(year=iy:fy, iter=1:it))
# # We want to fill up our multi_rec_residuals FLQuant by randomly sampling from these log residuals
# sample_years <- sample(dimnames(residuals(sr.geomean.allSeries))$year, it * 16, replace = TRUE)
# multi_rec_residuals_geomean[] <- exp(residuals(sr.geomean.allSeries)[,sample_years])
# sr.geomean <- fmle(as.FLSR(window(stk, start=2000), model="geomean")) 
# 

}
 

#######################################################################################################################################

# Fixed objects
TAC <- FLQuant(NA, dimnames=list(TAC="all", year=vy, iter=1:it))
BB <- FLQuant(0, dimnames=list(TAC="all", year=vy, iter=1:it))

# # Prepare stock objects we need, with iterations and propagate towards final year
# DONE ABOVE            #stk <- iter(ANCHOVY2idx, 1)
                        # simulate new stock ased on a4a final fit
                        #sstk <- stk + simulate(fit, it)
                        #summary(sstk)
                        #plot(sstk)
                        
# short term forecast: start with a projection of F into the future to ny (16 yrs)
# this serves as a starting point for projecting the stock
pstk <- stf(sstk, ny, 3, 3)             # harvest is average last 3 years

landings.n(pstk) <- propagate(landings.n(pstk), it)           # NOT NEEDED already the right NB ITERS
discards.n(pstk) <- propagate(discards.n(pstk), it)

# # S/R residuals
# ##sr.res <- window(rec(pstk), iy, fy)
# ##sr.res[] <- sample(c(residuals(sr)), ny*it, replace=TRUE) 
# sr.res <- make.arma.resid(arima.fit, age = 0, years = iy:fy, nit = it)
# sr.res.bkpt.bpa <- make.arma.resid(arima.fit.bkpt.bpa, age = 0, years = iy:fy, nit = it)

# Prepare index object 
# Uncertainty around idx modelled as a lognormal distribution and propagated forward
# essentially it is calculating catchability with deviance around it

idx <- ids


# create an index of abundance using the error associated to the real index                   
skip.this <- T                                                           ### Should do that using the parameters estimated from SAM (catchability and observation variances), 
                                                                          ## corresponding to each iteration
if(skip.this)
{                                                                                             
          for (i in 1:length(idx)){
          cat(i)
            # align dimensions of the two FLQuant (index should represent stock in ideal world)
            # fill with NAs if not the same dims
            lst <- mcf(list(idx[[i]]@index, stock.n(stk)))
            # log of ratio between index and stock
            idx.lq <- log(lst[[1]]/lst[[2]])
            # if stock is 0 and index is a number it gives infinity so if the value is infinity put NA
            idx.lq[is.infinite(idx.lq)] <- NA # fix zeros
            # select the first iteration of stock.n and assign it to idx.qmu and idx.qsig to give the same dimensions to the FLQuant
            idx.qmu <- idx.qsig <- stock.n(iter(pstk,1))
            # year means of the log ratios
            idx.qmu[] <- yearMeans(idx.lq)
            # cv of log ratio : variance /means
            idx.qsig[] <- log((sqrt(yearVars(idx.lq))/yearMeans(idx.lq))^2 + 1)
            # make two empty FLQuants of stock.n dims
            idx.q <- idx_temp <- FLQuant(NA, dimnames=dimnames(stock.n(pstk)))
            # propagate exp(log ratios) into future for all iterations
            idx.q[,ac(y0:dy)] <- propagate(exp(idx.lq[,ac(y0:dy)]), it)
            # where it is not NA, generate random deviates from the log normal distribution of idx.qmu and and idx.qsig
            idx.q[!is.na(idx.qmu)] <- rlnorm(it, idx.qmu[!is.na(idx.qmu)], idx.qsig[!is.na(idx.qmu)])
            plot(idx.q)
            # idx is, in effect, your catchability, so you multipy it by the stock to get the index
            idx_temp <- idx.q * stock.n(pstk)
            # creqate new FLIndex object
            idx[[i]] <- FLIndex(index=idx_temp, index.q=idx.q)
            # set the timing of the survey to 0 - WHY?????
            #  range(idx[[i]])[c("startf", "endf")] <- c(0, 0)
            plot(index(idx[[i]]))
          }
}

## DO it using SAM formalism

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
      
      fleet.type <- sam.ctrl@fleets[n.]
      # fill in the catchability slot for age structured surveys
      if (fleet.type == 2 )
      {
      agesQ     <-  names(sam.ctrl@catchabilities[n.,]) [!is.na(sam.ctrl@catchabilities[n.,])]
      agesSdObs <-  names(sam.ctrl@obs.vars[n.,]) [!is.na(sam.ctrl@obs.vars[n.,])]
      
      for (its in 1:it) 
          {
          iter( index.q(x)[agesQ,]   , its )     <- FLQuant( matrix(rep(exp(logQ[its,sam.ctrl@catchabilities[n.,agesQ]]),dim(index.q(x)[agesQ,])[2]),nrow = dim(index.q(x)[agesQ,])[1]), dimnames = dimnames(iter( index.q(x)[agesQ,] , its )))
          iter( index.var(x)[agesSdObs,] , its ) <- FLQuant( matrix(rep(exp(logSdObs[its,sam.ctrl@obs.vars[n.,agesSdObs]]),dim(index.var(x)[agesSdObs,])[2]),nrow = dim(index.var(x)[agesSdObs,])[1]), dimnames = dimnames(iter( index.var(x)[agesSdObs,] , its )))
          }
      }
      # fill in the catchability slot for SSB surveys
      if (fleet.type == 3 )
      {
      for (its in 1:it) 
          {
          iter( index.q(x)  , its )     <- FLQuant( matrix(rep(exp(logQSsb[its,sam.ctrl@catchabilities[n.,agesQ]]),dim(index.q(x)[agesQ,])[2]),nrow = dim(index.q(x)[agesQ,])[1]), dimnames = dimnames(iter( index.q(x)[agesQ,] , its )))
          iter( index.var(x)[agesSdObs,] , its ) <- FLQuant( matrix(rep(exp(logSdObs[its,sam.ctrl@obs.vars[n.,agesSdObs]]),dim(index.var(x)[agesSdObs,])[2]),nrow = dim(index.var(x)[agesSdObs,])[1]), dimnames = dimnames(iter( index.var(x)[agesSdObs,] , its )))
          }
      }      
      
      })




#------------------------------------------------------------------------------
# 1a. (1a) STATUS QUO SCENARIO - SEGMENTED STOCK RECRUITMENT WITH BPT AT MEAN SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- 45936
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
  # pstk is at the beginning only populated into the future for F
  # the rest is only 1975-2016 but as the loop progresses through the projection
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
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
ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0 <- pstk
ANE_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0 <- stk0
plot(ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
save(list = c("ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0", "ANE_Opt1a_mse.stk0.GFCM_segregmeanSSB_Fsq0"), file = "Results/ANE/ANE_Opt1a_mse_SegregmeanSSB_statusQuo_250it.RData")
png("ANE_Opt1a_mse_SegregmeanSSB_statusQuo_250it.png", width=700, height=700)
plot(ANE_Opt1a_mse.pstk.GFCM_segregmeanSSB_Fsq0)
dev.off()














# #------------------------------------------------------------------------------
# # 1b. STATUS QUO - SEGMENTED STOCK RECRUITMENT WITH BPT AT BPA
# #------------------------------------------------------------------------------
# 
# # Set up the Btrigger (in this case halfway between Blim and Bpa)
# blim <- min(ssb(stk))
# bpa <- blim*2
# Btrig <- blim+((bpa-blim)/2)
# idx0 <- idx
# dt <- date()
# 
# #########################################################
# # go fish
# #Btrig <- blim+((bpa-blim)/2)
# idx0 <- idx
# dt <- date()
# 
# # go fish
# for(i in vy[-length(vy)]){
#   ## i <- vy[-length(vy)][1]
#   print(i)
#   gc()
#   ay <- an(i)   # an is equivalent to as.numeric
#   cat(i, "\n")
#   vy0 <- 1:(ay-y0+1) 
#   sqy <- (ay-y0-3+1+1):(ay-y0+1)
#   stk0 <- pstk[,vy0]
#   catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
#   ## note that vy0 is changing below so index is being updated
#   for (index_counter in 1:length(idx)){
#     idx0[[index_counter]] <- idx[[index_counter]][,vy0]
#     index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
#   }
#   ##
#   qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
#   fmod4 <- ~ s(age, k = 5) + factor(year)
#   rmodel4 <- ~s(year, k=20)
#   fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
#   stk0 <- stk0 + fit
#   # fwd control
#   fsq0 <- yearMeans(fbar(stk0)[,38:40]) # status quo 2012-2014
#   dnms <- list(iter=1:it, year=c(ay, ay + 1, ay + 1), c("min", "val", "max"))
#   arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
#   ## ftrg.vec <- rep(ftrg, it) ## original
#   refpt <- data.frame(harvest = 1)
#   ftrg.vec <- an(fsq0) # Ftarget = status quo
#   #Bescape <- blim
#   arr0[,,"val"] <- c(fsq0, ftrg.vec, rep(NA, it))
#   #arr0[,,"min"] <- c(rep(NA, 2 * it), rep(Bescape, it))
#   #arr0 <- aperm(arr0, c(2,3,1))
#   ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
#   ctrl@trgtArray <- arr0
#   ## 
#   stkTmp <- stf(stk0, 2)
#   stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.bkpt.bpa, sr.residuals = exp(sr.res.bkpt.bpa[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
#   TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
#   # OM proj
#   ctrl@target <- ctrl@target[2,]
#   ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
#   pstk <- fwd(pstk, ctrl=ctrl, sr=sr.bkpt.bpa, sr.residuals = exp(sr.res.bkpt.bpa[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
# }
# 
# return(val)
# date()
# ANE_Opt1b_mse.pstk.GFCM_segreg_Fsq0_bkptBpa <- pstk
# ANE_Opt1b_mse.stk0.GFCM_segreg_Fsq0_bkptBpa <- stk0
# plot(ANE_Opt1b_mse.pstk.GFCM_segreg_Fsq0_bkptBpa)
# save(list = c("ANE_Opt1b_mse.pstk.GFCM_segreg_Fsq0_bkptBpa", "ANE_Opt1b_mse.stk0.GFCM_segreg_Fsq0_bkptBpa"), file = "ANE_Opt1b_mse_Segreg_statusQuo_100it_bkptBpa.RData")
# png("ANE_Opt1_mse_Segreg_statusQuo_100it_bkptBpa.png", width=700, height=700)
# plot(ANE_Opt1b_mse.pstk.GFCM_segreg_Fsq0_bkptBpa)
# dev.off()



# #------------------------------------------------------------------------------
# # 1c. STATUS QUO WITH GEOM MEAN
# #------------------------------------------------------------------------------

# # Set up the Btrigger (in this case halfway between Blim and Bpa)
# blim <- min(ssb(stk))
# bpa <- blim*2
# Btrig <- blim+((bpa-blim)/2)
# idx0 <- idx
# dt <- date()

# #########################################################
# # go fish
# #Btrig <- blim+((bpa-blim)/2)
# idx0 <- idx
# dt <- date()

# # go fish
# for(i in vy[-length(vy)]){
#   ## i <- vy[-length(vy)][1]
#   print(i)
#   gc()
#   ay <- an(i)   # an is equivalent to as.numeric
#   cat(i, "\n")
#   vy0 <- 1:(ay-y0) # data years (positions vector)
#   sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
#   stk0 <- pstk[,vy0]
#   catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
#   ## note that vy0 is changing below so index is being updated
#   for (index_counter in 1:length(idx)){
#     idx0[[index_counter]] <- idx[[index_counter]][,vy0]
#     index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
#   }
#   ##
#   qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
#   fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
#   rmodel4 <- ~ s(year, k=20)  
#   fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5)
#   stk0 <- stk0 + fit
#   # fwd control
#   fsq0 <- yearMeans(fbar(stk)[,38:40]) # status quo 2012-2014
#   dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
#   arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
#   ## ftrg.vec <- rep(ftrg, it) ## original
#   refpt <- data.frame(harvest = 1)
#   ftrg.vec <- rep(an(fsq0),it) # Ftarget = status quo
#   #Bescape <- blim
#   arr0[,,"val"] <- c(rep(an(fsq0),it), ftrg.vec)
#   #arr0[,,"min"] <- c(rep(NA, 2 * it), rep(Bescape, it))
#   #arr0 <- aperm(arr0, c(2,3,1))
#   ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
#   ctrl@trgtArray <- arr0
#   ## 
#   stkTmp <- stf(stk0, 2)
#   stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.geomean, sr.residuals = multi_rec_residuals_geomean[,ac(ay:(ay+1))], sr.residuals.mult = TRUE) #
#   TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
#   # OM proj
#   ctrl@target <- ctrl@target[2,]
#   ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
#   pstk <- fwd(pstk, ctrl=ctrl, sr=sr.geomean, sr.residuals = multi_rec_residuals_geomean[,ac(ay:(ay+1))], sr.residuals.mult = TRUE) #
# }

# return(val)
# date()
# Opt1b_mse.pstk.GFCM_geomean_Fsq0 <- pstk
# Opt1b_mse.stk0.GFCM_geomean_Fsq0 <- stk0
# plot(Opt1b_mse.pstk.GFCM_geomean_Fsq0)
# save(list = c("Opt1b_mse.pstk.GFCM_geomean_Fsq0", "Opt1b_mse.stk0.GFCM_geomean_Fsq0"), file = "PIL_Opt1b_mse_geomean_statusQuo_250it_UncertAllTimeSeries.RData")
# png("PIL_Opt1b_mse_Geomean_statusQuo_250it_UncertAllTimeSeries.png", width=700, height=700)
# plot(Opt1b_mse.pstk.GFCM_geomean_Fsq0)
# dev.off()

#------------------------------------------------------------------------------
# 1d. (1c) STATUS QUO SCENARIO - SEGMENTED RESTRITED STOCK RECRUITMENT WITH BPT AT MEAN SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- 45936
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
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
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.short, sr.residuals = exp(sr.short.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  # update pstk with stkTmp
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.short, sr.residuals = exp(sr.short.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
}


return(val)
date()

# breakpoint mean(SSB)
ANE_Opt1d_mse.pstk.GFCM_segregshort_meanSSB_Fsq0 <- pstk
ANE_Opt1d_mse.stk0.GFCM_segregshort_meanSSB_Fsq0 <- stk0
#plot(ANE_Opt1d_mse.pstk.GFCM_segregshort_meanSSB_Fsq0)
save(list = c("ANE_Opt1d_mse.pstk.GFCM_segregshort_meanSSB_Fsq0", "ANE_Opt1d_mse.stk0.GFCM_segregshort_meanSSB_Fsq0"), file = "Results/ANE/ANE_Opt1d_mse_SegregShort_meanSSB_statusQuo_250it.RData")
# png("ANE_Opt1d_mse_SegregShort_meanSSB_statusQuo_250it.png", width=700, height=700)
# plot(ANE_Opt1d_mse.pstk.GFCM_segregshort_meanSSB_Fsq0)
# dev.off()


#------------------------------------------------------------------------------
# 1e. (1d) STATUS QUO SCENARIO - SEGMENTED RESTRITED STOCK RECRUITMENT WITH BPT AT MAX SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- 45936
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
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
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.max, sr.residuals = exp(sr.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  # update pstk with stkTmp
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.max, sr.residuals = exp(sr.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
}


return(val)
date()

# breakpoint mean(SSB)
ANE_Opt1e_mse.pstk.GFCM_segregshort_maxSSB_Fsq0 <- pstk
ANE_Opt1e_mse.stk0.GFCM_segregshort_maxSSB_Fsq0 <- stk0
#plot(ANE_Opt1e_mse.pstk.GFCM_segregshort_maxSSB_Fsq0)
save(list = c("ANE_Opt1e_mse.pstk.GFCM_segregshort_maxSSB_Fsq0", "ANE_Opt1e_mse.stk0.GFCM_segregshort_maxSSB_Fsq0"), file = "Results/ANE/ANE_Opt1e_mse_SegregShort_maxSSB_statusQuo_250it.RData")
# png("ANE_Opt1e_mse_SegregShort_maxSSB_statusQuo_250it.png", width=700, height=700)
# plot(ANE_Opt1e_mse.pstk.GFCM_segregshort_maxSSB_Fsq0)
# dev.off()


#---------------------------------------------------------------------------------------
# 1f. (1b) STATUS QUO SCENARIO - SEGMENTED STOCK RECRUITMENT WITH BPT AT MAX SSB entire TS
#---------------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- 45936
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
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
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.all.max, sr.residuals = exp(sr.all.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  # update pstk with stkTmp
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.all.max, sr.residuals = exp(sr.all.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
}


return(val)
date()

# breakpoint mean(SSB)
ANE_Opt1f_mse.pstk.GFCM_segregmaxSSB_Fsq0 <- pstk
ANE_Optf_mse.stk0.GFCM_segregmaxSSB_Fsq0 <- stk0
#plot(ANE_Opt1f_mse.pstk.GFCM_segregmaxSSB_Fsq0)
save(list = c("ANE_Opt1f_mse.pstk.GFCM_segregmaxSSB_Fsq0", "ANE_Opt1f_mse.stk0.GFCM_segregmaxSSB_Fsq0"), file = "Results/ANE/ANE_Opt1f_mse_SegregmaxSSB_statusQuo_250it.RData")
# png("ANE_Opt1f_mse_SegregmaxSSB_statusQuo_250it.png", width=700, height=700)
# plot(ANE_Opt1f_mse.pstk.GFCM_segregmaxSSB_Fsq0)
# dev.off()


#------------------------------------------------------------------------------
# 2a. (5a) SCENARIO GFCM REGULATION 2013 Halfway Blim-Bpa = Btrigger - USING F INCREASING LINEARLY WHEN B > 
# AND WHEN F < 0.64 (with and without Bescape)
# entire TS bkpt mean SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.64
trgyF <- 2020 # you want your target F to be reached by 2020, no need for it to be reached straight away
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, trgyF=trgyF)
  #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
  ftrg.vec <- an(ftrg.q)
  # Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
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
return(val)
date()

# # with bescape
# Opt2_ANE_mse.pstk.GFCM_segreg_Btrig <- pstk
# Opt2_ANE_mse.stk0.GFCM_segreg_Btrig <- stk0
# plot(Opt2_ANE_mse.pstk.GFCM_segreg_Btrig)
# save(list = c("Opt2_ANE_mse.pstk.GFCM_segreg_Btrig", "Opt2_ANE_mse.stk0.GFCM_segreg_Btrig"), file = "Results/ANE/ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.RData")
# png("ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.png", width=700, height=700)
# plot(ANE_Opt2_mse.pstk.GFCM_segreg_Btrig)
# dev.off()

# without bescape
ANE_Opt2a_mse.pstk.GFCM_segreg_Btrig_noBesc <- pstk
ANE_Opt2a_mse.stk0.GFCM_segreg_Btrig_noBesc <- stk0
save(list = c("ANE_Opt2a_mse.pstk.GFCM_segreg_Btrig_noBesc", "ANE_Opt2a_mse.stk0.GFCM_segreg_Btrig_noBesc"), file = "Results/ANE/ANE_Opt2a_mse_Segreg_GFCM-HCR_Btrig_250it_noBesc.RData")
# png("ANE_Opt2a_mse_Segreg_GFCM-HCR_Btrig_250it_noBesc.png", width=700, height=700)
# plot(ANE_Opt2a_mse.pstk.GFCM_segreg_Btrig_noBesc)
# dev.off()

#------------------------------------------------------------------------------
# 2b. (5b) SCENARIO GFCM REGULATION 2013 Halfway Blim-Bpa = Btrigger - USING F INCREASING LINEARLY WHEN B > 
# AND WHEN F < 0.64 (with and without Bescape)
# entire TS bkpt max SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.64
trgyF <- 2020 # you want your target F to be reached by 2020, no need for it to be reached straight away
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, trgyF=trgyF)
  #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
  ftrg.vec <- an(ftrg.q)
  # Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
  ## 
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.all.max, sr.residuals = exp(sr.all.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.all.max, sr.residuals = exp(sr.all.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}

##
return(val)
date()

# # with bescape
# Opt2_ANE_mse.pstk.GFCM_segreg_Btrig <- pstk
# Opt2_ANE_mse.stk0.GFCM_segreg_Btrig <- stk0
# plot(Opt2_ANE_mse.pstk.GFCM_segreg_Btrig)
# save(list = c("Opt2_ANE_mse.pstk.GFCM_segreg_Btrig", "Opt2_ANE_mse.stk0.GFCM_segreg_Btrig"), file = "Results/ANE/ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.RData")
# png("ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.png", width=700, height=700)
# plot(ANE_Opt2_mse.pstk.GFCM_segreg_Btrig)
# dev.off()

# without bescape
ANE_Opt2b_mse.pstk.GFCM_segreg_maxSSB_Btrig_noBesc <- pstk
ANE_Opt2b_mse.stk0.GFCM_segreg_maxSSB_Btrig_noBesc <- stk0
save(list = c("ANE_Opt2b_mse.pstk.GFCM_segreg_maxSSB_Btrig_noBesc", "ANE_Opt2b_mse.stk0.GFCM_segreg_maxSSB_Btrig_noBesc"), file = "ANE_Opt2b_mse_Segreg_maxSSB_GFCM-HCR_Btrig_250it_noBesc.RData")
# png("ANE_Opt2b_mse_Segreg_maxSSB_GFCM-HCR_Btrig_250it_noBesc.png", width=700, height=700)
# plot(ANE_Opt2b_mse.pstk.GFCM_segreg_maxSSB_Btrig_noBesc)
# dev.off()


#------------------------------------------------------------------------------
# 2c. (5c) SCENARIO GFCM REGULATION 2013 Halfway Blim-Bpa = Btrigger - USING F INCREASING LINEARLY WHEN B > 
# AND WHEN F < 0.64 (with and without Bescape)
# reduced TS bkpt mean SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.64
trgyF <- 2020 # you want your target F to be reached by 2020, no need for it to be reached straight away
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, trgyF=trgyF)
  #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
  ftrg.vec <- an(ftrg.q)
  # Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
  ## 
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.short, sr.residuals = exp(sr.short.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.short, sr.residuals = exp(sr.short.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}



##
return(val)
date()

# # with bescape
# Opt2_ANE_mse.pstk.GFCM_segreg_Btrig <- pstk
# Opt2_ANE_mse.stk0.GFCM_segreg_Btrig <- stk0
# plot(Opt2_ANE_mse.pstk.GFCM_segreg_Btrig)
# save(list = c("Opt2_ANE_mse.pstk.GFCM_segreg_Btrig", "Opt2_ANE_mse.stk0.GFCM_segreg_Btrig"), file = "Results/ANE/ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.RData")
# png("ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.png", width=700, height=700)
# plot(ANE_Opt2_mse.pstk.GFCM_segreg_Btrig)
# dev.off()

# without bescape
ANE_Opt2c_mse.pstk.GFCM_segreg_red_meanSSB_Btrig_noBesc <- pstk
ANE_Opt2c_mse.stk0.GFCM_segreg_red_meanSSB_Btrig_noBesc <- stk0
save(list = c("ANE_Opt2c_mse.pstk.GFCM_segreg_red_meanSSB_Btrig_noBesc", "ANE_Opt2c_mse.stk0.GFCM_segreg_red_meanSSB_Btrig_noBesc"), file = "ANE_Opt2c_mse_Segreg_red_meanSSB_GFCM-HCR_Btrig_250it_noBesc.RData")
#png("ANE_Opt2c_mse_Segreg_red_meanSSB_GFCM-HCR_Btrig_250it_noBesc.png", width=700, height=700)
# plot(ANE_Opt2c_mse.pstk.GFCM_segreg_red_meanSSB_Btrig_noBesc)
# dev.off()



#------------------------------------------------------------------------------
# 2d. (5d) SCENARIO GFCM REGULATION 2013 Halfway Blim-Bpa = Btrigger - USING F INCREASING LINEARLY WHEN B > 
# AND WHEN F < 0.64 (with and without Bescape)
# reduced TS bkpt max SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.64
trgyF <- 2020 # you want your target F to be reached by 2020, no need for it to be reached straight away
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, trgyF=trgyF)
  #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
  ftrg.vec <- an(ftrg.q)
  # Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
  ## 
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.max, sr.residuals = exp(sr.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.max, sr.residuals = exp(sr.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}

##
return(val)
date()

# # with bescape
# Opt2_ANE_mse.pstk.GFCM_segreg_Btrig <- pstk
# Opt2_ANE_mse.stk0.GFCM_segreg_Btrig <- stk0
# plot(Opt2_ANE_mse.pstk.GFCM_segreg_Btrig)
# save(list = c("Opt2_ANE_mse.pstk.GFCM_segreg_Btrig", "Opt2_ANE_mse.stk0.GFCM_segreg_Btrig"), file = "Results/ANE/ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.RData")
# png("ANE_Opt2_mse_Segreg_GFCM-HCR_Btrig_250it.png", width=700, height=700)
# plot(ANE_Opt2_mse.pstk.GFCM_segreg_Btrig)
# dev.off()

# without bescape
ANE_Opt2d_mse.pstk.GFCM_segreg_red_maxSSB_Btrig_noBesc <- pstk
ANE_Opt2d_mse.stk0.GFCM_segreg_red_maxSSB_Btrig_noBesc <- stk0
save(list = c("ANE_Opt2d_mse.pstk.GFCM_segreg_red_maxSSB_Btrig_noBesc", "ANE_Opt2d_mse.stk0.GFCM_segreg_red_maxSSB_Btrig_noBesc"), file = "Results/ANE/ANE_Opt2d_mse_Segreg_red_maxSSB_GFCM-HCR_Btrig_250it_noBesc.RData")
# png("ANE_Opt2d_mse_Segreg_red_maxSSB_GFCM-HCR_Btrig_250it_noBesc.png", width=700, height=700)
# plot(ANE_Opt2d_mse.pstk.GFCM_segreg_red_maxSSB_Btrig_noBesc)
# dev.off()



#------------------------------------------------------------------------------
# 3a. SCENARIO GFCM REGULATION 2013 - ALLOWING F TO INCREASE WHEN B > BPA BUT NOT WHEN
# F BIGGER THAN 0.64 (CORRESPONDS TO E=0.4)
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
Ftar <- 0.64
dt <- date()
trgyF <- 2020


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
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, trgyF=trgyF)
  #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
  ftrg.vec <- an(ftrg.q)
  Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
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
return(val)
date()
ANE_Opt3a_mse.pstk.GFCM_segreg_BtrigBpa <- pstk
ANE_Opt3a_mse.stk0.GFCM_segreg_BtrigBpa <- stk0
save(list = c("Pil_Opt3a_mse.pstk.GFCM_segreg_BtrigBpa", "Pil_Opt3a_mse.stk0.GFCM_segreg_BtrigBpa"), file = "Results/ANE/Pil_Opt3_mse_Segreg_GFCM-HCR_Fincrease_BtrigBpa.RData")
png("Pil_Opt3_mse_Segreg_GFCM-HCR_Fincrease_BtrigBpa.png", width=700, height=700)
plot(ANE_Opt3a_mse.pstk.GFCM_segreg_BtrigBpa)
dev.off()


stock(pstk)

# check out single iterations just to see if the rule is being applied most of the time
# it won't work every time because it takes into account the uncertainty in the data

# take random iterations and save them in p
p <- pstk[,,,,,(ceiling(runif(5, 1, 50)))]
stock <- as.data.frame(stock(p[,,,,,1]))[,c(2,7)]
stock$i2 <-as.data.frame(stock(p[,,,,,2]))[,7]
stock$i3 <-as.data.frame(stock(p[,,,,,3]))[,7]
stock$i4 <-as.data.frame(stock(p[,,,,,4]))[,7]
stock$i5 <-as.data.frame(stock(p[,,,,,5]))[,7]

ssb <- as.data.frame(ssb(p[,,,,,1]))[,c(2,7)]
ssb$i2 <-as.data.frame(ssb(p[,,,,,2]))[,7]
ssb$i3 <-as.data.frame(ssb(p[,,,,,3]))[,7]
ssb$i4 <-as.data.frame(ssb(p[,,,,,4]))[,7]
ssb$i5 <-as.data.frame(ssb(p[,,,,,5]))[,7]

fbar <- as.data.frame(fbar(p[,,,,,1]))[,c(2,7)]
fbar$i2 <-as.data.frame(fbar(p[,,,,,2]))[,7]
fbar$i3 <-as.data.frame(fbar(p[,,,,,3]))[,7]
fbar$i4 <-as.data.frame(fbar(p[,,,,,4]))[,7]
fbar$i5 <-as.data.frame(fbar(p[,,,,,5]))[,7]

landings <- as.data.frame(landings(p[,,,,,1]))[,c(2,7)]
landings$i2 <-as.data.frame(landings(p[,,,,,2]))[,7]
landings$i3 <-as.data.frame(landings(p[,,,,,3]))[,7]
landings$i4 <-as.data.frame(landings(p[,,,,,4]))[,7]
landings$i5 <-as.data.frame(landings(p[,,,,,5]))[,7]

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 4a. (7a) SCENARIO BASED BOTH ON F AND BIO - CHOOSE THE MORE PRECAUTIONARY OPTION
# entire TS bkpt mean SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
Fmsy <- 0.554
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
ANE_Opt4a_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it <- pstk
ANE_Opt4a_mse.stk0.GFCM_segreg_BtarFtar_trgy_250it <- stk0
save(list = c("ANE_Opt4a_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it", "ANE_Opt4a_mse.stk0.GFCM_segreg_BtarFtar_trgy_250it"), file = "Results/ANE/ANE_Opt4a_mse_Segreg_FtarBtar_250it.RData")
# png("ANE_Opt4a_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it.png", width=700, height=700)
# plot(window(ANE_Opt4a_mse.pstk.GFCM_segreg_BtarFtar_trgy_250it, end=2030)) + ggtitle("Anchovy - Target Fmsy & recover to Bpa") + theme_light(base_size=12)
# dev.off()

#-----------------------------------------------------------------------------------------------------
# 4b. (7b) SCENARIO BASED BOTH ON F AND BIO - CHOOSE THE MORE PRECAUTIONARY OPTION 
# SEGREG max SSB
#-----------------------------------------------------------------------------------------------------

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
  stkTmp10 <- fwd(stkTmp[,,,,,], ctrl=ctrl, sr=sr.all.max, sr.residuals = exp(sr.all.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #,exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2:3,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2:3,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.all.max, sr.residuals = exp(sr.all.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}


##
return(val)
date()
ANE_Opt4b_mse.pstk.GFCM_segreg_maxSSB_BtarFtar_trgy_75it <- pstk
ANE_Opt4b_mse.stk0.GFCM_segreg_maxSSB_BtarFtar_trgy_75it <- stk0
save(list = c("ANE_Opt4b_mse.pstk.GFCM_segreg_maxSSB_BtarFtar_trgy_75it", "ANE_Opt4b_mse.stk0.GFCM_segreg_maxSSB_BtarFtar_trgy_75it"), file = "Results/ANE/ANE_Opt4b_mse_Segreg_maxSSB_FtarBtar_75it.RData")
# png("ANE_Opt4a_mse.pstk.GFCM_segreg_maxSSB_BtarFtar_trgy_75it.png", width=700, height=700)
# plot(window(ANE_Opt4b_mse.pstk.GFCM_segreg_maxSSB_BtarFtar_trgy_75it, end=2030)) + ggtitle("Anchovy - Target Fmsy & recover to Bpa") + theme_light(base_size=12)
# dev.off()

#------------------------------------------------------------------------------
# 4c. (7c) SCENARIO BASED BOTH ON F AND BIO - CHOOSE THE MORE PRECAUTIONARY OPTION
# reduced TS bkpt mean SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)

Btrig <- bpa
idx0 <- idx
Fmsy <- 0.554
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
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
  stkTmp10 <- fwd(stkTmp[,,,,,], ctrl=ctrl, sr=sr.short, sr.residuals = exp(sr.short.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #,exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2:3,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2:3,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.short, sr.residuals = exp(sr.short.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}

##
return(val)
date()
ANE_Opt4c_mse.pstk.GFCM_segreg_redmeanSSB_BtarFtar_trgy_250it <- pstk
ANE_Opt4c_mse.stk0.GFCM_segreg_redmeanSSB_BtarFtar_trgy_250it <- stk0
save(list = c("ANE_Opt4c_mse.pstk.GFCM_segreg_redmeanSSB_BtarFtar_trgy_250it", "ANE_Opt4c_mse.stk0.GFCM_segreg_redmeanSSB_BtarFtar_trgy_250it"), file = "Results/ANE/ANE_Opt4a_mse_Segreg_redmeanSSB_FtarBtar_250it.RData")
# png("ANE_Opt4c_mse.pstk.GFCM_segreg_redmeanSSB_BtarFtar_trgy_250it.png", width=700, height=700)
# plot(window(ANE_Opt4c_mse.pstk.GFCM_segreg_redmeanSSB_BtarFtar_trgy_250it, end=2030)) + ggtitle("Anchovy - Target Fmsy & recover to Bpa") + theme_light(base_size=12)
# dev.off()


#------------------------------------------------------------------------------
# 4d. (7d) SCENARIO BASED BOTH ON F AND BIO - CHOOSE THE MORE PRECAUTIONARY OPTION
# reduced TS bkpt max SSB
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)

Btrig <- bpa
idx0 <- idx
Fmsy <- 0.554
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
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
  stkTmp10 <- fwd(stkTmp[,,,,,], ctrl=ctrl, sr=sr.max, sr.residuals = exp(sr.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #,exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
  # OM proj
  ctrl@target <- ctrl@target[2:3,]
  ## original was catch
  ##ctrl@target[,"quantity"] <- "catch"
  ctrl@trgtArray <- ctrl@trgtArray[2:3,,,drop=FALSE]
  ## original was catch
  ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr.max, sr.residuals = exp(sr.max.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}

##
return(val)
date()
ANE_Opt4d_mse.pstk.GFCM_segreg_redmaxSSB_BtarFtar_trgy_250it <- pstk
ANE_Opt4d_mse.stk0.GFCM_segreg_redmaxSSB_BtarFtar_trgy_250it <- stk0
save(list = c("ANE_Opt4d_mse.pstk.GFCM_segreg_redmaxSSB_BtarFtar_trgy_250it", "ANE_Opt4d_mse.stk0.GFCM_segreg_redmaxSSB_BtarFtar_trgy_250it"), file = "Results/ANE/ANE_Opt4d_mse_Segreg_redmaxSSB_FtarBtar_250it.RData")
# png("ANE_Opt4d_mse.pstk.GFCM_segreg_redmaxSSB_BtarFtar_trgy_250it.png", width=700, height=700)
# plot(window(ANE_Opt4d_mse.pstk.GFCM_segreg_redmaxSSB_BtarFtar_trgy_250it, end=2030)) + ggtitle("Anchovy - Target Fmsy & recover to Bpa") + theme_light(base_size=12)
# dev.off()


#------------------------------------------------------------------------------
# 5. (8) STATUS QUO TO FMSY
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
dt <- date()
Fmsy <- 0.554
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  ftrg.q <- Fmsy 
  ftrg.vec <- rep(Fmsy,it)
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- c(an(fsq0), ftrg.vec) #rep(NA, it)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
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
ANE_Opt5_mse.pstk.GFCM_segreg_Fmsy_SQ_250it <- pstk
ANE_Opt5_mse.stk0.GFCM_segreg_Fmsy_SQ_250it <- stk0
plot(ANE_Opt5_mse.pstk.GFCM_segreg_Fmsy_SQ_250it)
save(list = c("ANE_Opt5_mse.pstk.GFCM_segreg_Fmsy_SQ_250it", "ANE_Opt5_mse.stk0.GFCM_segreg_Fmsy_SQ_250it"), file = "Results/ANE/ANE_Opt5_mse_segreg_Fmsy_SQ_250it.RData")
png("ANE_Opt5_mse_Segreg_Fmsy_SQ_250it.png", width=700, height=700)
plot(ANE_Opt5_mse.pstk.GFCM_segreg_Fmsy_SQ_250it)
dev.off()



# #------------------------------------------------------------------------------
# # 6a. ENFORCEMENT OF FISHING DAYS REDUCTION - REDUCTION OF 9% OF F IN 2015
# #------------------------------------------------------------------------------
# 
# # Set up the Btrigger (in this case Bpa)
# Btrig <- bpa
# idx0 <- idx
# dt <- date()
# 
# #########################################################
# # go fish
# for(i in vy[-length(vy)]){
#   ## i <- vy[-length(vy)][1]
#   print(i)
#   gc()
#   ay <- an(i)   # an is equivalent to as.numeric
#   cat(i, "\n")
#   vy0 <- 1:(ay-y0) # data years (positions vector)
#   sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
#   stk0 <- pstk[,vy0]
#   catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
#   ## note that vy0 is changing below so index is being updated
#   for (index_counter in 1:length(idx)){
#     idx0[[index_counter]] <- idx[[index_counter]][,vy0]
#     index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
#   }
#   ##
#   qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
#   fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
#   rmodel4 <- ~ s(year, k=20)  
#   fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5) #, qmodel=qmod5)
#   stk0 <- stk0 + fit
#   # fwd control
#   fsq0 <- yearMeans(fbar(stk0)[,c(sqy)])
#   dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
#   arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
#   ## ftrg.vec <- rep(ftrg, it) ## original
#   #refpt <- data.frame(ssb = 1, harvest = 1)
#   #ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, refpt = refpt, Btrig = Btrig, Fmin = 0, Blim = blim, Bpa=bpa)
#   #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
#   #ftrg.vec <- an(ftrg.q)
#   #Bescape <- blim
#   ftrg.q <- fbar(stk)[,"2014",,,,] * 0.91 
#   ftrg.vec <- rep(an(ftrg.q),it)
#   arr0[,,"val"] <- c(fsq0, ftrg.vec) #rep(NA, it)
#   arr0[,,"min"] <- c(rep(NA, 2 * it))
#   arr0 <- aperm(arr0, c(2,3,1))
#   ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
#   ctrl@trgtArray <- arr0
#   #future_catch <- c(catch(stk0)[,"2013"]) * 0.9
#   #ctrl_catch <- fwdControl(data.frame(year=an(ay:(ay+1)), quantity = "catch", val=future_catch))
#   #ctrl_target <- ctrl_target[order(ctrl_target$year),]
#   #ctrl <- fwdControl(ctrl_catch)
#   #ctrl <- fwdControl(data.frame(year=c(ay, ay+1, ay + 3), quantity=c('f', 'f', 'ssb'), val=NA))
#   #ctrl@trgtArray <- arr0
#   ## 
#   stkTmp <- stf(stk0, 2)
#   stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
#   ## USING F
#   TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
#   # OM proj
#   ctrl@target <- ctrl@target[2,]
#   ## original was catch
#   ##ctrl@target[,"quantity"] <- "catch"
#   ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
#   ## original was catch
#   ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
#   pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
#   #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
# }
# 
# 
# 
# ##
# 
# date()
# PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it <- pstk
# PIL_Opt6a_mse.stk0.GFCM_segreg_9perc_red_250it <- stk0
# plot(PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it)
# save(list = c("PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it", "PIL_Opt6a_mse.stk0.GFCM_segreg_9perc_red_250it"), file = "Results/ANE/PIL_Opt6a_mse_segreg_9perc_red_250it.RData")
# png("PIL_Opt6a_mse_Segreg__9perc_red_250it_250it.png", width=700, height=700)
# plot(PIL_Opt6a_mse.pstk.GFCM_segreg_9perc_red_250it)
# dev.off()
# 
# 
# 
# #------------------------------------------------------------------------------
# # 6b. ENFORCEMENT OF FISHING DAYS REDUCTION - REDUCTION OF 9% OF CATCH IN 2015
# #------------------------------------------------------------------------------
# 
# # Set up the Btrigger (in this case Bpa)
# Btrig <- bpa
# idx0 <- idx
# dt <- date()
# 
# #########################################################
# # go fish
# for(i in vy[-length(vy)]){
#   ## i <- vy[-length(vy)][1]
#   print(i)
#   gc()
#   ay <- an(i)   # an is equivalent to as.numeric
#   cat(i, "\n")
#   vy0 <- 1:(ay-y0) 
#   sqy <- (ay-y0-nsqy+1):(ay-y0)
#   stk0 <- pstk[,vy0]
#   catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
#   ## note that vy0 is changing below so index is being updated
#   for (index_counter in 1:length(idx)){
#     idx0[[index_counter]] <- idx[[index_counter]][,vy0]
#     index(idx[[index_counter]])[,i] <- stock.n(pstk)[,i]*index.q(idx[[index_counter]])[,i] + 1
#   }
#   ##
#   qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
#   fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
#   rmodel4 <- ~ s(year, k=20)  
#   fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5) #, qmodel=qmod5)
#   stk0 <- stk0 + fit
#   # fwd control
#   fsq0 <- yearMeans(fbar(stk0)[,sqy])
#   dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
#   arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
#   ## ftrg.vec <- rep(ftrg, it) ## original
#   #refpt <- data.frame(ssb = 1, harvest = 1)
#   #ftrg.q <- hcr.nocheck.GFCM.f(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, refpt = refpt, Btrig = Btrig, Fmin = 0, Blim = blim, Bpa=bpa)
#   #ftrg.q <- hcr.nocheck(ssb(stk0)[, ac(an(i) - 1)], refpt = refpt, Ftar = ftrg, Btrig = bpa, Fmin = 0, Blim = blim)
#   #ftrg.vec <- an(ftrg.q)
#   #Bescape <- blim
#   catch.q <- yearMeans(landings(stk)[,"2014"]) * 0.91
#   arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
#   arr0[,,"val"] <- c(yearMeans(landings(stk0)[,sqy]), rep(an(catch.q),it))
#   #arr0[,,"min"] <- c(rep(NA, 2 * it), SSBtrg.vec)
#   arr0 <- aperm(arr0, c(2,3,1))
#   ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('catch', 'catch'), val=NA))
#   ctrl@trgtArray <- arr0
#   #ctrl@trgtArray <- arr0
#   ## 
#   stkTmp <- stf(stk0, 2)
#   stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
#   ## USING F
#   TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
#   # OM proj
#   ctrl@target <- ctrl@target[2,]
#   ## original was catch
#   ##ctrl@target[,"quantity"] <- "catch"
#   ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
#   ## original was catch
#   ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
#   pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
#   #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
# }
# 
# 
# ##
# 
# date()
# ANE_Opt6b_mse.pstk.GFCM_segreg_9perc_red_onCatch <- pstk
# PIL_Opt6b_mse.stk0.GFCM_segreg_9perc_red_onCatch <- stk0
# plot(PIL_Opt6b_mse.pstk.GFCM_segreg_9perc_red_onCatch)
# save(list = c("PIL_Opt6b_mse.pstk.GFCM_segreg_9perc_red_onCatch", "PIL_Opt6b_mse.stk0.GFCM_segreg_9perc_red_onCatch"), file = "Results/ANE/PIL_Opt6b_mse_segreg_9perc_red_250it_onCatch.RData")
# png("PIL_Opt6b_mse_Segreg__9perc_red_250it_onCatch.png", width=700, height=700)
# plot(PIL_Opt6b_mse.pstk.GFCM_segreg_9perc_red_onCatch)
# dev.off()
# 
# 


#------------------------------------------------------------------------------
# 7a. (10a) SCENARIO BASED ON TAC = Catch(2014) = 31085 tonnes
#------------------------------------------------------------------------------
# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
Fmsy <- 0.554
Ftar <- Fmsy
dt <- date()
trgyF <- 2020
trgyB <- 2018
Fmin <- 0.1
TrgtC <- 31085


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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  stk0 <- stk0 + fit
  # fwd control
  Csq0 <- catch(stk0)[,ac(ay-1)]
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  
  ctrg <- hcr.max.catch(ssb(stk0)[, ac(an(i) - 1)], Csq0=Csq0, refpt = data.frame(ssb = 1, harvest = 1), Btrig = Btrig, Fmin = 0.1, Blim = blim, TrgtC=TrgtC)
  ctrg.vec <- an(ctrg)
  
  # 
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- c(Csq0,ctrg.vec)
  arr0 <- aperm(arr0, c(2,3,1))
  
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c("catch", "catch"), val=NA))
  # ctrl <- fwdControl(data.frame(year=i, quantity='catch', val=NA),
  #                    trgtArray=array(c(rbind(NA, NA, 50000)), dim=c(1, 3, it),
  #                                    dimnames=list(1, c('min','val','max'), iter=seq(it))))
  #ctrl@trgtArray <- ctrl@trgtArray[1,,,drop=FALSE]
  ctrl@trgtArray <- arr0
  
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- TrgtC
  
  
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}


##
return(val)
date()
plot(pstk)


ANE_Opt7_mse.pstk.GFCM_segreg_bkptmeanSSB_TAC <- pstk
ANE_Opt7_mse.stk0.GFCM_segreg_bkptmeanSSB_TAC <- stk0
plot(ANE_Opt7_mse.pstk.GFCM_segreg_bkptmeanSSB_TAC)
save(list = c("ANE_Opt7_mse.pstk.GFCM_segreg_bkptmeanSSB_TAC", "ANE_Opt7_mse.stk0.GFCM_segreg_bkptmeanSSB_TAC"), file = "Results/ANE/ANE_Opt7_mse_Segreg_bkptmeanSSB_TAC_250it.RData")
# png("ANE_Opt7_mse_Segreg_bkptmeanSSB_TAC_250it.png", width=700, height=700)
# plot(ANE_Opt7_mse.pstk.GFCM_segreg_bkptmeanSSB_TAC)
# dev.off()
# #

#------------------------------------------------------------------------------
# 7b. (10b) SCENARIO BASED ON TAC = minimum Catch = 5000 tonnes
#------------------------------------------------------------------------------
# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
Fmsy <- 0.554
Ftar <- Fmsy
dt <- date()
trgyF <- 2020
trgyB <- 2018
Fmin <- 0.1
TrgtC <- 5000


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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  stk0 <- stk0 + fit
  # fwd control
  Csq0 <- catch(stk0)[,ac(ay-1)]
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  
  ctrg <- hcr.max.catch(ssb(stk0)[, ac(an(i) - 1)], Csq0=Csq0, refpt = data.frame(ssb = 1, harvest = 1), Btrig = Btrig, Fmin = 0.1, Blim = blim, TrgtC=TrgtC)
  ctrg.vec <- an(ctrg)
  
  # 
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- c(Csq0,ctrg.vec)
  arr0 <- aperm(arr0, c(2,3,1))
  
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c("catch", "catch"), val=NA))
  # ctrl <- fwdControl(data.frame(year=i, quantity='catch', val=NA),
  #                    trgtArray=array(c(rbind(NA, NA, 50000)), dim=c(1, 3, it),
  #                                    dimnames=list(1, c('min','val','max'), iter=seq(it))))
  #ctrl@trgtArray <- ctrl@trgtArray[1,,,drop=FALSE]
  ctrl@trgtArray <- arr0
  
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- TrgtC
  
  
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}


##
return(val)
date()
plot(pstk)


ANE_Opt7b_mse.pstk.GFCM_segreg_bkptmeanSSB_mincatchTAC <- pstk
ANE_Opt7b_mse.stk0.GFCM_segreg_bkptmeanSSB_mincatchTAC <- stk0
plot(ANE_Opt7_mse.pstk.GFCM_segreg_bkptmeanSSB_TAC)
save(list = c("ANE_Opt7b_mse.pstk.GFCM_segreg_bkptmeanSSB_mincatchTAC", "ANE_Opt7b_mse.stk0.GFCM_segreg_bkptmeanSSB_mincatchTAC"), file = "Results/ANE/ANE_Opt7b_mse_Segreg_bkptmeanSSB_mincatchTAC_250it.RData")
# png("ANE_Opt7b_mse_Segreg_bkptmeanSSB_mincatchTAC_250it.png", width=700, height=700)
# plot(ANE_Opt7b_mse.pstk.GFCM_segreg_bkptmeanSSB_mincatchTAC)
# dev.off()
# #

#------------------------------------------------------------------------------
# 7c. (10c) SCENARIO BASED ON TAC = catch at Fmsy = 15327.0 tonnes
#------------------------------------------------------------------------------
# Set up the Btrigger (in this case Bpa)
Btrig <- bpa
idx0 <- idx
Fmsy <- 0.554
Ftar <- Fmsy
dt <- date()
trgyF <- 2020
trgyB <- 2018
Fmin <- 0.1
TrgtC <- 15327


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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  stk0 <- stk0 + fit
  # fwd control
  Csq0 <- catch(stk0)[,ac(ay-1)]
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  
  ctrg <- hcr.max.catch(ssb(stk0)[, ac(an(i) - 1)], Csq0=Csq0, refpt = data.frame(ssb = 1, harvest = 1), Btrig = Btrig, Fmin = 0.1, Blim = blim, TrgtC=TrgtC)
  ctrg.vec <- an(ctrg)
  
  # 
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- c(Csq0,ctrg.vec)
  arr0 <- aperm(arr0, c(2,3,1))
  
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c("catch", "catch"), val=NA))
  # ctrl <- fwdControl(data.frame(year=i, quantity='catch', val=NA),
  #                    trgtArray=array(c(rbind(NA, NA, 50000)), dim=c(1, 3, it),
  #                                    dimnames=list(1, c('min','val','max'), iter=seq(it))))
  #ctrl@trgtArray <- ctrl@trgtArray[1,,,drop=FALSE]
  ctrl@trgtArray <- arr0
  
  stkTmp <- stf(stk0, 2)
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
  ## USING F
  TAC[,ac(ay+1)] <- TrgtC
  
  
  ctrl@target <- ctrl@target[2,]
  ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
  pstk <- fwd(pstk, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}


##
return(val)
date()
plot(pstk)


ANE_Opt7c_mse.pstk.GFCM_segreg_bkptmeanSSB_fmsycatchTAC <- pstk
ANE_Opt7c_mse.stk0.GFCM_segreg_bkptmeanSSB_fmsycatchTAC <- stk0
plot(ANE_Opt7c_mse.pstk.GFCM_segreg_bkptmeanSSB_fmsycatchTAC)
save(list = c("ANE_Opt7c_mse.pstk.GFCM_segreg_bkptmeanSSB_fmsycatchTAC", "ANE_Opt7c_mse.stk0.GFCM_segreg_bkptmeanSSB_fmsycatchTAC"), file = "Results/ANE/ANE_Opt7c_mse_Segreg_bkptmeanSSB_fmsycatchTAC_250it.RData")
# png("ANE_Opt7c_mse_Segreg_bkptmeanSSB_fmsycatchTAC_250it.png", width=700, height=700)
# plot(ANE_Opt7c_mse.pstk.GFCM_segreg_bkptmeanSSB_fmsycatchTAC)
# dev.off()
# #

#------------------------------------------------------------------------------
# 9. (3) Changing F by different percentages according to the year (2015 - 2018)
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.554
trgyF <- 2019 # you want your measures to stop at trgyF (2019) and be equal to status quo from then onaight away
dt <- date()
# construct a vector of percentage reductions to apply on Fsq
perc <- c(1, 0.91, 1, 0.91)

#########################################################
# go fish
for(i in vy[-length(vy)]){
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  ay <- an(i)   # an is equivalent to as.numeric
  trgVec <- ay
  perc <- c(1, 0.91, 1, 0.91)
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
  fmod4 <- ~ s(age, k = 5) + factor(year)
  rmodel4 <- ~s(year, k=20)
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod4, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- yearMeans(fbar(stk0)[,sqy])
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.red.f.year(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, trgyF=trgyF, trgVec, perc)
  ftrg.vec <- an(ftrg.q)
  Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
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
return(val)
date()
# with bescape
ANE_Opt9_mse.pstk.red.F.yr <- pstk
ANE_Opt9_mse.stk0.red.F.yr <- stk0
save(list = c("ANE_Opt9_mse.pstk.red.F.yr", "ANE_Opt9_mse.stk0.red.F.yr"), file = "Results/ANE/ANE_Opt9_mse_Segreg_red.F.yr_250it.RData")
png("ANE_Opt9_mse_Segreg_red.F.yr_250it.RData.png", width=700, height=700)
plot(ANE_Opt9_mse.pstk.red.F.yr)
dev.off()

#------------------------------------------------------------------------------
# 12. (12) ENFORCEMENT OF FISHING DAYS REDUCTION - REDUCTION OF 20% OF F IN 2015
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
  ftrg.q <- fbar(stk)[,"2015",,,,] * 0.8
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
ANE_Op12_mse.pstk.GFCM_segreg_20perc_red_250it <- pstk
ANE_Opt12_mse.stk0.GFCM_segreg_20perc_red_250it <- stk0
#plot(ANE_Op12_mse.pstk.GFCM_segreg_20perc_red_250it)
save(list = c("ANE_Op12_mse.pstk.GFCM_segreg_20perc_red_250it", "ANE_Opt12_mse.stk0.GFCM_segreg_20perc_red_250it"), file = "Results/ANE/ANE_Opt12_mse_segreg_20perc_red_250it.RData")
# png("ANE_Opt12_mse_Segreg__20perc_red_250it_250it.png", width=700, height=700)
# plot(ANE_Op12_mse.pstk.GFCM_segreg_20perc_red_250it)
# dev.off()




#------------------------------------------------------------------------------
# 13. (13) ENFORCEMENT OF FISHING DAYS REDUCTION - REDUCTION OF 50% OF F IN 2015
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
  ftrg.q <- fbar(stk)[,"2015",,,,] * 0.5 
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
ANE_Op13_mse.pstk.GFCM_segreg_9perc_red_250it <- pstk
ANE_Opt13_mse.stk0.GFCM_segreg_9perc_red_250it <- stk0
#plot(ANE_Opt13_mse.pstk.GFCM_segreg_9perc_red_250it)
save(list = c("ANE_Op13_mse.pstk.GFCM_segreg_9perc_red_250it", "ANE_Opt13_mse.stk0.GFCM_segreg_9perc_red_250it"), file = "Results/ANE/ANE_Opt13_mse_segreg_9perc_red_250it.RData")
# png("ANE_Opt13_mse_Segreg__9perc_red_250it_250it.png", width=700, height=700)
# plot(ANE_Opt13_mse.pstk.GFCM_segreg_9perc_red_250it)
# dev.off()

#------------------------------------------------------------------------------
# 14. (14) ENFORCEMENT OF FISHING DAYS REDUCTION - REDUCTION OF 90% OF F IN 2015
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
  ftrg.q <- fbar(stk)[,"2015",,,,] * 0.1 
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
ANE_Op14_mse.pstk.GFCM_segreg_90perc_red_250it <- pstk
ANE_Opt14_mse.stk0.GFCM_segreg_90perc_red_250it <- stk0
#plot(ANE_Opt13_mse.pstk.GFCM_segreg_9perc_red_250it)
save(list = c("ANE_Op14_mse.pstk.GFCM_segreg_90perc_red_250it", "ANE_Opt14_mse.stk0.GFCM_segreg_90perc_red_250it"), file = "Results/ANE/ANE_Opt14_mse_segreg_90perc_red_250it.RData")
# png("ANE_Opt14_mse_Segreg__90perc_red_250it.png", width=700, height=700)
# plot(ANE_Op14_mse.pstk.GFCM_segreg_90perc_red_250it)
# dev.off()


#------------------------------------------------------------------------------
# 15. (2) SIMULATION OF THE EFFECT OF CLOSED AREAS ON RECRUITMENT
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
blim <- min(ssb(stk))
bpa <- blim*2
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
dt <- date()

#--------------------------------------------
# Reduce the F of ages 0 and 1 in the harvest slot of the
# pstk for the first year of simulations (2016) according
# to Ernesto's email of 04/03/2017:

pstk.n <- pstk

harvest(pstk.n)[1,42,,,,] <- harvest(pstk.n)[1,42,,,,] * 0.9
harvest(pstk.n)[2,42,,,,] <- harvest(pstk.n)[2,42,,,,] * 0.95



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
  # define stock0 from pstk.n until the last populated year
  # pstk.n is at the beginning only populated into the future fro F
  # the rest is only 1975-2015 but as the loop progresses through the projection
  # years the object is populated with projected numbers
  stk0 <- pstk.n[,vy0]
  # add 1 to everything to avoid zeros
  catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
  ## note that vy0 is changing below so index is being updated
  for (index_counter in 1:length(idx)){
    idx0[[index_counter]] <- idx[[index_counter]][,vy0]
    index(idx[[index_counter]])[,i] <- stock.n(pstk.n)[,i]*index.q(idx[[index_counter]])[,i] + 1
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
  # update pstk.n with stkTmp
  pstk.n <- fwd(pstk.n, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #
}


return(val)
date()

# breakpoint mean(SSB)
ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec <- pstk.n
ANE_Opt15_mse.stk0.GFCM_segregmeanSSB_redrec <- stk0
plot(ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec)
save(list = c("ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec", "ANE_Opt15_mse.stk0.GFCM_segregmeanSSB_redrec"), file = "Results/ANE/ANE_Opt15_mse_SegregmeanSSB_redrec_250it.RData")
png("ANE_Opt15_mse_SegregmeanSSB_redrec_250it.png", width=700, height=700)
plot(ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec)
dev.off()

#------------------------------------------------------------------------------
# 15. (4) SIMULATION OF THE EFFECT OF CLOSED AREAS ON RECRUITMENT
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.554
trgyF <- 2016 # you want your measures to stop at trgyF (2016) and be equal to status quo from then on
dt <- date()

#--------------------------------------------
# Reduce the F of ages 0 and 1 in the harvest slot of the
# pstk for the first year of simulations (2016) according
# to Ernesto's email of 04/03/2017:

pstk.n <- pstk

harvest(pstk.n)[1,42,,,,] <- harvest(pstk.n)[1,42,,,,] * 0.9
harvest(pstk.n)[2,42,,,,] <- harvest(pstk.n)[2,42,,,,] * 0.95

Fsq15 <- yearMeans(fbar(pstk)[,39:41])
Fpstkn <- fbar(pstk.n[,"2016",,,,])

#########################################################
# go fish
for(i in vy[-length(vy)]){
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  ay <- an(i)   # an is equivalent to as.numeric
  trgVec <- ay
  cat(i, "\n")
  vy0 <- 1:(ay-y0) # data years (positions vector)
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
  stk0 <- pstk.n[,vy0]
  catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
  ## note that vy0 is changing below so index is being updated
  for (index_counter in 1:length(idx)){
    idx0[[index_counter]] <- idx[[index_counter]][,vy0]
    index(idx[[index_counter]])[,i] <- stock.n(pstk.n)[,i]*index.q(idx[[index_counter]])[,i] + 1
  }
  ##
  qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
  fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
  rmodel4 <- ~ s(year, k=20)  
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- Fsq15
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.fbar.mod01(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa, Fpstkn=Fpstkn)
  ftrg.vec <- an(ftrg.q)
  Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
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
  pstk.n <- fwd(pstk.n, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}

return(val)
date()

# breakpoint mean(SSB)
ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec <- pstk.n
ANE_Opt15_mse.stk0.GFCM_segregmeanSSB_redrec <- stk0
plot(ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec)
save(list = c("ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec", "ANE_Opt15_mse.stk0.GFCM_segregmeanSSB_redrec"), file = "Results/ANE/ANE_Opt15_mse_SegregmeanSSB_redrec_250it.RData")
png("ANE_Opt15_mse_SegregmeanSSB_redrec_250it.png", width=700, height=700)
plot(ANE_Opt15_mse.pstk.GFCM_segregmeanSSB_redrec)
dev.off()

#------------------------------------------------------------------------------
# 16. SIMULATION OF ALL F REDUCTIONS
#------------------------------------------------------------------------------

# Set up the Btrigger (in this case halfway between Blim and Bpa)
Btrig <- blim+((bpa-blim)/2)
idx0 <- idx
Ftar <- 0.554
trgyF <- 2019 # you want your measures to stop at trgyF (2016) and be equal to status quo from then on
dt <- date()
# construct a vector of percentage reductions to apply on Fsq
perc <- c(1, 0.91, 1, 0.91)
#--------------------------------------------
# Reduce the F of ages 0 and 1 in the harvest slot of the
# pstk for the first year of simulations (2016) according
# to Ernesto's email of 04/03/2017:

pstk.n <- pstk

harvest(pstk.n)[1,42,,,,] <- harvest(pstk.n)[1,42,,,,] * 0.9
harvest(pstk.n)[2,42,,,,] <- harvest(pstk.n)[2,42,,,,] * 0.95

Fsq15 <- yearMeans(fbar(pstk)[,39:41])
Fpstkn <- fbar(pstk.n[,"2016",,,,])

#########################################################
# go fish
for(i in vy[-length(vy)]){
  ## i <- vy[-length(vy)][1]
  print(i)
  gc()
  ay <- an(i)   # an is equivalent to as.numeric
  trgVec <- ay
  cat(i, "\n")
  vy0 <- 1:(ay-y0) # data years (positions vector)
  sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector)
  stk0 <- pstk.n[,vy0]
  catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
  ## note that vy0 is changing below so index is being updated
  for (index_counter in 1:length(idx)){
    idx0[[index_counter]] <- idx[[index_counter]][,vy0]
    index(idx[[index_counter]])[,i] <- stock.n(pstk.n)[,i]*index.q(idx[[index_counter]])[,i] + 1
  }
  ##
  qmod5 <- list(~s(age, k=5) + s(year, k=4), ~s(age, k=5) + s(year, k=4))
  fmod8 <- ~ s(age, k = 5) + s(year, k=18) + te(age, year, k = c(4,5)) 
  rmodel4 <- ~ s(year, k=20)  
  fit <- sca(stk0, FLIndices(idx0), fmodel=fmod8, srmodel=rmodel4, qmodel=qmod5)
  #   rmodel <- ~s(year, k=10)
  # #fmod <- ~ s(age, k=4) + s(year, k = 10) + te(age, year, k = c(5,3))
  # fmod <- ~ s(age, k = 5) + s(year, k=12) 
  # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
  stk0 <- stk0 + fit
  # fwd control
  fsq0 <- Fsq15
  dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  ## ftrg.vec <- rep(ftrg, it) ## original
  refpt <- data.frame(ssb = 1, harvest = 1)
  ftrg.q <- hcr.allF(ssb(stk0)[, ac(an(i) - 1)], Fsq0=fsq0, Ftar=Ftar, refpt = refpt, Btrig = Btrig, Fmin = 0.025, Blim = blim, Bpa=bpa,trgyF=trgyF,trgVec=trgVec,  Fpstkn=Fpstkn,perc=perc)
  ftrg.vec <- an(ftrg.q)
  Bescape <- blim
  arr0[,,"val"] <- c(fsq0, ftrg.vec)
  arr0[,,"min"] <- c(rep(NA, 2 * it))
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'f'), val=NA))
  ctrl@trgtArray <- arr0
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
  pstk.n <- fwd(pstk.n, ctrl=ctrl, sr=sr, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
  #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
}

return(val)
date()

# breakpoint mean(SSB)
ANE_Opt16_mse.pstk.GFCM_segregmeanSSB_allF <- pstk.n
ANE_Opt16_mse.stk0.GFCM_segregmeanSSB_allF <- stk0
plot(ANE_Opt16_mse.pstk.GFCM_segregmeanSSB_allF)
save(list = c("ANE_Opt16_mse.pstk.GFCM_segregmeanSSB_allF", "ANE_Opt16_mse.stk0.GFCM_segregmeanSSB_allF"), file = "Results/ANE/ANE_Opt16_mse_SegregmeanSSB_allF_250it.RData")
png("ANE_Opt16_mse_SegregmeanSSB_allF_250it.png", width=700, height=700)
plot(ANE_Opt16_mse.pstk.GFCM_segregmeanSSB_allF)
dev.off()
