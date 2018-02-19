####################################################################
####################################################################
# MSE Functions
# Modified from Piera Carpi, Agurtzane Urtizberea Ijurco & Miguel Bernal
# 20 January 2016
# Changed for Anchovy MSE for GFCM 
# 
# Notes: new functions allow to test for different HCRs
# 1. Status quo (F constant as average last 3 years) 
# 2. HCR as in Rec GFCM 2013
# 3. HCR for gradual reduction of F between Bpa and Blim
# 4. HCR considering both F and B (if F < FMSY, then use B;
#    if F > FMSY, then use F. In any case, if B < Blim, F must be equal to 0)
# 5. HCR Emergency measure: reducing F by a fixed percentage
####################################################################
####################################################################

yearCumsum <- function(x, ...){
  x[] <- apply(x, c(1,3:6), cumsum)
  return(x)
}

yearDiffPerc <- function(x, ...){
  #x[,-1] <- x[,-1]/x[,-ncol(x)]-1
  x[,-1] <- x[,-1]/c(x[,1])-1
  x[,1] <- 0
  return(x)
}

as.table.FLQuants <- function(x){
  x <- mcf(x)
  df0 <- as.data.frame(do.call("cbind", lapply(x, c)))
  row.names(df0) <- dimnames(x)[[2]]
  df0
}

iterQuantiles <- function(x, prob=0.5, ...){
  return(apply(x, c(1:5), quantile, prob=prob, na.rm = FALSE))
}

an <- function(x, ...) as.numeric(x, ...)

cbind.FLQuant <- function(x, y){
  lst <- mcf(list(x, y))
  res <- lst[[1]]
  res[,dimnames(y)[[2]]] <- y
  res
} 

make.arma.resid <- function(arima.fit, age, years, nit){
  ##----------------------------
  ## simulate structured recruitment residuals on log-scale
  ## not used as moved to time-varying parameters
  ## 
  ## arima fit is a fitted arima model
  ## standard deviation of the log recruitment  
  ## age is a numeric scalar age at recruitment
  ## years is a numeric vector of years for the residuals
  ## nit is a numeric scalar for the number of iterations
  ## Note: only works for ARMA components here so no drif, integration
  ##----------------------------
  nyr <- length(years)
  ## coefficients
  theta <- coef(arima.fit)
  ## model set up for simulations, not fully general
  model.list <- list()
  if("ar1" %in% names(theta)){
    model.list$ar <- an(theta[grep("ar", names(theta))])
  }
  if("ma1" %in% names(theta)){
    model.list$ma <- an(theta[grep("ma", names(theta))])
  }
  ## standard deviations of the innovations
  sd.innov <- sqrt(arima.fit$sigma2)
  ## container matrix
  srDev.mat <- matrix(NA, nrow = nit, ncol = nyr)
  ## sample
  for(i in 1:nit){
    srDev.mat[i, ] <- an(arima.sim(n = nyr,
                                   model = model.list,
                                   sd = sd.innov))    
  }
  ## create the FLQuant
  srDev <- FLQuant(an(t(srDev.mat)),
                   dimnames = list(year = years, age = age, iter = 1:nit))  
  return(srDev)
}

## hcr that doesn't check the margins
## from FLBRP
hcr.nocheck <- function (SSB, refpt, Ftar = 0.8, Btrig = 0.75, Fmin = 0.025, 
                         Blim = 0.25) {
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  Ftar = refpt[1, "harvest"] * Ftar
  a = FLPar((Ftar - Fmin)/(Btrig - Blim))
  b = FLPar(Ftar - a * Btrig)
  val = qmax(qmin(sweep(sweep(SSB, 6, a, "*", check.margin = FALSE), 6, b, "+", check.margin = FALSE), 
                  Ftar), Fmin)
  return(val)
}
environment(hcr.nocheck) <- environment(hcr)

## HARVEST CONTROL RULE EQUAL TO GFCM MANAGEMENT 
## from FLBRP
# hcr.nocheck.GFCM <- function (SSB, Fsq0, refpt, Btrig = 0.75, Fmin = 0.025, 
#                          Blim = 0.25, Bpa=0.5) {
#   if (Blim >= Btrig) 
#     stop("Btrig must be greater than Blim")
#   Btrig = refpt[1, "ssb"] * Btrig
#   Blim = refpt[1, "ssb"] * Blim
#   Fmin = refpt[1, "harvest"] * Fmin
#   a = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), 1)
#   a[a < 0.01] <- 0.0
#   a = FLPar(a)
#   val = sweep(Fsq0, 6, a, "*", check.margin = FALSE) 
#   print(SSB)
#   print(a)
#   print(Fsq0)
#   return(val)
# }
# environment(hcr.nocheck.GFCM) <- environment(hcr)


## HARVEST CONTROL RULE EQUAL TO GFCM MANAGEMENT WITH FREE F
# CHECKED FROM SONIA!!! OK!!!
## from FLBRP
hcr.nocheck.GFCM.f <- function (SSB, Fsq0, refpt, Btrig = 0.75, Fmin = 0.025, 
                                Blim = 0.25, Bpa=0.5, Ftar, trgyF) {
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  a = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), ifelse(SSB > Bpa & Fsq0 < Ftar, 0, 1))
  b = ifelse(SSB < Btrig, 0, ifelse(SSB > Bpa & Fsq0 < Ftar, Ftar, 0))
  #a = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), 1)
  #a = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), 
  #           ifelse(Fsq0 > Fmsy, Fmsy/Fsq0, 1))
  a[a < 0.0] <- 0.0
  a = FLPar(a)
  val = qmax(qmin(sweep(sweep(Fsq0, 6, a, "*", check.margin = FALSE), 6, b, "+", check.margin = FALSE), 
                  Ftar), Fmin)
  valF_min = ifelse(val<Fmin, Fmin, val)
  valF = ifelse(Fsq0 > valF_min, Fsq0-((Fsq0-valF_min)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))), Fsq0+((valF_min-Fsq0)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))))
  return(valF)
}
environment(hcr.nocheck.GFCM.f) <- environment(hcr)
## ------------------------------------------------------
## ------------------------------------------------------


## HARVEST CONTROL RULE EQUAL TO GFCM MANAGEMENT WITH FREE F, REDUCTION APPLIED
## TO CATCHES 
hcr.nocheck.GFCM.csq <- function (SSB, Csq0, refpt, Btrig = 0.75, Fmin = 0.025, 
                                  Blim = 0.25, Bpa=0.5, trgyC, Cmin) {
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  Fmin = refpt[1, "harvest"] * Fmin
  Cmin = refpt[1, "catch"] * Cmin
  a = (as.numeric(SSB - Blim))/(Bpa - Blim)
  #b = ifelse(SSB < Btrig, 0, )
  #a = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), 1)
  #a = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), 
  #           ifelse(Fsq0 > Fmsy, Fmsy/Fsq0, 1))
  a[a < 0.00] <- 0.0
  a = FLPar(a)
  val = sweep(Csq0, 6, a, "*", check.margin = FALSE)
  valC_min = ifelse(val<Cmin, Cmin, val)
  valC = ifelse(Csq0 > valC_min, Csq0-((Csq0-valC_min)/(ifelse(trgyC - ay < 1, 1, trgyC - ay))), Csq0+((valC_min-Csq0)/(ifelse(trgyC - ay < 1, 1, trgyC - ay))))
  #print(SSB)
  #print(a)
  #print(b)
  #print(Csq0)
  return(valC_min)
}
environment(hcr.nocheck.GFCM.csq) <- environment(hcr)



## HARVEST CONTROL RULE BASED ON F AND B (if B < Bpa, B ==> Bpa; if F > FMSY, F ==> FMSY)
# CHECKED WITH SONIA: IT'S OK!!! 
# hcr.Ftar.Btar <- function (SSB, refpt, Ftar = 0.8, Btrig = 0.75, Fmin = 0.1, 
#                            Blim = 0.25, Fsq0 = fsq0, trgyF = trgyF) {
#   if (Blim >= Btrig) 
#     stop("Btrig must be greater than Blim")
#   Btrig = refpt[1, "ssb"] * Btrig
#   Blim = refpt[1, "ssb"] * Blim
#   Fmin = refpt[1, "harvest"] * Fmin
#   Ftar = refpt[1, "harvest"] * Ftar
#   
#   aBio = FLPar((SSB - Blim)/(Btrig - Blim))
#   aBio[aBio < 0.01] <- 0.0
#   valBio = sweep(Fsq0, 6, aBio, "*", check.margin = FALSE)
#   ifelse(valBio > Ftar, valBio-((valBio-Ftar)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))), valBio+((Ftar-valBio)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))))
#   
#   valF = ifelse(Fsq0 > Ftar, Fsq0-((Fsq0-Ftar)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))), Fsq0+((Ftar-Fsq0)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))))
#   valMin <- qmin(valBio, valF)
#   valMin_min <- ifelse(valMin<Fmin, Fmin, valMin)
#   return(valMin_min)
#   }
# environment(hcr.Ftar.Btar) <- environment(hcr)


hcr.Ftar.Btar <- function (SSB, refpt, Ftar = 0.8, Btrig = 0.75, Fmin = 0.1, 
                           Blim = 0.25, Bpa=0.5, Fsq0 = fsq0, trgyF = trgyF) {
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  Ftar = refpt[1, "harvest"] * Ftar
  aBio = ifelse(SSB < Btrig, (as.numeric(SSB - Blim))/(Bpa - Blim), ifelse(SSB > Bpa & Fsq0 < Ftar, 0, 1))
  bBio = ifelse(SSB < Btrig, 0, ifelse(SSB > Bpa & Fsq0 < Ftar, Ftar, 0))
  aBio[aBio < 0.0] <- 0.0
  aBio = FLPar(aBio)
  valBio = sweep(sweep(Fsq0, 6, aBio, "*", check.margin = FALSE), 6, bBio, "+", check.margin = FALSE)
  ifelse(valBio > Ftar, valBio-((valBio-Ftar)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))), valBio+((Ftar-valBio)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))))
  valF = ifelse(Fsq0 > Ftar, Fsq0-((Fsq0-Ftar)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))), Fsq0+((Ftar-Fsq0)/(ifelse(trgyF - ay < 1, 1, trgyF - ay))))
  valMin <- qmin(valBio, valF)
  valMin_min <- ifelse(valMin<Fmin, Fmin, valMin)
  return(valMin_min)
}
environment(hcr.Ftar.Btar) <- environment(hcr)
