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




make.arma.resid.lst <- function(arima.fit.list, age, years){

#    arima.fit.list= arima.fit.lst ;  age= 0 ; years= arima.fit.list

  ##----------------------------
  ## simulate structured recruitment residuals on log-scale
  ## fseparately for each iteration of the stock
  ## 
  ## arima.fit.list is a list of fitted arima models (each corresponding to an iteration
  ## standard deviation of the log recruitment  
  ## age is a numeric scalar age at recruitment
  ## years is a numeric vector of years for the residuals
  ## nit is a numeric scalar for the number of iterations
  ## Note: only works for ARMA components here so no drif, integration
  ##----------------------------
  nyr <- length(years)
  nit <- length(arima.fit.list)
  ## container matrix
  srDev.mat <- matrix(NA, nrow = nit, ncol = nyr)
  
  for(i in 1:nit){
  arima.fit <-arima.fit.list[[i]]
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
  
  ## sample
  
    srDev.mat[i, ] <- an(arima.sim(n = nyr,
                                   model = model.list,
                                   sd = sd.innov))    
  }
  ## create the FLQuant
  srDev <- FLQuant(an(t(srDev.mat)),
                   dimnames = list(year = years, age = age, iter = 1:nit))  
  return(srDev)
}











need.it <- F

if (need.it)
         {
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


## Changing F by different percentages according to the year (2015 - 2018)

hcr.red.f.year <- function (SSB, Fsq0, refpt, Btrig = 0.75, Fmin = 0.025, 
                            Blim = 0.25, Bpa=0.5, Ftar, trgyF, trgVec, perc) {
  
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  valF = NA
  
  
  if(trgyF-trgVec==4){valF = Fsq0*perc[1]} else {
    if(trgyF-trgVec==3){valF = Fsq0*perc[2]} else {
      if (trgyF-trgVec==2){valF = Fsq0*perc[3]} else {
        if(trgyF-trgVec==1){valF = Fsq0*perc[4]} else {
          valF = (fbar(stk0[,"2018",,,,]))*perc[4]}}}}
  
  #  if(trgyF-trgVec==4){valF = Fsq0*1} else {
  # if(trgyF-trgVec==3){valF = Fsq0*0.91} else {
  #   if (trgyF-trgVec==2){valF = Fsq0*1} else {
  #     if(trgyF-trgVec==1){valF = Fsq0*0.91} else {
  #       valF = (fbar(stk0[,"2018",,,,]))*perc}}}}
  
  
  return(valF)
}

environment(hcr.red.f.year) <- environment(hcr)


## TAC

hcr.max.catch <- function (SSB, Csq0, refpt, Btrig = 0.75, Fmin = 0.025, 
                           Blim = 0.25, Bpa=0.5, TrgtC) {
  
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  valC = NA
  
  
  valC = ifelse(Csq0 < TrgtC, Csq0, TrgtC)
  
  
  #  if(trgyF-trgVec==4){valF = Fsq0*1} else {
  # if(trgyF-trgVec==3){valF = Fsq0*0.91} else {
  #   if (trgyF-trgVec==2){valF = Fsq0*1} else {
  #     if(trgyF-trgVec==1){valF = Fsq0*0.91} else {
  #       valF = (fbar(stk0[,"2018",,,,]))*perc}}}}
  
  
  return(valC)
}

environment(hcr.max.catch) <- environment(hcr)

## ------------------------------------------------------
## ------------------------------------------------------

## Changing Fbar to include modifications of ages 0 and 1
## from 2016 included onwards

hcr.fbar.mod01 <- function (SSB, Fsq0, refpt, Btrig = 0.75, Fmin = 0.025, 
                            Blim = 0.25, Bpa=0.5, Ftar, Fpstkn) {
  
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  valF = NA
  
  
  if(trgVec==2015){valF = Fsq0} else {
    valF = Fpstkn}
  
  return(valF)
}

environment(hcr.fbar.mod01) <- environment(hcr)

## ------------------------------------------------------
## ------------------------------------------------------

## Applying all F reductions:
## 1.Changing Fbar to include modifications of ages 0 and 1 (10% and 5% respectivey) in 2016 +
## 9% reduction on Fbar in 2016
## 9% reduction on Fbar in 2018
## Maintain 2018 from 2019 onwards

hcr.allF <- function (SSB, Fsq0, refpt, Btrig = 0.75, Fmin = 0.025, 
                      Blim = 0.25, Bpa=0.5, Ftar,trgyF,trgVec, Fpstkn, perc) {
  
  if (Blim >= Btrig) 
    stop("Btrig must be greater than Blim")
  
  Btrig = refpt[1, "ssb"] * Btrig
  Blim = refpt[1, "ssb"] * Blim
  Fmin = refpt[1, "harvest"] * Fmin
  valF = NA
  
  
  if(trgyF-trgVec==4){valF = Fsq15} else {
    if(trgyF-trgVec==3){valF = Fpstkn*perc[2]} else {
      if (trgyF-trgVec==2){valF = Fpstkn*perc[3]} else {
        if(trgyF-trgVec==1){valF =  Fpstkn*perc[4]} else {
          valF = (fbar(pstk.n[,"2018",,,,]))*perc[4]}}}}
  
  return(valF)
}

environment(hcr.allF) <- environment(hcr)




               }




#-----------------------------------------------------------------------------------------------------------
# modified mcmc function to resample stock from the SAM assessment.   Thomas Brunel 2018
#-----------------------------------------------------------------------------------------------------------
    
# modified ADMB FLSAM's monteCarloStock function to get it to output the parameters    
    
monteCarloStock2 <- 
function (stck, sam, realisations,seed_number) 
{
  require(MASS)
  ctrl <- sam@control
  mcstck <- propagate(stck, iter = realisations)
  mcstck <- window(mcstck, start = range(sam)["minyear"], end = range(sam)["maxyear"])
  mcstck@stock.n[] <- NA
  mcstck@harvest[] <- NA
  set.seed(seed_number)
  random.param <<- mvrnorm(realisations, sam@params$value, sam@vcov)
  #save(random.param, file = file.path(run.dir, "random.param.RData"))
  n.states <- length(unique(ctrl@states[names(which(ctrl@fleets == 
                                                      0)), ]))
  yrs <- dims(sam)$minyear:dims(sam)$maxyear
  ages <- dims(sam)$age
  u <- random.param[, which(colnames(random.param) == "U")]
  ca <- random.param[, which(colnames(random.param) == "logCatch")]
  idxNs <- c(mapply(seq, from = seq(1, ncol(u), ages + n.states), 
                    to = seq(1, ncol(u), ages + n.states) + ages - 1, by = 1))
  idxFs <- c(mapply(seq, from = seq(1, ncol(u), ages + n.states) + 
                      ages, to = seq(1, ncol(u), ages + n.states) + n.states + 
                      ages - 1, by = 1))
  mcstck@stock.n[] <- exp(aperm(array(u[, idxNs], dim = c(realisations, 
                                                          ages, length(yrs))), perm = c(2, 3, 1)))
  mcstck@harvest[] <- exp(aperm(array(u[, idxFs], dim = c(realisations, 
                                                          n.states, length(yrs))), perm = c(2, 3, 1)))[ctrl@states[names(which(ctrl@fleets == 
                                                                                                                                 0)), ], , ]
  mcstck@catch[] <- exp(aperm(array(ca, dim = c(realisations, 
                                                length(yrs))), perm = c(2, 1)))
  return(mcstck)
}


# modified the TMB new FLSAM's monteCarloStock function to get it to output the parameters    

monteCarloStockTMB <-
function (stck, tun, sam, realisations, return.sam = FALSE, ...) 
{
    
#    stck <- ANCHOVY2 ; tun  <- ANCHOVY2.tun   ;  sam  <- ANCHOVY2.sam;    realisations <- 2
    
    require(doParallel)
    ctrl <- sam@control
    if (class(stck) == "FLStocks") {
        idxStck <- which.max(unlist(lapply(NSHs, function(x) {
            x@range["maxyear"]
        })))
        mcstck <- propagate(stck[[idxStck]], iter = realisations)
        mcstck <- window(mcstck, start = range(sam)["minyear"], 
            end = range(sam)["maxyear"])
        mcstck@stock.n[] <- NA
        mcstck@harvest[] <- NA
    }
    else {
        mcstck <- propagate(stck, iter = realisations)
        mcstck <- window(mcstck, start = range(sam)["minyear"], 
            end = range(sam)["maxyear"])
        mcstck@stock.n[] <- NA
        mcstck@harvest[] <- NA
    }
    sam@control@simulate <- TRUE
    sam@control@residuals <- FALSE
    object <- FLSAM(stck, tun, sam@control, return.fit = T)
    simdat <- replicate(realisations, c(object$data[names(object$data) != 
        "logobs"], object$obj$simulate(unlist(object$pl))["logobs"]), 
        simplify = FALSE)
    ncores <- detectCores() - 1
    ncores <- ifelse(realisations < ncores, realisations, ncores)
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library(FLSAM))
    clusterEvalQ(cl, library(stockassessment))
    registerDoParallel(cl)
    for (i in 1:realisations) {
        if (length(which(is.na(object$data$logobs))) > 0) 
            simdat[[i]]$logobs[which(is.na(object$data$logobs))] <- NA
    }
    runs <- foreach(i = 1:realisations) %dopar% try(sam.fitfast(simdat[[i]], 
        object$conf, object$pl, silent = T, ...))
    stopCluster(cl)
    if (return.sam) {
        resSAM <- list()
        for (i in 1:realisations) {
            if (!is.na(unlist(res[[i]]$sdrep)[1])) {
                resSAM[[i]] <- SAM2FLR(res[[i]], sam@control)
            }
            else {
                resSAM[[i]] <- NA
            }
        }
        resSAM <- as(resSAM, "FLSAMs")
    }
    if ("doParallel" %in% (.packages())) 
        detach("package:doParallel", unload = TRUE)
    if ("foreach" %in% (.packages())) 
        detach("package:foreach", unload = TRUE)
    if ("iterators" %in% (.packages())) 
        detach("package:iterators", unload = TRUE)
    if (!return.sam) {
        samRuns <- list()
        for (i in 1:realisations) {
            if (!is.na(unlist(runs[[i]]$sdrep)[1])) {
                samRuns[[i]] <- SAM2FLR(runs[[i]], sam@control)
                mcstck@stock.n[, , , , , i] <- samRuns[[i]]@stock.n
                mcstck@harvest[, , , , , i] <- samRuns[[i]]@harvest
            }
        }
    }
   
    
   pars  <- lapply( samRuns , function(x) params(x)$value)  
   random.param <<- matrix(unlist(pars)    , byrow= T, nrow = realisations ,dimnames =list(NULL, params(samRuns[[1]])$name ))
    
    
    if (return.sam) 
        ret <- resSAM
    if (!return.sam) 
        ret <- mcstck
    return(ret)
}
