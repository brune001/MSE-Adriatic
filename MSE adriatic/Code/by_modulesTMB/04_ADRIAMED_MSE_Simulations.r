
load(sname)
biol  <- biolsave
cat("running",run,"MSE\n")
 # sc <- "F.low"
 # sc <- "F.msy"
 # sc <- "Fmsy2025"
 # sc <- "C2014"
 # sc <- "Chistmin"
 # sc <- "C50red"
 # sc <- "GFCM.HCR" 
 # sc <- "Bpa.Fmsy2020"


# define the management options for this scenario
mgt.target <- management.scenarios[[sc]][["target"]]
HCR        <- management.scenarios[[sc]][["HCR"]]
sp.closure <- management.scenarios[[sc]][["spatial.closure"]]
addFred    <- management.scenarios[[sc]][["additionnal.F.reduction"]]

# if scenario with spatial closure, adjust the selection pattern to be applied for the OM
if (!is.null(sp.closure) && sp.closure){
  sel.change <- c(0.9,0.95,1.1,1.1,1.1)
  harvest(biol)[,vy]       <-  sweep (harvest(biol)[,vy] , c(1,3:6) , sel.change , "*" )
}


#-------------------------------------------------------------------------------
# Run the MSE
#-------------------------------------------------------------------------------

escapeRuns <- numeric()
for(i in vy[-length(vy)]){
  print(i)
  iay                   <- an(i)
  vy0                   <- 1:(iay-y0)
  sqy                   <- (iay-y0-nsqy+1):(iay-y0) # status quo years

  #-----------------------------------------------------------------------------
  # Create perceived stock from underlying 'true population' biol
  #-----------------------------------------------------------------------------
  stock                 <- biol[,vy0]
  # Apply uncertainty to the catch
  catch.n(stock)        <- catch.n(stock) * catch.dev[,vy0] # avoid zeros
  #- Correct for the cases where you have a zero-catch
  catch.n(stock)@.Data[catch.n(stock)<100] <- -1

  #-----------------------------------------------------------------------------
  # Create survey indices
  #-----------------------------------------------------------------------------

  idx0                  <- window(idx , end = iay-1)
  # Compute indices taking parameter uncertainty in Q and obs.var into account
  # NOTE THAT HISTORICAL PART OF THE TIME SERIES ARE NOT MODIFIED


  if (iay > iy) # don't do that for the first year because it correspond to the actual year the current assessment was carried out
  {
      for (idx.nb in 1:length(idx))
      {
        n.              <- name(idx[[idx.nb]])
        # define necessary dimensions
        idsty           <-range(idx[[idx.nb]])["minyear"]                       # first year of the survey
        idx0[[idx.nb]]  <- idx[[idx.nb]][,ac(idsty:(iay-1))]
        type.idx        <- sam.ctrl@fleets[n.]                                  # type of index (abund or ssb)
        age.surv        <- dimnames(index(idx0[[idx.nb]]))$age                  # age range
        time.surv       <- mean(c(range(idx0[[idx.nb]])[c("startf","endf")]))   # time of the year when survey is carried out

        if (type.idx == 2)
        {
        mod.n           <- stock.n(biol) * exp (- time.surv * ( harvest(biol) + m(biol)))
        index(idx0[[idx.nb]])[age.surv,ac(iy:(iay-1))] <- mod.n[age.surv,ac(iy:(iay-1))]*index.q(idx[[idx.nb]])[age.surv,ac(iy:(iay-1))]     # modelled index
        index(idx0[[idx.nb]])[age.surv,ac(iy:(iay-1))] <- index(idx0[[idx.nb]])[age.surv,ac(iy:(iay-1))] * index.var(idx[[idx.nb]])[age.surv,ac(iy:(iay-1))]   # adding uncertainty
        }

        if (type.idx == 3)
        {
        index(idx0[[idx.nb]])[,ac(iy:(iay-1))] <- ssb(biol)[,ac(iy:(iay-1))]*index.q(idx[[idx.nb]])[,ac(iy:(iay-1))]           # SSB index
        index(idx0[[idx.nb]])[,ac(iy:(iay-1))] <- index(idx0[[idx.nb]])[,ac(iy:(iay-1))] *  index.var(idx[[idx.nb]])[,ac(iy:(iay-1))]    # adding uncertainty
        }
      }
  }

  #-----------------------------------------------------------------------------
  # Do the assessment
  #-----------------------------------------------------------------------------
  sam0.ctrl             <- sam.ctrl
  sam0.ctrl@range["maxyear"]  <- iay-1
  
  continueRuns          <- which(!(1:dims(stock)$iter) %in% escapeRuns)
  if (iay == iy)  res   <- FLSAM.MSE(stock,idx0,sam0.ctrl,return.sam=T)
  if (iay > iy)   res   <- FLSAM.MSE(stock,idx0,sam0.ctrl,starting.sam=res,return.sam=T)
  
  trouble               <- data.frame(iter = 1:it , failure =  unlist(lapply(res , function(x) is.na(x))))
  trouble               <- trouble[trouble$failure == T,]
  if (dim(trouble)[1] == 1){
    tres                <-  try(FLSAM(iter(stock,trouble$iter),
                                FLIndices(lapply(idx0 , function(x) iter(x,trouble$iter))),
                                sam0.ctrl,silent=T))
    if(class(tres)=="try-error"){
      res[[trouble$iter]] <- new("FLSAM")
    } else {
      res[[trouble$iter]] <- tres
    }
  }
  if (dim(trouble)[1] > 1){
    resTrouble            <- FLSAM.MSE(iter(stock,trouble$iter),
                                       FLIndices(lapply(idx0,function(x) iter(x,trouble$iter))),
                                       sam0.ctrl,return.sam=T)
    counter           <- 1
    for(ii in trouble$iter){
      if(is.na(resTrouble[[counter]])){
        res[[ii]]     <- new("FLSAM")
      } else {
        res[[ii]]     <-  resTrouble[[counter]]
      }
      counter         <- counter + 1
    }
  }
  for(ii in 1:it){
    if(!is.na(res[[ii]]@harvest[1,1,drop=T])){
      iter(stock@harvest,ii) <- res[[ii]]@harvest
      iter(stock@stock.n,ii) <- res[[ii]]@stock.n
    } else {
      escapeRuns      <- sort(unique(c(escapeRuns,ii)))
    }
  }

  #-----------------------------------------------------------------------------
  # Run an Short Term Forecast
  #-----------------------------------------------------------------------------
  mean_rec <- exp(yearMeans(log(rec(stock)[,ac(iay-c(1:3))])))
  if(!exists("srSTF")) srSTF <- fmle(as.FLSR(stock, model="geomean"))
  #if(!exists("srSTF")) srSTF <- fmle(as.FLSR(stock, model="segreg"), fixed=list(b=mean(ssb(stock[,ac(1975:2016)],na.rm=T)))
  for (its in 1:it)  params(srSTF)["a",its] <- iter(mean_rec,its)
 
  # compute the target to achieve in the advice year.
  # target has a slots  as a fwdControl object
  #                   - val for its value
  #                   - quant to describe if we are manageing F, catches or biomass
  #                   - rel to express if the quantity and values apply as a multiplier to the value in a given year

  # doing the advice separately for each iteration in order to allow for 
  # different types of target to be set for different iterations

  stkTmpAll           <- stf(stock, 2)
  srSTFAll            <- list()
  for(its in 1:it)
    srSTFAll[[its]]   <- iter(srSTF,its)
    
  #- Parallel computation
  require(doParallel);  ncores <- detectCores()-1;  ncores <- ifelse(dims(stock)$iter<ncores,dims(stock)$iter,ncores); cla <- makeCluster(ncores); clusterEvalQ(cla,library(FLash)); registerDoParallel(cla)

  #- Set a target, evaluate if target can be reached, if not, set a target of F = 10 as a max
  target <- foreach(its = 1:it) %dopar% HCR(stock[,,,,,its],mgt.target)
  for(its in 1:it){
    if (!is.na(addFred))  target[[its]]$val$y2     <- (1-addFred) * target[[its]]$val$y2
  }
  ctrl                <- foreach(its = 1:it) %dopar% fwdControl(data.frame(year=c(iay, iay+1), quantity= target[[its]]$quant, val=c((target[[its]]$val$y1), (target[[its]]$val$y2)), rel.year = target[[its]]$rel))
  stkTmp              <- foreach(its = 1:it) %dopar% fwd(stkTmpAll[,,,,,its],ctrl=ctrl[[its]],sr=srSTFAll[[its]],maxF=10)
  if(ctrl[[1]]@target[2,"quantity"]=="ssb"){
    ssbAchieved       <- unlist(lapply(stkTmp,function(x){ssb(x)[,ac(iay+1)]}))
    ssbTarget         <- unlist(lapply(ctrl,function(y)y@target[2,"val"]))
    notAchieved       <- which(ssbAchieved < ssbTarget)
    if(length(notAchieved)>0){
      ctrl[notAchieved]  <- foreach(its = notAchieved) %dopar% fwdControl(data.frame(year=c(iay, iay+1), quantity= c(target[[its]]$quant[1],"f"), val=c((target[[its]]$val$y1), 0), rel.year = target[[its]]$rel))
    }
    stkTmp            <- foreach(its = 1:it) %dopar% fwd(stkTmpAll[,,,,,its],ctrl=ctrl[[its]],sr=srSTFAll[[its]],maxF=10)
  }

  #- Update objects
  for (its in 1:it){
    TAC[,ac(iay+1),,,,its]    <- catch(stkTmp[[its]])[,ac(iay+1)]
    SSBad[,ac(iay+1),,,,its]  <- ssb(stkTmp[[its]])[,ac(iay+1)]
    Fad[,ac(iay+1),,,,its]    <- fbar(stkTmp[[its]])[,ac(iay+1)]
  }
  if("doParallel" %in% (.packages()))
    detach("package:doParallel",unload=TRUE)
  if("foreach" %in% (.packages()))
    detach("package:foreach",unload=TRUE)
  if("iterators" %in% (.packages()))
    detach("package:iterators",unload=TRUE)
  stopCluster(cla)

  #-----------------------------------------------------------------------------
  # Update the biology
  #-----------------------------------------------------------------------------
  dnms                <- list(year=c(iay), c("min", "val", "max"),iter=1:it)
  arr0                <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,"val",]       <- TAC[,ac(iay)]
  ctrlOM              <- fwdControl(data.frame(year=c(iay), quantity=c('catch'), val=c(iterMeans(TAC[,ac(iay)])@.Data)))
  ctrlOM@trgtArray    <- arr0
  biol                <- fwd(biol, ctrl=ctrlOM, sr=sr, sr.residuals = exp(sr.res[,ac(iay)]), sr.residuals.mult = TRUE , maxF = 10)
 
  #- Apply process error per cohort age
  for(idxAge in 2:dims(biol)$age){
    biol@stock.n[idxAge,ac(i)] <- biol@stock.n[idxAge,ac(i)]*varProccError[idxAge-1,ac(an(ac(i))-(idxAge-1))]
  }

  #-----------------------------------------------------------------------------
  # save at each time step
  #-----------------------------------------------------------------------------
  restosave           <- list(biol = biol,Fad=Fad,SSBad=SSBad,TAC=TAC,trouble=trouble,stock=stock)
  save(restosave,file = paste0("./Results/",species,"/simres/",sc,"_",it,"its_",fy,"test.RData"))

}  # end of year loops