
load(sname)
pstk  <- pstksave
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
if (!is.null(sp.closure) && sp.closure)
{
sel.change <- c(0.9,0.95,1.1,1.1,1.1)
harvest(pstk)[,vy]       <-  sweep (harvest(pstk)[,vy] , c(1,3:6) , sel.change , "*" )
}



#########################################################
# go fish!
escapeRuns <- numeric()
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
  sam0.ctrl@range["maxyear"]  <- iay-1
  
  continueRuns        <- which(!(1:dims(stk0)$iter) %in% escapeRuns)
  if (iay == iy)  res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,return.sam=T)
  if (iay > iy)   res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,starting.sam=res,return.sam=T)
  
  trouble <- data.frame(iter = 1:it , failure =  unlist(lapply(res , function(x) is.na(x))))
  trouble <- trouble[trouble$failure == T,]
  #if (dim(trouble)[1] > 0) for (ii in trouble$iter) res[[ii]] <-  FLSAM(iter(stk0,ii),FLIndices(lapply(idx0 , function(x) iter(x,ii))),sam0.ctrl)
  if (dim(trouble)[1] == 1){
    tres    <-  try(FLSAM(iter(stk0,trouble$iter),
                    FLIndices(lapply(idx0 , function(x) iter(x,trouble$iter))),
                    sam0.ctrl,silent=T))
    if(class(tres)=="try-error"){
      res[[trouble$iter]] <- new("FLSAM")
    } else {
      res[[trouble$iter]] <- tres
    }
  }
  if (dim(trouble)[1] > 1){
    resTrouble          <- FLSAM.MSE(iter(stk0,trouble$iter),
                                     FLIndices(lapply(idx0,function(x) iter(x,trouble$iter))),
                                     sam0.ctrl,return.sam=T)
    counter <- 1
    for(ii in trouble$iter){
      if(is.na(resTrouble[[counter]])){
        res[[ii]] <- new("FLSAM")
      } else {
        res[[ii]] <-  resTrouble[[counter]]
      }
      counter <- counter + 1
    }
  }
#

  for(ii in 1:it){
    if(!is.na(res[[ii]]@harvest[1,1,drop=T])){
      iter(stk0@harvest,ii) <- res[[ii]]@harvest
      iter(stk0@stock.n,ii) <- res[[ii]]@stock.n
    } else {
      escapeRuns <- sort(unique(c(escapeRuns,ii)))
    }
  }
### STF ON THE PERCIEVED STOCK TO PRODUCE AND ADVICE
# 
#  define the recruitment assumption to use in for the short term
mean_rec <- exp(yearMeans(log(rec(stk0)[,ac(iay-c(1:3))])))
if(!exists("srSTF")) srSTF <- fmle(as.FLSR(stk0, model="geomean"))
for (its in 1:it)  params(srSTF)["a",its] <- iter(mean_rec,its)
 
 
# compute the target to achieve in the advice year. 
# target has a slots  as a fwdControl object
#                   - val for its value 
#                   - quant to describe if we are manageing F, catches or biomass
#                   - rel to express if the quantity and values apply as a multiplier to the value in a given year
#
  
  # doing the advice separately for each iteration in order to allow for 
  # different types of target to be set for different iterations

  stkTmpAll <- stf(stk0, 2)
  require(doParallel)
  ncores <- detectCores()-1
  ncores <- ifelse(dims(stk0)$iter<ncores,dims(stk0)$iter,ncores)
  cla <- makeCluster(ncores) #set up nodes
  clusterEvalQ(cla,library(FLash))
  registerDoParallel(cla)

  target <- foreach(its = 1:it) %dopar% HCR(stk0[,,,,,its],mgt.target)
  for(its in 1:it){
    if (!is.na(addFred))  target[[its]]$val$y2     <- (1-addFred) * target[[its]]$val$y2
  }
  ctrl  <- foreach(its = 1:it) %dopar% fwdControl(data.frame(year=c(iay, iay+1), quantity= target[[its]]$quant, val=c((target[[its]]$val$y1), (target[[its]]$val$y2)), rel.year = target[[its]]$rel))
  stkTmp <- foreach(its = 1:it) %dopar% fwd(stkTmpAll[,,,,,its],ctrl=ctrl[[its]],sr=iter(srSTF,its),maxF=10)
  for (its in 1:it){
         TAC[,ac(iay+1),,,,its]   <- catch(stkTmp[[its]])[,ac(iay+1)]
         SSBad[,ac(iay+1),,,,its] <- ssb(stkTmp[[its]])[,ac(iay+1)]
         Fad[,ac(iay+1),,,,its]   <- fbar(stkTmp[[its]])[,ac(iay+1)]
  }
  if("doParallel" %in% (.packages()))
    detach("package:doParallel",unload=TRUE)
  if("foreach" %in% (.packages()))
    detach("package:foreach",unload=TRUE)
  if("iterators" %in% (.packages()))
    detach("package:iterators",unload=TRUE)

  stopCluster(cla)

#  start.time <- Sys.time()
#  for (its in 1:it)
#        {
#        target <-  HCR(iter(stk0,its) , mgt.target )
#                  # add extra reduction in F when part of the scenario
#        if (!is.na(addFred))  target$val$y2     <- (1-addFred) * target$val$y2
#
#        # create the control object
#        ctrl <- fwdControl(data.frame(year=c(iay, iay+1), quantity= target$quant, val=c((target$val$y1), (target$val$y2)), rel.year = target$rel))
#        # populate the iteration specific values
#        ## Short term forecast object 2 years of stk0
#        stkTmp              <- stkTmpAll[,,,,,its]
#        # project forward with the control you want and the SR rel you defined above, with residuals
#        stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=iter(srSTF,its)  ,maxF = 10)
#
#
#         # update objects storing the basis for the advice
#         TAC[,ac(iay+1),,,,its]   <- catch(stkTmp)[,ac(iay+1)]
#         SSBad[,ac(iay+1),,,,its] <- ssb(stkTmp)[,ac(iay+1)]
#         Fad[,ac(iay+1),,,,its]   <- fbar(stkTmp)[,ac(iay+1)]
#       }
#     difftime(Sys.time(),start.time)
 ### UPDATE THE OM BASED ON THE TAC ADVICE (with on year lag)
 dnms <- list(year=c(iay), c("min", "val", "max"),iter=1:it)
 arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
 arr0[,"val",] <- TAC[,ac(iay)]
 ctrlOM <- fwdControl(data.frame(year=c(iay), quantity=c('catch'), val=c(iterMeans(TAC[,ac(iay)])@.Data)))
 ctrlOM@trgtArray <- arr0
 # update pstk with stkTmp
  pstk <- fwd(pstk, ctrl=ctrlOM, sr=sr, sr.residuals = exp(sr.res[,ac(iay)]), sr.residuals.mult = TRUE , maxF = 10)


# save at each time step
restosave <- list(pstk = pstk,Fad=Fad,SSBad=SSBad,TAC=TAC)
save(restosave,file = paste0("./Results/",species,"/simres/",sc,"_",it,"its_",fy,".RData"))

}  # end of year loops

