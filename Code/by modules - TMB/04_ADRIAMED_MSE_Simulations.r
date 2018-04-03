
go_fish <- function(sc)                
{
 # sc <- "Fmsy"
 # sc <- "Fmsy2020"
 # sc <- "C2014"
 # sc <- "Chistmin"
 # sc <- "GFCM.HCR" 
 
 
 
 
load(file=paste0("./Results/",species,"/MSE_",assess.name,"_blank_objects_MSE_",".RData"))



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

  if (iay == iy) 
      {
      res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,return.sam=T)
      for(i in 1:it)
        {
        stk0@harvest <- res[[i]]@harvest
        stk0@stock.n <- res[[i]]@stock.n
        }
       }
  
  if (iay > iy) 
      { 
      res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,starting.sam=res,return.sam=T)
      for(i in 1:it)
        {
        stk0@harvest <- res[[i]]@harvest
        stk0@stock.n <- res[[i]]@stock.n
        }
      }


### STF ON THE PERCIEVED STOCK TO PRODUCE AND ADVICE
 
# compute the target to achieve in the advice year. 
# target has a slots  as a fwdControl object
#                   - val for its value 
#                   - quant to describe if we are manageing F, catches or biomass
#                   - rel to express if the quantity and values apply as a multiplier to the value in a given year
#
  target <-  HCR(stk0 , mgt.target )
  
  # create the control object
  ctrl <- fwdControl(data.frame(year=c(iay, iay+1), quantity= target$quant, val=c(mean(target$val$y1), mean(target$val$y2)), rel.year = target$rel))
  
  # populate the iteration specific values
#  dnms <- list(year=c(iay,(iay+1)), c("min", "val", "max"),iter=1:it)
  dnms <- list(1:2, c("min", "val", "max"),iter=1:it)  
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[1,"val",] <- c(target$val$y1)            #intermediate year in  the STF
  arr0[2,"val",] <- c(target$val$y2)      # advice year in the STF                                 # changed as above
  ctrl@trgtArray <- arr0
  
  ## Short term forecast object 2 years of stk0
  stkTmp <- stf(stk0, 2)
  # project forward with the control you want and the SR rel you defined above, with residuals
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr  ,Fmax = 5) 
  
  # update objects storing the basis for the advice
  TAC[,ac(iay+1)] <- catch(stkTmp)[,ac(iay+1)]
  SSBad[,ac(iay+1)] <- ssb(stkTmp)[,ac(iay+1)]
  Fad[,ac(iay+1)] <- fbar(stkTmp)[,ac(iay+1)]


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
 pstk <- fwd(pstk, ctrl=ctrlOM, sr=sr, sr.residuals = exp(sr.res[,ac(iay)]), sr.residuals.mult = TRUE , Fmax = 5) #
}  # end of year loops


res <- list(pstk = pstk,Fad=Fad,SSBad=SSBad,TAC=TAC)

save(res,file = paste0("./Results/",species,"/simres/",sc,"_",it"its_",fy,".RData"))
}