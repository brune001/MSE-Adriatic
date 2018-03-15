


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

