
load(fname)
cat("running",run,"MSE\n")
 # sc <- "F.low"
 # sc <- "F.msy"
 # sc <- "Fmsy2025"
 # sc <- "C2014"
 # sc <- "Chistmin"
 # sc <- "C5red"
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


### UPDATE THE OM BASED ON THE TAC ADVICE (with on year lag)
 
# instead we have to find the F to applied, based on the TAC, assuming catch = TAC
 dnms <- list(year=c(iay), c("min", "val", "max"),iter=1:it)
 arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
 arr0[,"val",] <- TAC[,ac(iay)]
 ctrlOM <- fwdControl(data.frame(year=c(iay), quantity=c('catch'), val=c(iterMeans(TAC[,ac(iay)])@.Data)))
 ctrlOM@trgtArray <- arr0
 # update pstk with stkTmp
 pstk <- fwd(pstk, ctrl=ctrlOM, sr=sr, sr.residuals = exp(sr.res[,ac(iay)]), sr.residuals.mult = TRUE , maxF = 10) #



# CREATE PERCEIVED STOCK
  stk0 <- pstk[,vy0]



### STF ON THE PERCIEVED STOCK TO PRODUCE AND ADVICE
#  define the recruitment assumption to use in for the short term
mean_rec <- exp(yearMeans(log(rec(stk0)[,ac(iay-c(1:3))])))
if(!exists("srSTF"))   srSTF <- fmle(as.FLSR(stk0, model="geomean"))
for (its in 1:it)  params(srSTF)["a",its] <- iter(mean_rec,its)
  
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
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=srSTF  ,maxF = 10) 
  
  # update objects storing the basis for the advice
  TAC[,ac(iay+1)] <- catch(stkTmp)[,ac(iay+1)]
  SSBad[,ac(iay+1)] <- ssb(stkTmp)[,ac(iay+1)]
  Fad[,ac(iay+1)] <- fbar(stkTmp)[,ac(iay+1)]




# save at each time step
restosave <- list(pstk = pstk,Fad=Fad,SSBad=SSBad,TAC=TAC)
save(restosave,file = paste0("./Results/",species,"/simres/",sc,".NO_SAM_",it,"its_",fy,".RData"))

}  # end of year loops