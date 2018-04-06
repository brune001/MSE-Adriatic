##############################################################################
#
#EJ(20150119)
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
#
#   NEW VERSION  MARCH 2018 : Thomas BRUNEL (WMR)
#   SAM is fully integrated (assessment in the loop, uncertainty on starting conditions and model parameters incorporated in the OM)
#
###############################################################################




#-------------------------------------------------------------------------------
# 1 : setting dimensions
#-------------------------------------------------------------------------------

it <- number.replicates.stock       # 250                      # iterations - should be 250
y0 <- range(stk)["minyear"]         # year zero (initial) = 1975
ny <- number.years.simulated        # number of years to project - Usually 20
# In order for this code to run iy = dy
dy <- range(stk)["maxyear"]         # data year
ay <- dy                            # assessment year
iy <- ay+1                          # initial projections year (also intermediate)
fy <- iy + ny -1                    # final year
py <- y0:dy                         # past years
vy <- ac(iy:fy)                     # future years
nsqy <- 3                           # number of SQ years upon which to average results

#mny <- 2020                         #2016 # min year to get to trg
#mxy <- 2020                         # 2016 # max year to get to trg
#
#-------------------------------------------------------------------------------
# 2 : Create stock object & use vcov for new realisations
#-------------------------------------------------------------------------------

sstk                            <- monteCarloStockTMB ( stk , ids , sam , it)



#--------------------------------------------------------------------------------------------
# 4 :  S/R
#--------------------------------------------------------------------------------------------

############################################################################
#      4A. hockey stick with break point at meanSSB over the entire time series:
############################################################################
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
#plot(sr)
sr.res <- residuals(sr)
#plot(sr.res)


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
sr.res <- make.arma.resid.lst(arima.fit.lst, age = 0, years = iy:(iy + ny-1) , rec.res)

# for plotting
pl.res <- window(residuals(sr),end = iy+ny-1 )
pl.res[,ac(iy:(iy + ny-1))] <- sr.res
pl.res <- exp(pl.res)
pl.res <- pl.res[,,,,,1:min(it,5)]
png(paste0("Results/",species,"/PLOTS/recruiment deviations.png"), width=700, height=700)
print(ggplot(pl.res, aes (x=year , y =data, colour = iter) ) + geom_line() +ggtitle("ARIMA recruitment deviations by iteration") + geom_vline(xintercept = iy-1)    +  geom_hline(yintercept = 1) )
dev.off()

# WITHOUT AUTOCORRELATION:
##sr.res1[] <- sample(c(residuals(sr)), ny*it, replace=TRUE)

# PLOTS
# Confidence interval plot
ssb <- apply(matrix(as.vector(stk@stock.n*stk@mat*stk@stock.wt), nrow=1+stk@range[[2]] ), 2, sum)       ### this does seem a little strange? why computing SSB this way (and not taking spawning time into account
ssb <- c(ssb(stk)@.Data)

segreg.meanssb  <- function(ab, ssb) log(ifelse(ssb >= mean(ssb), ab$a*mean(ssb), ab$a*ssb))
fit.meanssb <- eqsr_fit(stk,nsamp=2000, models = c("segreg.meanssb"))




plot.rec.res <- F

if (plot.rec.res == T)
{

png(paste0("Results/",species,"/PLOTS/segreg_bkptmeanSSB_confidence.png"), width=700, height=700)
eqsr_plot(fit.meanssb)
dev.off()
# GENERAL PLOTS
# Labeled plot (with years)
fitted <- as.data.frame(fitted(sr)); names(fitted)[7] <- "fitted"
Ssb    <- as.data.frame(ssb(sr))   ; names(Ssb)[7]    <- "ssb"
Obs    <- as.data.frame(rec(sr))   ; names(Obs)[7]    <- "obs"



df.sr.bktp.meanssb <- cbind(fitted[,-1] , data.frame(ssb = Ssb[,7] ,obs = Obs[,7]))
df.sr.bktp.meanssb.med <- aggregate(cbind(fitted,ssb,obs)~year , df.sr.bktp.meanssb , median)

png(paste0("Results/",species,"/PLOTS/SR.png"), width=700, height=700)
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

png(paste0("Results/",species,"/PLOTS/SSB.png"), width=700, height=700)
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

png(paste0("Results/",species,"/PLOTS/REC.png"), width=700, height=700)
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

 }







#######################################################################################################################################

# Fixed objects
TAC <- SSBad <- Fad <- FLQuant(NA, dimnames=list(TAC="all", year=vy, iter=1:it))


# CONDITIONNING FOR PROJECTING THE OPERATING MODEL DURING THE FIRST YEAR OF THE SIMS
# assumption here is that the catches in 2017 for each species are equal to the catch limit 
# established at the 2014 level
TAC[,(ac(iy))] <- catch(sstk)[,"2014"] #  WHAT IS A REALISTIC VALUE , needs to go somewhere else



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

logQ        <- random.param[,colnames(random.param)=="logFpar"]
logSdObs    <- random.param[,colnames(random.param)=="logSdLogObs"]


# check that the configuration is accoring to the convention used to write this script
Q<-unique(c(sam.ctrl@catchabilities ))
Q<-Q[Q>=0]
if(  max(Q)>(min(Q)+(length(Q)-1)))  stop("sorry, for this script to work, values in sam.ctrl@catchabilities must start from 0, increment by 1 and end at the number of parameters minus 1")
O<-unique(c(sam.ctrl@obs.vars ))
O<-O[O>=0]
if( max(O)>(min(O)+(length(O)-1)))  stop("sorry, for this script to work, values in sam.ctrl@obs.vars must start from 0, increment by 1 and end at the number of parameters minus 1")



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
     
      agesQ     <-  names(sam.ctrl@catchabilities[n.,]) [!is.na(sam.ctrl@catchabilities[n.,]) & sam.ctrl@catchabilities[n.,]!= -1]
      agesSdObs <-  names(sam.ctrl@obs.vars[n.,]) [!is.na(sam.ctrl@obs.vars[n.,]) & sam.ctrl@obs.vars[n.,] != -1]
      
      agesQq <- agesQ
      agesSdObsq  <-  agesSdObs
      
      if (fleet.type == 3 )    agesQq <-  agesSdObsq  <- "all"
      
      for (its in 1:it) 
          {
          # each iteration the estimated catchability in the slot index.q for each survey :
          iter( index.q(x)[agesQq,]   , its )     <- FLQuant( matrix(rep(exp(logQ[its,c(1+sam.ctrl@catchabilities[n.,agesQ])]),dim(index.q(x)[agesQq,])[2]),nrow = dim(index.q(x)[agesQq,])[1]), dimnames = dimnames(iter( index.q(x)[agesQq,] , its )))
          #  the observation SD in the slot index.var
          idx.dev<- FLQuant( matrix(rep(exp(logSdObs[its,(sam.ctrl@obs.vars[n.,agesSdObs]+1)]),dim(index.var(x)[agesSdObsq,])[2]),nrow = dim(index.var(x)[agesSdObsq,])[1]), dimnames = dimnames(iter( index.var(x)[agesSdObsq,] , its )))
          # we want to replace that by actual deviations that we can just multiply by the modeled index
          iter( index.var(x)[agesSdObsq,] , its )  <-  exp(apply(idx.dev , c(1,2) , function (x) rnorm(1,0,x)))
          # we don't want to apply a deviation on the historical data, so impose a value of 1
          index.var(x)[agesSdObsq,ac(y0idx:ay),,,,its]   <- 1
          
      

############ this is how it should be done if the survey was to be continued in the future, but apparently, this SSB survey stopped after 2012,
# so replace all the above by NAs          
          if (fleet.type == 3 ) index.var(x)[,ac(iy:fy),,,,its]   <- NA  # so that we generate NAs in the future ,
                                                   # we could just not include this survey in the part where new index values are generate (loops),
                                                   # but in case one day this survey constinue, maybe it is better to keep it in the framework.         
          }
          
       
    return(x)  
      })





## similarly generate a future errors for the catch matrix, deviations to be use as multipliers
      agesSdObs <-  names(sam.ctrl@obs.vars["catch unique",]) [!is.na(sam.ctrl@obs.vars["catch unique",]) & sam.ctrl@obs.vars["catch unique",]!= -1]
      catch.dev <- window(catch.n(sstk) , start=range(sstk)["minyear"],end=fy)

      for (its in 1:it) 
        {
        cdev <- matrix(rep(logSdObs[its,c(1+sam.ctrl@obs.vars["catch unique",agesSdObs])],dim(catch.dev)[2]),nrow = dim(catch.dev)[1])
        cdev <- exp(cdev)
        cdev <- FLQuant(cdev, dimnames = dimnames(iter( catch.dev , its )))
        cdev <- apply(cdev , c(1,2) , function (x) rnorm(1,0,x))
        cdev[,ac(y0:ay)]    <- 0
        cdev                <- exp(cdev) # just transform in the dimension in which it can be used as a multiplier on the catches
        iter(catch.dev,its) <- cdev
        }




pstksave <- pstk









