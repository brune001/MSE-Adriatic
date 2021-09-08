###### STECF EWG 17-09
# Split, 23 - 29 September 2017

# source("http://flr-project.org/R/instFLR.R")

library(FLCore)
library(FLEDA)
library(FLXSA)
library(FLAssess)
library(FLash)
library(FLBRP) 
library(FLSAM)
# library(SQLiteFL)
library(doBy)
library(reshape)

### ============================================================================
### Misc
### ============================================================================
data.source <- file.path("C:/Users/s.angelini/CNR/STECF_small pelagics 2017/LAST STOCK OBJ LIGAS/PIL_sa/")    #Data source, not code or package source!!!

### ============================================================================
### Prepare stock object for assessment
### ============================================================================
#Load object
SARDINE <- readFLStock(file.path(data.source, "Sar17_18.ndx"),no.discards=TRUE)

#Catch is calculated from: catch.wt * catch.n, however, the reported landings are
#normally different (due to SoP corrections). Hence we overwrite the calculate landings
SARDINE@catch           <- SARDINE@landings
units(SARDINE)[1:17]    <- as.list(c(rep(c("tonnes","thousands","tonnes"),4), 
                                     rep("NA",2),"f",rep("NA",2)))

#Set object details
SARDINE@name                              <- "SARDINE GSA 17_18"
range(SARDINE)[c("minfbar","maxfbar")]    <- c(1,3)
SARDINE                                   <- setPlusGroup(SARDINE,SARDINE@range["max"])

save(SARDINE, file='SARDINE.RData') # save the stock object

### ============================================================================
### Prepare index object for assessment
### ============================================================================
#Load and modify all numbers at age data
# load("C:/Users/PC09/Documents/Work/Adriatic/GFCM_2016/Anchovy/PC/Run_2idx_same2015/ANCHOVY2idx.RData)
SARDINE.tun   <- readFLIndices(file.path(data.source,"Sar17_18_fl.dat"))
SARDINE.tun   <- lapply(SARDINE.tun,function(x) {x@type <- "number"; return(x)})
# names(SARDINE.tun) <- c("Echo West", "Echo East", 'Echo East Biomass')
names(SARDINE.tun) <- c("Echo West", "Echo East")
SARDINE.tun[["Echo West"]]@range["plusgroup"] <- NA
SARDINE.tun[["Echo East"]]@range["plusgroup"] <- NA
# SARDINE.tun[["Echo East Biomass"]]@range["plusgroup"] <- NA

# Creates the FLIndex for the Echo East Biomass
dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2012,unit="unique",season="all",area="unique",iter=1))
#dmns     <- FLQuant(NA,dimnames=list(age="all",year=2003:2016,unit="unique",season="all",area="unique",iter=1))
biomSurf   <- FLIndex(index=dmns)
biomSurf@index[]  <- c(213410, 213477, 107902, 246593, 136907, 131542, 231809, 125031, 79372, 
                       89329)
range(biomSurf)["startf"] <- 0.75
range(biomSurf)["endf"] <- 0.83
type(biomSurf)  <- "biomass"
EchoEastBiomass <- FLIndices(EchoEastBiomass=biomSurf)

# ANCHOVY.tun <- FLIndices(c(ANCHO.tun[["Echo West"]], ANCHO.tun[["Echo East"]],EchoEastBiomass[["EchoEastBiomass"]])) 

SARDINE.tun <- FLIndices(list(SARDINE.tun[["Echo West"]], SARDINE.tun[["Echo East"]],EchoEastBiomass[["EchoEastBiomass"]]))
names(SARDINE.tun)<-c("Echo West","Echo East","Echo East Biomass")

save(SARDINE.tun, file='SARDINEtun.RData') # save the stock object

### ============================================================================
### Apply plusgroup to all data sets
### ============================================================================
pg <- 4

#- This function already changes the stock and landings.wts correctly
SARDINE <- setPlusGroup(SARDINE,pg)


######## ONLY ONE SURVEY
#ANCHOVY1idx.tun   <- readFLIndices(file.path(data.source,"Ane17_18_fl_1idx.dat"))
#ANCHOVY1idx.tun   <- lapply(ANCHOVY1idx.tun,function(x) {x@type <- "number"; return(x)})
#names(ANCHOVY1idx.tun) <- c("Echo WestEast")
#ANCHOVY1idx.tun[["Echo WestEast"]]@range["plusgroup"] <- NA
#ANCHOVY.tun[["Echo West18"]]@range["plusgroup"] <- NA


#################################
#### set up the control object
#################################

# library(FLSAM)

SARDINE.ctrl <- FLSAM.control(SARDINE,SARDINE.tun)

#Set the variances. Separate variance for recruitment and plus group
SARDINE.ctrl@logN.vars[]      <- c(1, rep(2,dims(SARDINE)$age-1))
#SARDINE.ctrl@f.vars["catch",] <- c(rep(1,3) ,rep(2,2))

SARDINE.ctrl@f.vars["catch",ac(0:2)] <- 11
SARDINE.ctrl@f.vars["catch",ac(3:4)] <- 12

#All fishing mortality states are free except
#oldest ages to ensure stablity
SARDINE.ctrl@states["catch",] <- seq(dims(SARDINE)$age)
SARDINE.ctrl@states["catch",ac(4)] <- 4

#Group observation variances of catches to ensure stability
#SARDINE.ctrl@obs.vars["catch",ac(0:1)]  <- 201
#SARDINE.ctrl@obs.vars["catch",ac(2:3)]  <- 202
#SARDINE.ctrl@obs.vars["catch",ac(4)]  <- 203
SARDINE.ctrl@obs.vars["catch",ac(0)]  <- 201 # new setting
SARDINE.ctrl@obs.vars["catch",ac(1)]  <- 202 # new setting
SARDINE.ctrl@obs.vars["catch",ac(2)]  <- 203 # new setting
SARDINE.ctrl@obs.vars["catch",ac(3:4)]  <- 204 # new setting

#SARDINE.ctrl@obs.vars["Echo West",ac(0:2)]    <- 204
#SARDINE.ctrl@obs.vars["Echo West",ac(3:4)]    <- 203
#SARDINE.ctrl@obs.vars["Echo East",ac(0:2)]    <- 205
#SARDINE.ctrl@obs.vars["Echo East",ac(3)]  <- 206
#SARDINE.ctrl@obs.vars["Echo East",ac(4)]  <- 207
SARDINE.ctrl@obs.vars["Echo West",ac(0)]    <- 205 # new setting
SARDINE.ctrl@obs.vars["Echo West",ac(1)]    <- 206 # new setting
SARDINE.ctrl@obs.vars["Echo West",ac(2:3)]    <- 207 # new setting
SARDINE.ctrl@obs.vars["Echo West",ac(4)]    <- 208 # new setting

SARDINE.ctrl@obs.vars["Echo East",ac(0)]    <- 209 # new setting
SARDINE.ctrl@obs.vars["Echo East",ac(1)]  <- 210 # new setting
SARDINE.ctrl@obs.vars["Echo East",ac(2)]  <- 211 # new setting


#Group catchability parametesr
#SARDINE.ctrl@catchabilities["Echo West",ac(0:4)]  <- c(rep(501,1),rep(502,2), rep(503, 1), rep(504, 1))
#SARDINE.ctrl@catchabilities["Echo East",ac(0:2)]  <- c(rep(501,1),rep(502,2)) # setting for age 0-2
#SARDINE.ctrl@catchabilities["Echo East",ac(0:4)]  <- c(rep(501,1),rep(502,2), rep(503, 1), rep(504, 1))
SARDINE.ctrl@catchabilities["Echo West",ac(0:4)]  <- c(rep(501,1),rep(502,1), rep(503, 1), rep(504, 1), rep(505,1)) # new setting
SARDINE.ctrl@catchabilities["Echo East",ac(0:2)]  <- c(rep(506,1),rep(507,1),rep(508,1)) # new setting


#Finalise
SARDINE.ctrl@name <- "Final Assessment"
SARDINE.ctrl <- update(SARDINE.ctrl)




### ============================================================================
### Run the assessment
### ============================================================================

#Perform assessment
SARDINE.sam <- FLSAM(SARDINE,SARDINE.tun,SARDINE.ctrl, run.dir="c:\\samtemp\\")
name(SARDINE.sam) <- "SARDINE GSA1718 2017"

#Update stock object
SARDINE      <- SARDINE + SARDINE.sam
SARDINE@stock <- computeStock(SARDINE)

# Save results
output.dir          <-  file.path("C:/Users/s.angelini/CNR/STECF_small pelagics 2017/LAST STOCK OBJ LIGAS/PIL_sa/results/")
save(SARDINE,SARDINE.tun,SARDINE.ctrl,SARDINE.sam,file=file.path(output.dir,paste(name(SARDINE),".RData",sep="")))


#ane15 <- load("C:/Users/PC09/Documents/Work/Adriatic/GFCM_SAF/Results/Anchovy_SAMassessmentGFCM2014/ANCHOVY1718_Run6.RData")
#ane15sam <- get(ane15[1])

##comparison <- FLStocks(ane15=ane15sam, ane16=ANCHOVY)
#tiff("PC/comparison15-16.tiff", width = 1800, height = 2000, res=200)
##plot(comparison) + theme_light(base_size = 16)
#dev.off()



###### Comparison between the assessment
# ane15 <- load("C:/Users/PC09/Documents/Work/Adriatic/GFCM_SAF/Results/Anchovy_SAMassessmentGFCM2014/ANCHOVY1718_Run6.RData")
# ane15sam <- get(ane15[1])

# ane16_2idx <- ANCHOVY
# ane16_2idx_tun <- ANCHOVY.tun
# ane16_2idx_sam <- ANCHOVY.sam
# ane16_2idx_ctrl <- ANCHOVY.ctrl

# comparison <- FLStocks(ane15=ane15sam, ane16=ANCHOVY, ane16_1idx=ANCHOVY1idx)
# tiff("PC/comparison15-16.tiff", width = 1800, height = 2000, res=200)
# plot(comparison) + theme_light(base_size = 16)
# dev.off()


######GRAFOVI-FUNCTIONS:###################
### Data exploration plots
### ======================================================================================================

#Ratio of mature and immature biomass              
mat.immat.ratio <- function(stk,...){
  default.args <- list(x=data~year, data=as.data.frame(bmass(stk)),
                       groups=expression(qname),
                       type="l",
                       main="Mature - Immature biomass ratio",
                       key=simpleKey(text=c("Mature", "Immature"), points=F, lines=T),
                       ylab="Relative biomass")
  arg.ins      <- list(...)
  call.args    <- default.args
  call.args[names(arg.ins)]  <- arg.ins
  do.call(xyplot,call.args)
}


# Plot time series of any slot in an FLR object (with class FLQuant) (added 18-03-2010 by NTH)
timeseries <- function(stck.,slot.){
  assign("stck.",stck.,envir=.GlobalEnv);assign("slot.",slot.,envir=.GlobalEnv);
  print(xyplot(data~year,data=slot(stck.,slot.),
               groups=age,
               auto.key=list(space="right",points=FALSE,lines=TRUE,type="b"),
               type="b",
               xlab="Year",ylab=paste("Time series of",slot.,ifelse(units(slot(stck.,slot.))=="NA","",paste("(",units(slot(stck.,slot.)),")",sep=""))),
               main=paste(stck.@name,"timeseries of",slot.),
               par.settings=list(superpose.symbol=list(pch=as.character(0:8),cex=1.25))))}

#Anomaly plots - primarily oriented towards weight at age anomalies, but should work for other things too
anom.plot <- function(x,...) {
  #Calculate anomalies
  means <- rowMeans(x)
  sds   <- apply(x,1,sd)
  anoms <- (x-means)/drop(sds@.Data)
  weight.anoms <- as.data.frame(t(drop(anoms@.Data)))
  #Generate plot
  yrs <- as.numeric(colnames(anoms))
  matplot(yrs,weight.anoms,pch=rownames(x),...)
  grid()
  smoother <- loess(unlist(weight.anoms) ~ rep(yrs,ncol(weight.anoms)),span=0.2)
  predict.yrs <- seq(min(yrs),max(yrs),length.out=100)
  smoothed <- predict(smoother,predict.yrs,se=TRUE)
  polygon(c(predict.yrs,rev(predict.yrs)),
          c(smoothed$fit+smoothed$se.fit*1.96,rev(smoothed$fit-smoothed$se.fit*1.96)),
          col="lightgrey")
  matpoints(yrs,weight.anoms,type="p",pch=rownames(x))
  lines(predict.yrs,smoothed$fit,lwd=2,col="black")
}

#Growth anomaly plots
growth.anom.plot <- function(y,...) {
  #Calculate growth
  flc <- FLCohort(y)
  growth <- apply(flc,2,diff)
  #Coerce back to an FLQuant
  x <- flc[1:(dims(flc)$age-1),]
  x@.Data[1:dims(x)$age,1:dims(x)$cohort,1,1,1,1] <-growth
  flq   <- flc2flq(x)
  flq <- flq[,2:(dims(flq)$year-1)]    #Drop the first and last year, where we don't have any information
  #Plot anomaly
  anom.plot(flq,...)
}

#Overlay time series by cohort and align by cohort
overlayTimeseries <- function(x,nyrs,ages){
  require(doBy)
  validObject(x)
  if(class(x)!="FLQuants") stop("Object is not an FLQuants")
  if(any(is.na(names(x))==T)) stop("Each FLQuant must have a name")
  require(reshape)
  lng   <- length(x)
  dmns  <- list()
  for(i in 1:lng)
    dmns[[i]] <- dimnames(x[[i]])
  
  ags   <- unique(unlist(lapply(dmns,function(x){return(x$age)})))
  if("all" %in% ags){ x <- x[-which(ags == "all")]; dmns <- dmns[-which(ags == "all")]}
  lng   <- length(x);
  idx   <- lapply(dmns,function(x){any(x$age %in% ages)})
  if(length(idx)>0){
    x   <- x[unlist(idx)]; dmns <- dmns[unlist(idx)]}
  lng   <- length(x)
  yrs   <- range(unlist(lapply(x,function(y){dims(y)[c("minyear","maxyear")]})))
  
  stk   <- data.frame()
  for(i in 1:lng)
    stk   <- rbind(stk,cbind(as.data.frame(rescaler(trim(window(x[[i]],start=(max(an(yrs))-(nyrs-1)),end=max(an(yrs))),age=c(max(ages[1],dims(x[[i]])$min):min(rev(ages)[1],dims(x[[i]])$max))))),qname=names(x)[i]))
  stk$track <- stk$year - stk$age
  
  stk <- orderBy(~age+qname+track,data=stk)
  xyplot(data ~ track,data=stk,groups=qname,type="l",
         prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},xlab="Cohort",ylab="Standardized residuals",
         auto.key=list(space="right",points=FALSE,lines=TRUE,type="l"),
         panel = panel.superpose,
         panel.groups = function(...) {
           res <- list(...)
           lng <- length(res$x)/nyrs
           for(i in 1:lng){
             panel.grid(v=-1,h=-1,lty=3)
             panel.xyplot(res$x[(nyrs*i-nyrs+1):(nyrs*i)],res$y[(nyrs*i-nyrs+1):(nyrs*i)],lty=i,type="l",col=res$col.line)
             panel.text(res$x[(nyrs*i-nyrs+1):(nyrs*i)],res$y[(nyrs*i-nyrs+1):(nyrs*i)],labels=stk$age[res$subscript[(nyrs*i-nyrs+1):(nyrs*i)]],col=res$col.line,cex=0.8)
           }
         },
         scales=list(alternating=1,y=list(relation="free",rot=0)))
  
}
#Plot all survey timeseries by age
surveyTimeseries <- function(x){
  if(class(x) != "FLIndices") stop("Object must be of class 'FLIndices'")
  validObject(x)
  
  lst <- lapply(x,index)
  names(lst) <- names(x)
  
  inds <- mcf(lapply(ANCHOVY.tun,index))
  # scale
  indsN01 <- lapply(inds, function(x){
    arr <- apply(x@.Data, c(1,3,4,5,6), scale)
    arr <- aperm(arr, c(2,1,3,4,5,6))
    # small trick to fix an apply "feature"
    dimnames(arr) <- dimnames(x)
    x <- FLQuant(arr)
    #x <- x[ac(2:10)]
  })
  
  indsN01 <- FLQuants(indsN01)
  # fine tune
  pfun <- function(x,y,...){
    panel.grid(v=-1,h=-1,col="grey",lty=3)
    panel.xyplot(x,y, ...)
  }
  assign("pfun", pfun, envir = .GlobalEnv)
  
  # plot
  print(xyplot(data~year|factor(age), data=indsN01, type="l",
               xlab="", ylab="", auto.key=list(space="right",columns=1,type="l",lines=T,points=F),panel=pfun,
               groups=qname,as.table=TRUE))
}

#Time series scaled between 0 and 1 and stacked
stacked.area.plot <- function(x,data,groups,...) {
  #Input arguments
  in.args <-  list(...)
  
  #Stack function
  stck <- function(x,y,subscripts,groups) {
    #Build data frame for splitting with NAs to zero
    dat <- data.frame(x,groups=groups[subscripts],y)
    dat$y[is.na(dat$y)] <- 0
    #Cumsums of y for a given x
    cumsums.l <- lapply(split(dat,dat$x),function(tmp) {
      tmp$cumsums <- cumsum(tmp$y)
      return(tmp)})
    dat   <- do.call(rbind,cumsums.l)
    return(dat)
  }
  
  #Panel function
  pfun <- function(x,y,subscripts,groups,...) {
    panel.grid(h=-1,v=-1)
    #Generate cumsums
    dat <- stck(x,y,subscripts,groups)
    #Prepare the colour vector
    grps   <- unique(dat$groups)
    pfun.args <- list(...)
    cols    <- rep(pfun.args$cols,length.out=length(grps))
    pfun.args["col"] <- NULL
    pfun.args["cols"] <- NULL
    #Now plot each group!
    for(i in 1:length(grps)) {        #For loops aren't sexy, but they allow us to move through the colours as well
      poly.dat  <- subset(dat,dat$groups==grps[i])
      poly.plot <- data.frame(x=c(poly.dat$x,rev(poly.dat$x)),
                              y=c(poly.dat$cumsums,rev(poly.dat$cumsums-poly.dat$y)))
      do.call(panel.polygon,c(list(x=poly.plot$x,y=poly.plot$y,col=cols[i]),pfun.args))
    }
  }
  
  #Key default definition, from the "lattice" book, Figure 5.6
  grps     <-  subset(data,select=groups)[,1]
  n.grps   <-  length(unique(grps))
  cols     <-  if(is.null(in.args$col)) {rainbow(n.grps)} else { in.args$col}
  key.default <- list(right = list(fun = draw.colorkey,
                                   args = list(key = list(col = rep(cols,length.out=n.grps),
                                                          at = (1:(n.grps+1))-0.5,
                                                          labels=list(labels=as.character(unique(grps)),at=1:n.grps)),
                                               draw = FALSE)))
  
  #Setup default arguments
  default.args <- list(x,data,groups=as.formula(paste("~",groups)),
                       panel=pfun,
                       prepanel=function(x,y,subscripts,groups,...) {
                         dat <- stck(x,y,subscripts,groups)
                         return(list(ylim=range(pretty(c(0,dat$cumsums)))))
                       },
                       cols=cols,
                       scales=list(alternating=1),
                       legend=key.default)
  
  #Add in optional input args, and do plot
  plot.args <- default.args
  plot.args[names(in.args)] <- in.args
  do.call(xyplot,plot.args)
}

#Set penality function so that we don't get any scientific notation
options("warn.FPU"=FALSE)

# bubbleRes <- function(x, ...){
  #require(FLCore)
  #res <- x@residuals
  #res <- tapply(res$std.res, list(res$age, res$year, res$fleet), unique)
  ## make FLQuant obj
  #flq <- FLQuant(res, dim=dim(res))
  #dimnames(flq)[1:3] <- dimnames(res)
  ## bubble plot of the residuals
  #bubbles(quant ~ year | unit, data=flq, ...)
#}

#Setup plots
pdf(file.path(output.dir,paste(name(SARDINE),".pdf",sep="")))
# png(file.path(output.dir,"figures - %02d.png"),units = "px", height=800,width=672, bg = "white")
### ===========

### ============================================================================
### Input data
### ============================================================================

# Plot the overlay of tuning series
print(overlayTimeseries(lapply(SARDINE.tun,index),nyrs=10,ages=0:4))

# Plot the overlay by year and age
print(surveyTimeseries(SARDINE.tun))

# Cohort plot
cohort.df <- as.data.frame(catch.n(SARDINE), cohort = TRUE)
ggplot(cohort.df, aes(x = age, y = data/1000)) +
  geom_line() + geom_point() + facet_wrap( ~ cohort) + ylab("No. of Ind * 10^3") + xlab("Age") +
  theme(legend.position = "none") + theme_bw() + ggtitle("By-cohort year (yearclass)")

# Plot catch at age versus each other (internal consistency)
catchn.fli <- FLIndex(FLQuant(NA, dimnames=list(year = dimnames(catch.n(SARDINE))$year, age = dimnames(catch.n(SARDINE))$age)))
index(catchn.fli) <- catch.n(SARDINE)
plot(catchn.fli, type="internal", main='Internal consistency - catch')

#Plot survey index versus each other (internal consistency)
plot(SARDINE.tun[["Echo West"]],type="internal", main='Internal consistency - Echo West')
plot(SARDINE.tun[["Echo East"]],type="internal", main='Internal consistency - Echo East')

# Plot the proportion of catch and weight in numbers and weight to see if the catch is representative for the stock build-up
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE@catch.n)),groups="age",main="Proportion of Catch numbers at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE@stock.wt)),groups="age",main="Proportion of Stock weight at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE@catch.wt)),groups="age",main="Proportion of Catch weight at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))

# Plot the proportion of catch in numbers in the indices to see if the indices are having specific yearclass trends
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE.tun[["Echo West"]]@index)),groups="age",main="Proportion of Acoustic index at age - Echo West",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE.tun[["Echo East"]]@index)),groups="age",main="Proportion of Acoustic index at age - Echo East",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))
# print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE.tun[["Echo East Biomass"]]@index)),groups="age",main="Proportion of Acoustic index at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))

# Plot the proportion of natural mortality
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE@m)),groups="age",main="Proportion of natural at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))

# Plot the time series of weight in the stock and catch in the stock
timeseries(window(SARDINE,1975,range(SARDINE)["maxyear"]),slot="stock.wt")
timeseries(window(SARDINE,1975,range(SARDINE)["maxyear"]),slot="catch.wt")
timeseries(window(SARDINE,1975,range(SARDINE)["maxyear"]),slot="harvest")
timeseries(window(SARDINE,1975,range(SARDINE)["maxyear"]),slot="mat")
timeseries(window(SARDINE,1975,range(SARDINE)["maxyear"]),slot="m")

# Plot the time series of the surveys
timeseries(SARDINE.tun[["Echo West"]],slot="index")
timeseries(SARDINE.tun[["Echo East"]],slot="index")
# timeseries(SARDINE.tun[["Echo East Biomass"]],slot="index")

#Time series of west by cohort
west.by.cohort      <- as.data.frame(FLCohort(window(SARDINE@stock.wt,1980,range(SARDINE)["maxyear"])))
west.by.cohort      <- subset(west.by.cohort,!is.na(west.by.cohort$data))
west.by.cohort$year <- west.by.cohort$age + west.by.cohort$cohort
west.cohort.plot    <- xyplot(data~year,data=west.by.cohort,
                              groups=cohort,
                              auto.key=list(space="right",points=FALSE,lines=TRUE,type="b"),
                              type="b",
                              xlab="Year",ylab="Weight in the stock (kg)",
                              main=paste(SARDINE@name,"Weight in the stock by cohort"),
                              par.settings=list(superpose.symbol=list(pch=as.character(unique(west.by.cohort$cohort)%%10),cex=1.25)),
                              panel=function(...) {
                                panel.grid(h=-1,v=-1)
                                panel.xyplot(...)
                              })
print(west.cohort.plot)

### ============================================================================
### Model fit
### ============================================================================

#Survey fits
# residual.diagnostics(SARDINE.sam)

#Survey fits
residual.diagnostics9<-function (x, title = x@name) 
{
  index.res <- x@residuals
  index.res$obs <- exp(index.res$log.obs)
  index.res$mdl <- exp(index.res$log.mdl)
  index.res.l <- split(index.res, list(index.res$age, index.res$fleet), 
                       drop = TRUE)
  oldpar <- par(mfrow = c(3, 2), las = 0, oma = c(0, 0, 3, 
                                                  0), mgp = c(1.75, 0.5, 0), mar = c(3, 3, 2.5, 1), cex.main = 1, 
                tck = -0.01)
  for (i in 1:9) {
    #for (i in 1:length(index.res.l)) {
    ind.age <- unlist(strsplit(names(index.res.l[i]), "\\."))
    ttl.fmt <- ifelse(x@control@fleets[ind.age[2]] %in% c(3, 
                                                          4), "Diagnostics - %s", "Diagnostics - %s, age %s")
    ttl <- sprintf(ttl.fmt, ind.age[2], ind.age[1])
    ttl <- paste(title, ttl)
    idx.rng <- range(c(index.res.l[[i]]$obs, index.res.l[[i]]$mdl), 
                     na.rm = TRUE)
    idx.lim <- range(c(0, idx.rng))
    idx.exp <- floor(log10(max(pretty(idx.lim)))/3) * 3
    idx.div <- 10^idx.exp
    idx.label <- ifelse(idx.exp > 1, paste("values ", "[", 
                                           "10^", idx.exp, "]", sep = ""), paste("values "))
    index.res.l[[i]]$obs <- index.res.l[[i]]$obs/idx.div
    index.res.l[[i]]$mdl <- index.res.l[[i]]$mdl/idx.div
    idx.rng <- idx.rng/idx.div
    idx.lim <- idx.lim/idx.div
    idx.min <- min(c(index.res.l[[i]]$obs, index.res.l[[i]]$mdl))
    idx.max <- max(c(index.res.l[[i]]$obs, index.res.l[[i]]$mdl))
    idx.lim <- c(idx.min/1.15, ceiling(idx.max/0.12))
    res.range <- abs(range(index.res.l[[i]]$std.res))
    res.lim <- c(-max(res.range), max(res.range))
    plot(obs ~ year, index.res.l[[i]], log = "y", ylim = idx.lim, 
         xlab = "Year", ylab = idx.label, pch = 16)
    points(mdl ~ year, index.res.l[[i]], pch = 4)
    points(mdl ~ year, index.res.l[[i]], type = "l")
    legend("topleft", c("Observed", "Fitted"), pch = c(16, 
                                                       4), lty = c(NA, 1), horiz = TRUE)
    title("a) Observed and fitted values time series")
    plot(obs ~ mdl, index.res.l[[i]], log = "xy", ylim = idx.lim, 
         xlim = idx.lim, pch = 16, ylab = paste("Observed ", 
                                                idx.label, sep = ""), xlab = paste("Fitted ", 
                                                                                   idx.label, sep = ""))
    abline(0, 1, col = "black")
    legend("topleft", "1:1 line", lty = 1, horiz = TRUE)
    title("b) Observed vs fitted values")
    plot(std.res ~ year, index.res.l[[i]], ylim = res.lim, 
         ylab = "Standardised Residuals", xlab = "Year")
    points(std.res ~ year, index.res.l[[i]], type = "h")
    points(std.res ~ year, index.res.l[[i]], pch = 19, cex = 0.75)
    abline(h = 0)
    title("c) Standardised residuals over time")
    plot(std.res ~ mdl, index.res.l[[i]], ylab = "Standardised Residuals", 
         ylim = res.lim, xlim = idx.lim, pch = 19, cex = 0.75, 
         xlab = paste("Fitted ", idx.label, sep = ""), log = "x")
    plt.dat <- index.res.l[[i]]
    abline(h = 0)
    smoother <- loess.smooth(plt.dat$mdl, plt.dat$std.res)
    lines(smoother, col = "red")
    title("d) Tukey-Anscombe plot")
    qqnorm(index.res.l[[i]]$std.res, ylim = res.lim, xlab = "Quantiles of the Normal Distribution", 
           ylab = "Standardised Residuals", pch = 19, main = "")
    qqline(index.res.l[[i]]$std.res, col = "red")
    abline(0, 1, lty = 2)
    legend("topleft", c("qqline", "1:1 line"), lty = c(1, 
                                                       2), col = c("red", "black"), horiz = TRUE)
    title("e) Normal Q-Q plot")
    acf(as.ts(index.res.l[[i]]$std.res), ylab = "ACF", xlab = "Lag (yrs)", 
        type = c("partial"), ci.col = "black", main = "")
    legend("topright", legend = c("95% Conf. Int."), lty = c(2), 
           pch = c(NA), horiz = TRUE, box.lty = 0)
    title("f) Autocorrelation of Residuals")
    title(main = ttl, outer = TRUE)
  }
  par(oldpar)
}
###


residual.diagnostics10<-function (x, title = x@name) 
{
  index.res <- x@residuals
  index.res$obs <- exp(index.res$log.obs)
  index.res$mdl <- exp(index.res$log.mdl)
  index.res.l <- split(index.res, list(index.res$age, index.res$fleet), 
                       drop = TRUE)
  oldpar <- par(mfrow = c(3, 2), las = 0, oma = c(0, 0, 3, 
                                                  0), mgp = c(1.75, 0.5, 0), mar = c(3, 3, 2.5, 1), cex.main = 1, 
                tck = -0.01)
  for (i in 10:length(index.res.l)) {
    #for (i in 1:length(index.res.l)) {
    ind.age <- unlist(strsplit(names(index.res.l[i]), "\\."))
    ttl.fmt <- ifelse(x@control@fleets[ind.age[2]] %in% c(3, 
                                                          4), "Diagnostics - %s", "Diagnostics - %s, age %s")
    ttl <- sprintf(ttl.fmt, ind.age[2], ind.age[1])
    ttl <- paste(title, ttl)
    idx.rng <- range(c(index.res.l[[i]]$obs, index.res.l[[i]]$mdl), 
                     na.rm = TRUE)
    idx.lim <- range(c(0, idx.rng))
    idx.exp <- floor(log10(max(pretty(idx.lim)))/3) * 3
    idx.div <- 10^idx.exp
    idx.label <- ifelse(idx.exp > 1, paste("values ", "[", 
                                           "10^", idx.exp, "]", sep = ""), paste("values "))
    index.res.l[[i]]$obs <- index.res.l[[i]]$obs/idx.div
    index.res.l[[i]]$mdl <- index.res.l[[i]]$mdl/idx.div
    idx.rng <- idx.rng/idx.div
    idx.lim <- idx.lim/idx.div
    idx.min <- min(c(index.res.l[[i]]$obs, index.res.l[[i]]$mdl))
    idx.max <- max(c(index.res.l[[i]]$obs, index.res.l[[i]]$mdl))
    idx.lim <- c(idx.min/1.15, ceiling(idx.max/0.12))
    res.range <- abs(range(index.res.l[[i]]$std.res))
    res.lim <- c(-max(res.range), max(res.range))
    plot(obs ~ year, index.res.l[[i]], log = "y", ylim = idx.lim, 
         xlab = "Year", ylab = idx.label, pch = 16)
    points(mdl ~ year, index.res.l[[i]], pch = 4)
    points(mdl ~ year, index.res.l[[i]], type = "l")
    legend("topleft", c("Observed", "Fitted"), pch = c(16, 
                                                       4), lty = c(NA, 1), horiz = TRUE)
    title("a) Observed and fitted values time series")
    plot(obs ~ mdl, index.res.l[[i]], log = "xy", ylim = idx.lim, 
         xlim = idx.lim, pch = 16, ylab = paste("Observed ", 
                                                idx.label, sep = ""), xlab = paste("Fitted ", 
                                                                                   idx.label, sep = ""))
    abline(0, 1, col = "black")
    legend("topleft", "1:1 line", lty = 1, horiz = TRUE)
    title("b) Observed vs fitted values")
    plot(std.res ~ year, index.res.l[[i]], ylim = res.lim, 
         ylab = "Standardised Residuals", xlab = "Year")
    points(std.res ~ year, index.res.l[[i]], type = "h")
    points(std.res ~ year, index.res.l[[i]], pch = 19, cex = 0.75)
    abline(h = 0)
    title("c) Standardised residuals over time")
    plot(std.res ~ mdl, index.res.l[[i]], ylab = "Standardised Residuals", 
         ylim = res.lim, xlim = idx.lim, pch = 19, cex = 0.75, 
         xlab = paste("Fitted ", idx.label, sep = ""), log = "x")
    plt.dat <- index.res.l[[i]]
    abline(h = 0)
    #smoother <- loess.smooth(plt.dat$mdl, plt.dat$std.res)
    #lines(smoother, col = "red")
    title("d) Tukey-Anscombe plot")
    qqnorm(index.res.l[[i]]$std.res, ylim = res.lim, xlab = "Quantiles of the Normal Distribution", 
           ylab = "Standardised Residuals", pch = 19, main = "")
    qqline(index.res.l[[i]]$std.res, col = "red")
    abline(0, 1, lty = 2)
    legend("topleft", c("qqline", "1:1 line"), lty = c(1, 
                                                       2), col = c("red", "black"), horiz = TRUE)
    title("e) Normal Q-Q plot")
    acf(as.ts(index.res.l[[i]]$std.res), ylab = "ACF", xlab = "Lag (yrs)", 
        type = c("partial"), ci.col = "black", main = "")
    legend("topright", legend = c("95% Conf. Int."), lty = c(2), 
           pch = c(NA), horiz = TRUE, box.lty = 0)
    title("f) Autocorrelation of Residuals")
    title(main = ttl, outer = TRUE)
  }
  par(oldpar)
}

#output.dir          <-  file.path("C:/Users/s.angelini/CNR/STECF_small pelagics 2017/SARDINE/STECF_gfcm_2idx/results/other plots/")
#png(file.path(output.dir,"figures - %02d.png"),units = "px", height=800,width=672, bg = "white")

# RESIDUAL PLOTS:
residual.diagnostics9(SARDINE.sam)
residual.diagnostics10(SARDINE.sam)
# dev.off()

#Bubble plot of the residuals
# bubbleRes(SARDINE.sam, bub.scale=5)
# Log catchability residuals
bubbles(age ~ year|fleet, data = residuals(SARDINE.sam) , main = "Log catchability residuals")

# Plot the harvest pattern at age as a proportion over time & the stock.n as proportion over time
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE@harvest)),groups="age",main="Proportion of harvest pressure at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))
print(stacked.area.plot(data~year| unit, as.data.frame(pay(SARDINE@stock.n)),groups="age",main="Proportion of Stock numbers at age",ylim=c(-0.01,1.01),xlab="years",col=gray(9:0/9)))


#Plot result
#png('Fig38.png',units = "px", height=800,width=672, bg = "white")
print(plot(SARDINE.sam,futureYrs=T))
#dev.off()

#Plot uncertainties as a function of time
#png('Fig39.png',units = "px", height=800,width=672, bg = "white")
CV.yrs <- ssb(SARDINE.sam)$year
CV.dat <- cbind(SSB=ssb(SARDINE.sam)$CV,
                Fbar=fbar(SARDINE.sam)$CV,Rec=rec(SARDINE.sam)$CV)
matplot(CV.yrs,CV.dat,type="l",ylim=range(pretty(c(0,CV.dat))),yaxs="i",
        xlab="Year",ylab="CV of estimate",main="Uncertainties of key parameters")
legend("topleft",legend=colnames(CV.dat),lty=1:5,col=1:6,bty="n")
#dev.off()

#Plot catchabilities values
catch <- catchabilities(SARDINE.sam)
print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
             scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
             type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
             subset=fleet %in% c("Echo West"),
             main="Survey catchability parameters",ylab="Catchability",xlab="Age"))

#png('Fig41.png',units = "px", height=800,width=672, bg = "white")
catch <- catchabilities(SARDINE.sam)
print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
             scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
             type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
             subset=fleet %in% c("Echo East"),
             main="Survey catchability parameters",ylab="Catchability",xlab="Age"))
#dev.off()

#Plot obs_variance (weightings)
#png('Fig42.png',units = "px", height=800,width=672, bg = "white")
obv <- obs.var(SARDINE.sam)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv <- obv[order(obv$value),]
bp <- barplot(obv$value,ylab="Observation Variance",
              main="Observation variances by data source",col=factor(obv$fleet))
axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
#dev.off()

#png('Fig43.png',units = "px", height=800,width=672, bg = "white")
plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
     pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
#dev.off()

#Plot fishery selectivity pattern over time
# png('Fig44.png',units = "px", height=800,width=672, bg = "white")
sel.pat <- merge(f(SARDINE.sam),fbar(SARDINE.sam),
                 by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
print(xyplot(sel ~ age|sprintf("%i's",floor((year+2)/5)*5),sel.pat,
             groups=year,type="l",as.table=TRUE,
             scale=list(alternating=FALSE),
             main="Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar"))
# dev.off()

#Plot correlation matrix
cor.plot(SARDINE.sam)

#Plot otholith
# plot.otolith(SARDINE.sam,n=1000) #Warning, this takes very long!

dev.off()

### ======================================================================================================
### Reference points
### ======================================================================================================
library(FLBRP)
ref. <- brp(FLBRP(ANCHOVY,fbar=seq(0,1,length.out=101),nyears=3))
print(refpts(ref.))

ANCHOVY.SRR <- FLSR(
  rec = rec(ANCHOVY)[,ac((range(ANCHOVY)["minyear"]+1): range(ANCHOVY)["maxyear"])],
  ssb = ssb(ANCHOVY)[,ac((range(ANCHOVY)["minyear"])  :(range(ANCHOVY)["maxyear"]-1))],
  model='segreg')
ANCHOVY.SRR <- fmle(ANCHOVY.SRR)
plot(ANCHOVY.SRR)

newData <- predict(ANCHOVY.SRR,ssb=FLQuant(seq(0,max(ssb(ANCHOVY)),length.out=200)))
yrange  <- range(pretty(c(0,range(rec(ANCHOVY)))))/1e6; xrange <- range(pretty(c(0,range(ssb(ANCHOVY)))))/1e6
plot(y=newData/1e6,x=seq(0,max(ssb(ANCHOVY)),length.out=200)/1e6,type="l",lwd=2,
     xlab="SSB (million tonnes)",ylab="Recruitment (billions)",xlim=xrange,ylim=yrange,
     las=1,cex.lab=1.3,cex.axis=1.1,xaxs="i",yaxs="i")
points(y=rec(ANCHOVY)/1e6,x=ssb(ANCHOVY)/1e6)

dev.off()

### ======================================================================================================
### Reference points
### ======================================================================================================
library(FLBRP)
ref. <- brp(FLBRP(ANCHOVY,fbar=seq(0,1,length.out=101),nyears=3))
print(refpts(ref.))

ANCHOVY.SRR <- FLSR(
  rec = rec(ANCHOVY)[,ac((range(ANCHOVY)["minyear"]+1): range(ANCHOVY)["maxyear"])],
  ssb = ssb(ANCHOVY)[,ac((range(ANCHOVY)["minyear"])  :(range(ANCHOVY)["maxyear"]-1))],
  model='segreg')
ANCHOVY.SRR <- fmle(ANCHOVY.SRR)
plot(ANCHOVY.SRR)

newData <- predict(ANCHOVY.SRR,ssb=FLQuant(seq(0,max(ssb(ANCHOVY)),length.out=200)))
yrange  <- range(pretty(c(0,range(rec(ANCHOVY)))))/1e6; xrange <- range(pretty(c(0,range(ssb(ANCHOVY)))))/1e6
plot(y=newData/1e6,x=seq(0,max(ssb(ANCHOVY)),length.out=200)/1e6,type="l",lwd=2,
     xlab="SSB (million tonnes)",ylab="Recruitment (billions)",xlim=xrange,ylim=yrange,
     las=1,cex.lab=1.3,cex.axis=1.1,xaxs="i",yaxs="i")
points(y=rec(ANCHOVY)/1e6,x=ssb(ANCHOVY)/1e6)

dev.off()

### ======================================================================================================
### Document Assessment
### ======================================================================================================
##Create a standard table
stockSummaryTable <- cbind(rec(ANCHOVY.sam)$year,
                           rec(ANCHOVY.sam)$value,      rec(ANCHOVY.sam)$lbnd,    rec(ANCHOVY.sam)$ubnd,
                           tsb(ANCHOVY.sam)$value,      tsb(ANCHOVY.sam)$lbnd,    tsb(ANCHOVY.sam)$ubnd,
                           ssb(ANCHOVY.sam)$value,      ssb(ANCHOVY.sam)$lbnd,    ssb(ANCHOVY.sam)$ubnd,
                           catch(ANCHOVY.sam)$value,    catch(ANCHOVY.sam)$lbnd,  catch(ANCHOVY.sam)$ubnd,
                           catch(ANCHOVY.sam)$value / ssb(ANCHOVY.sam)$value, catch(ANCHOVY.sam)$lbnd / ssb(ANCHOVY.sam)$lbnd, catch(ANCHOVY.sam)$ubnd / ssb(ANCHOVY.sam)$ubnd,
                           fbar(ANCHOVY.sam)$value,     fbar(ANCHOVY.sam)$lbnd,   fbar(ANCHOVY.sam)$ubnd,
                           c(quantMeans(harvest(ANCHOVY.sam)[ac(0:1),])),
                           c(sop(ANCHOVY),NA))

colnames(stockSummaryTable) <-
  c("Year",paste(rep(c("Recruits Age 0 (Thousands)","Total biomass (tonnes)","Spawing biomass (tonnes)",
                       "Landings (tonnes)","Yield / SSB (ratio)","Mean F ages 2-6"),each=3),c("Mean","Low","High")),"Mean F ages 0-1","SoP (%)")
stockSummaryTable[nrow(stockSummaryTable),] <- NA
stockSummaryTable[nrow(stockSummaryTable),"Spawing biomass (tonnes) Mean"] <- 2271364
stockSummaryTable[nrow(stockSummaryTable),2:4] <- c(rec(ANCHOVY.sam)$value[nrow(rec(ANCHOVY.sam))],rec(ANCHOVY.sam)$lbnd[nrow(rec(ANCHOVY.sam))],rec(ANCHOVY.sam)$ubnd[nrow(rec(ANCHOVY.sam))])
write.csv(stockSummaryTable,file=file.path(output.dir,"stockSummaryTable.csv"))

### ============================================================================
### Finish
### ============================================================================
dev.off()


#RETROSPECTIVE
# dir.create("c:\\samtemp\\")
retro_pil <- retro(SARDINE,SARDINE.tun,SARDINE.ctrl, retro=1)
plot(retro_ane)
retro_sar <- retro_ane
save(retro_sar, file='retro_sar.Rdata')
getwd()


#### PLOT CATCH
png('catch_PIL.png',units = "px", height=700,width=800, bg = "white")
ggplot(data = catch(SARDINE), aes(year, data)) + geom_point() + geom_line() + 
  ylab("Catch (t)") + xlab("Year") + ggtitle('PIL GSA 17-18') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

cc <- read.csv('C:/Users/s.angelini/CNR/STECF_small pelagics 2017/SARDINE/catch by country.csv', header=T, sep=';')
head(cc)
png('catch_PIL.png', units = "cm", height=12,width=16, res= 600, bg = "white")
ggplot(cc, aes(year, Catch, group=Zone, fill=Zone, col=Zone))+
  geom_line(size=0.8) + theme_bw() + ggtitle('PIL GSA 17-18') + 
  xlab('Year') + ylab('Catch (tonnes)')+
  theme(plot.title = element_text(hjust = 0.5, face='bold', size=18),
        axis.text=element_text(size = 10),
        axis.title=element_text(size=12))
dev.off()

