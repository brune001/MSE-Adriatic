



pdf(paste(runname,"plot.pdf")) 
library(ggplotFL)
#### load previous assessment for comparison
assess2 <-assess+assess.sam 

units(assess2@harvest)  <- "f"

# figure - assessment result, spawning stock biomass, fishing mortality, recruitment
print(plot(assess.sam))
print(plot(assess2))


# figure - catchabilities at age from HERAS
catch <- catchabilities(assess.sam)
print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
       scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
       type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
       subset=fleet %in% names(assess.tun),
       main="Survey catchability parameters",ylab="Catchability",xlab="Age"))

# figure - variance by data source
obv <- obs.var(assess.sam)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv <- obv[order(obv$value),]
bp <- barplot(obv$value,ylab="Observation Variance",
              main="Observation variances by data source",col=factor(obv$fleet))
axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)

## figure - variance vs uncertainty for each data source
plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
     pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
#
# figure - fishing age selectivity per year
sel.pat <- merge(f(assess.sam),fbar(assess.sam),
                 by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
xyplot(sel ~ age|sprintf("%i's",floor((year+2)/5)*5),sel.pat,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar")

# figure - correlation matrix of model parameters
print(cor.plot(assess.sam))


if(assess.ctrl@residuals  ==T)
{
## figure - catch residuals per year per age
dat <- subset(residuals(assess.sam),fleet=="catch unique" & year >2001)
print(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year Catch",
       panel=function(...){
        lst <- list(...)
         panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=3*abs(lst$cex))
       }))
#
## figure - acosutic index residuals per year per age
for( su in 1:length(assess.tun))
{
dat <- subset(residuals(assess.sam),fleet==names(assess.tun)[su])
print(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main=paste("Residuals by year",names(assess.tun)[su]),
       panel=function(...){
         lst <- list(...)
         panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=3*abs(lst$cex))
       }))
}       
#       

# residual.diagnostics(assess.sam)
}
#library(stockassessment)
for (fl in 1:length(assess.ctrl@fleets))   print(fitplot(fit , fleet = fl)       )
#
detach()




## process error in terms of N
#print(procerr.plot(assess+assess.sam,weight="stock.wt",type="n",rel=T))
#
## process error in terms of additional mortality
#print(procerr.plot(assess+assess.sam,weight="stock.wt",type="mort",rel=F))
#

#print(plot(assess.retro,futureYrs=F))
#
## model parameters retrospective
#retroParams(assess.retro)
#
## stock trajectory retrospective
#retroSelectivity(assess.retro,2009:range(assess)["maxyear"])
#
#print(plot(assess.loo))
#
dev.off()