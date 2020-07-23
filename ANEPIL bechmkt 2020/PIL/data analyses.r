
library(FLCore)
library(lattice)
library(scales)
library(gplots)




load("./data/Sarstk.RData")


pdf(file.path("./Results/","input data new.pdf"))

ca<-catch.n(Sar)


library(ggplotFL)
print(plot(ca) + ggtitle("catch N@age"))
detach()


years <- an(dimnames(ca)$year)
ages <-  an(dimnames(ca)$age)
ltCol   <- length(ages)       # number of colours
cols <- rich.colors(ltCol)



CatARaw  <-(matrix(ca@.Data,ncol = length(years)))




plot(NA,NA,main="SARDINE: Catch at Age (N; Observed)", xlab="Year", ylab="Age",xlim=c(min(years),max(years)), ylim=c(0,max(ages)), yaxt="n")
axis(2,at=seq(0,12,by=1),labels=ac(0:12))
for (a in rev(ages)) {
  radius <- as.numeric((abs(CatARaw[a+1,]) / pi))
  points(years, rep(a,ncol(CatARaw)),cex=radius/100000, pch=21, col=1, bg=alpha(cols[a+1],0.8))
}
par(xpd=F)
for (cstrt in  c((min(years)-length(ages)):(max(years)+length(ages)))) lines(c((cstrt):(cstrt+length(ages)-1)),c(0:range(Sar)["max"]),lty=3)




ccplot(data~year, data=log(FLCohort(catch.n(Sar))), type="l",col=1)


ccplot(data~year, data=log(FLCohort(index(Sar.tun[[1]]))), type="l",col=1,main=name(Sar.tun[[1]]))
ccplot(data~year, data=log(FLCohort(index(Sar.tun[[2]]))), type="l",col=1,main=name(Sar.tun[[2]]))
ccplot(data~year, data=log(FLCohort(index(Sar.tun[[3]]))), type="l",col=1,main=name(Sar.tun[[3]]))
ccplot(data~year, data=log(FLCohort(index(Sar.tun[[4]]))), type="l",col=1,main=name(Sar.tun[[4]]))

sp<-FLIndex(index=catch.n(Sar))
name(sp)<-"f"
plot(sp,type="internal",main = "catch at age")




load("./data/PILtun.RData")

for (su in 1:4)
{

surv <- name(Sar.tun[[su]])
plot(Sar.tun[[su]],type="internal",main=surv)

#if (su == 1)  plot(trim(Sar.tun[[su]],year=2004:2014),type="internal",main=paste(surv,"pre2015"))
#if (su == 1)  plot(trim(Sar.tun[[su]],year=2015:2018),type="internal",main=paste(surv,"post2015"))
#
library(ggplotFL)
print(plot(Sar.tun[[su]]@index) + ggtitle(paste(surv, ": N@age")))
detach()
}



biomsurv <- data.frame(as.data.frame(Sar.tun[[5]]@index), survey = name(Sar.tun[[5]]) )

library(ggplot2)

ggplot( biomsurv  , aes(year , data , colour = survey ) ) + ggtitle("biomass indices") + geom_line()  + ylab("biomass index")

dev.off()


