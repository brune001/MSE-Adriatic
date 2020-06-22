
library(FLCore)
library(lattice)
library(scales)
library(gplots)




load("./data/ANEstk.RData")


pdf(file.path("./results/","input data base.pdf"))

ca<-catch.n(ANE)


library(ggplotFL)
print(plot(ca) + ggtitle("catch N@age"))
detach()


years <- an(dimnames(ca)$year)
ages <-  an(dimnames(ca)$age)
ltCol   <- length(ages)       # number of colours
cols <- rich.colors(ltCol)



CatARaw  <-(matrix(ca@.Data,ncol = length(years)))




plot(NA,NA,main="ANCHOVY: Catch at Age (N; Observed)", xlab="Year", ylab="Age",xlim=c(min(years),max(years)), ylim=c(0,max(ages)), yaxt="n")
axis(2,at=seq(0,12,by=1),labels=ac(0:12))
for (a in rev(ages)) {
  radius <- as.numeric((abs(CatARaw[a+1,]) / pi))
  points(years, rep(a,ncol(CatARaw)),cex=radius/100000, pch=21, col=1, bg=alpha(cols[a+1],0.8))
}
par(xpd=F)
for (cstrt in  c((min(years)-length(ages)):(max(years)+length(ages)))) lines(c((cstrt):(cstrt+length(ages)-1)),c(0:range(ANE)["max"]),lty=3)




ccplot(data~year, data=log(FLCohort(catch.n(ANE))), type="l",col=1)





sp<-FLIndex(index=catch.n(ANE))
name(sp)<-"f"
plot(sp,type="internal",main = "catch at age")




load("./data/ANEtun.RData")

for (su in 1:2)
{

surv <- name(ANE.tun[[su]])
plot(ANE.tun[[su]],type="internal",main=surv)

if (su == 1)  plot(trim(ANE.tun[[su]],year=2004:2014),type="internal",main=paste(surv,"pre2015"))
if (su == 1)  plot(trim(ANE.tun[[su]],year=2015:2018),type="internal",main=paste(surv,"post2015"))

library(ggplotFL)
print(plot(ANE.tun[[su]]@index) + ggtitle(paste(surv, ": N@age")))
detach()
}



biomsurv <- data.frame(as.data.frame(ANE.tun[[3]]@index), survey = name(ANE.tun[[3]]) )

library(ggplot2)

ggplot( biomsurv  , aes(year , data , colour = survey ) ) + ggtitle("biomass indices") + geom_line()  + ylab("biomass index")

dev.off()


