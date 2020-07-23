# function to do the jit runs
jitANE <- function (fit, nojit = 10, par = defpar(fit$data, fit$conf),
    sd = 0.25)
{
    parv <- unlist(par)
    pars <- lapply(1:nojit, function(i) relist(parv + rnorm(length(parv),
        sd = sd), par))
     fits <- lapply(pars, function(p) sam.fit(fit$data, fit$conf, p, silent = TRUE))

    attr(fits, "fit") <- fit
    attr(fits, "jitflag") <- 1
    class(fits) <- c("samset")
    res <- list(fits,pars)
    return(res)
}







# check  that the assessment can be reproduced using the stockassessment library
dats <-ANE10b$data
confs<-ANE10b$conf
ini.pars <- ANE10$pl

#reproduce ANE10b
fit<-sam.fit(dats,confs, defpar(dats,confs))      # works with default initial pars
fit<-sam.fit(dats,confs, ini.pars)               # doesnt work , says init pars have ot the wrong shape

# take the  default init par from ANE10b, put the values of ANE10 manually
init.pars2 <- defpar(dats,confs)

# fill in manually
init.pars2$logFpar  <-             ini.pars$logFpar
init.pars2$logSdLogFsta <-         ini.pars$logSdLogFsta
init.pars2$logSdLogN    <-         ini.pars$logSdLogN
init.pars2$logSdLogObs  <-         ini.pars$logSdLogObs
init.pars2$itrans_rho   <-         ini.pars$itrans_rho
#
fit<-sam.fit(dats,confs, init.pars2)
fit
ANE10b   # yes! it works   same likelihood





# do the jittering
njit <- 50

jit<- jitANE(ANE10b,par = init.pars2,nojit=njit,sd=1.5)


ANE10b.jit <- jit[[1]]  # the assessments
pars <- jit[[2]]        # the set of initial parameters

pars2 <- lapply(pars,function(x) {           # reformat the initial parameters
                                  x<-unlist(x)
                                  x<-x[1:dim(partable(fit))[1]]
                                  return(x)})
                                  
pars3<- data.frame(jit=sort(rep(1:dim(partable(fit))[1] ,njit)),par = names(unlist(pars2)),val =unlist(pars2))


                                  
library(ggplot2)
ggplot(pars3 , aes(x=par,y=val)) + geom_boxplot()              +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))   ### all the starting values tested

plot(ANE10b.jit)           ### all the runs are the same


lapply(ANE10b.jit,function(x) return(x))   # same lofg likelihood == same runs