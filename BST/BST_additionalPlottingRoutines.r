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

manualRetro <- function(BST.sam,BST,BST.tun,BST.ctrl){
  BST.retro <- list()
  BST.retro[[ac(2018)]] <- BST.sam
  for(iRetroYr in 2017:2014){
    rt.stck <- window(BST,end=iRetroYr)
    rt.tun  <- BST.tun
    for(iTun in names(rt.tun)){
      if(range(rt.tun[[iTun]])["maxyear"] >= iRetroYr)
        rt.tun[[iTun]] <- rt.tun[[iTun]][,-which(dimnames(rt.tun[[iTun]]@index)$year %in% ((iRetroYr+1):range(rt.tun[[iTun]])["maxyear"]))]
    }
    rt.ctrl <- BST.ctrl
    rt.ctrl@range["maxyear"] <- iRetroYr
    output <- capture.output(rt.sam <- try(FLSAM(rt.stck,rt.tun,rt.ctrl)))
    if(class(rt.sam)=="try-error"){rt.sam <- NA}
    BST.retro[[ac(iRetroYr)]] <- rt.sam
  }
  return(as(BST.retro,"FLSAMs"))}
