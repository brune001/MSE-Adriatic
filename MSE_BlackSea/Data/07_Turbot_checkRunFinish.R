#check which runs finished
setwd("D:/Max files for backup/Documents/Commitees/GFCM/Turbot MSE/MSE SAM/Results/Turbot/simres/")

fls <- dir()[grep(".RData",dir())]
scens <- unlist(lapply(as.list(fls),function(x){return(strsplit(x,".RData")[[1]][1])}))
finalYear <- matrix(NA,nrow=length(fls),ncol=1,dimnames=list(scenario=scens,finalYear="finalYear"))
for(iFile in fls){
  load(iFile)
  iters <- dims(restosave$pstk)$iter
  years <- dims(restosave$TAC)$minyear:dims(restosave$TAC)$maxyear
  NAyear <- which(is.na(restosave$TAC),arr.ind=T)
  iScen <- strsplit(iFile,".RData")[[1]]
  if(length(NAyear)>0){
    finalYear[iScen,] <- min(years[an(names(table(NAyear[,2])==iters))])
  } else {
    finalYear[iScen,] <- 2036
  }
}