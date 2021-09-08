#==============================================================================
# libraries
#==============================================================================
rm(list=ls())
library(FLash)
#library(FLasher)
library(FLAssess)
library(ggplotFL)
library(FLBRP)
library(FLSAM)
library(FLCore)
library(FLEDA)
library(doBy)
library(reshape)
library(devtools)
#install_github("ices-tools-prod/msy")
library(msy)


#==============================================================================
# running repertory and paths
#==============================================================================

species <- "ANCHOVY"  # "ANCHOVY" or "SARDINE"
assess.name <- "Anchovy GSA 17-18_tbmSAM"

#species <- "SARDINE"  # "ANCHOVY" or "SARDINE"
#assess.name <- "Sardine GSA 17-18_tbmSAM"
#

# path to local github repository
setwd("C:/Users/brune001/my git files/MSE-Adriatic/")

# source needed functions
run <- "full"

if(run == "full")  fname <-  paste0("./Results/",species,"/",assess.name,"_250iters_20yrs_blank_objects_MSE.RData")
if(run == "short") fname <-  paste0("./Results/",species,"/",assess.name,"_2iters_12yrs_blank_objects_MSE.RData")

load(fname)

source('./Code/by_modulesTMB/MSE_functions.R')
source('./Code/by_modulesTMB/03_ADRIAMED_MSE_BRPs_and_Scenarios.r')



sc <- names(management.scenarios)[c(1:11,24:30)]



#==============================================================================
# time periods
#==============================================================================

ST <- ac(2017:2021)
MT <- ac(2022:2027)
LT <- ac(2028:2035)


#==============================================================================
# extract output for each scenario
#==============================================================================

results <- lapply(sc , function(x) 
            {
            cat(x,"\n")
            # load data  and rename / reshapre  and compute what's needed
            load(file = paste0("./Results/",species,"/simres/",x,"_",it,"its_",fy,"V2.RData"))
            res<-restosave
            pstk <- window(res$pstk,end = range(res$pstk)["maxyear"] - 1)
            Fad <- window(res$Fad,end = range(res$pstk)["maxyear"] - 1)
            SSBad <-window(res$SSBad,end = range(res$pstk)["maxyear"] - 1)
            TAC <- window(res$TAC,end = range(res$pstk)["maxyear"] - 1)
            #assessment errors
            devSSB <-  (SSBad -  ssb(pstk)[,dimnames(SSBad)$year])  / ssb(pstk)[,dimnames(SSBad)$year]
            devF   <-  (Fad   - fbar(pstk)[,dimnames(Fad)$year])   / fbar(pstk)[,dimnames(Fad)$year]            
            #risk
            rsk <- iterSums(ssb(pstk) <= blim) / it
            
            

            # plot the output
                # stock trends
            png(paste0("Results/",species,"/PLOTS/",x,run,"_OM trends_V2.png"), width=6, height=10,units = "in" , res = 300)
            plot.iStk(pstk ,nits = 2 , title = paste("scenario",x))
            dev.off()
                # assessment errors
            png(paste0("Results/",species,"/PLOTS/",x,run,"_SSB error_V2.png"), width=6, height=4,units = "in" , res = 300)
            plot.iQuant(devSSB ,nits = 2 , title = "% error SSB Advice Year")
            dev.off()

            png(paste0("Results/",species,"/PLOTS/",x,run,"_Fbar error_V2.png"), width=6, height=4,units = "in" , res = 300)
            plot.iQuant(devF ,nits = 2 , title = "% error F Advice Year")
            dev.off()

                # risk

             time <- data.frame(time = c("Short Term","Medium Term","Long Term") , x = c(mean(an(ST)),mean(an(MT)),mean(an(LT))),y=0.9)

             P<- plot(rsk) 
             P<- P + geom_rect(data = as.data.frame(iter(rsk,1)),aes(xmin = min(an(ST)) ,xmax = max(an(ST)) , ymin =0 , ymax = 1 ) , fill = "blue" , alpha = 0.02) 
             P<- P + geom_rect(data = as.data.frame(iter(rsk,1)),aes(xmin = min(an(MT)) ,xmax = max(an(MT)) , ymin =0 , ymax = 1 ) , fill = "green" , alpha = 0.01)
             P<- P + geom_rect(data = as.data.frame(iter(rsk,1)),aes(xmin = min(an(LT)) ,xmax = max(an(LT)) , ymin =0 , ymax = 1 ) , fill = "red" , alpha = 0.02) 
             P<- P + xlim(ay,fy)+ ylim(0,1) +ylab("p(SSB<Blim)") + ggtitle("probability of falling below Blim")
             P<- P + geom_text(data=time , aes(x=x , y=y , label = time)) 
             png(paste0("Results/",species,"/PLOTS/",x,run,"_RiskBlim_V2.png"), width=6, height=4,units = "in" , res = 300)
             print(P)
             dev.off()

             # table of diagnostics
             ssb<- lapply (list(ST,MT,LT) , FUN = function(x) try(c(apply(yearMeans(ssb(pstk)[,x]),c(1:5),median)),silent=T))
             names(ssb) <- c("Short Term","Medium Term","Long Term")
             catch<- lapply (list(ST,MT,LT) , FUN = function(x) try(c(apply(yearMeans(catch(pstk)[,x]),c(1:5),median)),silent=T))
             names(catch) <- c("Short Term","Medium Term","Long Term")
             Fbar<- lapply (list(ST,MT,LT) , FUN = function(x) try(c(apply(yearMeans(fbar(pstk)[,x]),c(1:5),median)),silent=T))
             names(Fbar) <- c("Short Term","Medium Term","Long Term")
             risk<- lapply (list(ST,MT,LT) , FUN = function(x) try(c(apply(rsk[,x],c(1,3:6),FUN=max)),silent=T))
             names(risk) <- c("Short Term","Medium Term","Long Term")
           
             ssb<- unlist(ssb)
             ssb <- data.frame(scenario = x , time = names(ssb) , var = "ssb" , val = ssb)
             catch<- unlist(catch)
             catch <- data.frame(scenario = x , time = names(catch) , var = "catch" , val = catch)
             Fbar<- unlist(Fbar)
             Fbar <- data.frame(scenario = x , time = names(Fbar) , var = "Fbar" , val = Fbar)
             risk<- unlist(risk)
             risk <- data.frame(scenario = x , time = names(risk) , var = "risk" , val = risk)
            
             diags <- do.call (rbind , list(ssb,catch,Fbar,risk))
             
             # table of raw output for post hoc socio economic study
             ssb <- as.data.frame(ssb(pstk))
             ssb$variable <- "ssb"
             catch <- as.data.frame(catch(pstk))
             catch$variable <- "catch"
             fbar <- as.data.frame(fbar(pstk))
             fbar$variable <- "fbar"
             raw<-do.call(rbind, list(ssb,catch,fbar))
             raw$scenario <- x
             raw<-raw[,c("year","scenario", "variable" , "iter","data")]
             names(raw)[5] <- "value"
             
             output<-list(diagnostics = diags , detailed.output = raw)
             
             return(output)
             
             } )
             
             
diagnostics     <- do.call(rbind.data.frame, lapply(results , function(x) x[[1]]))             
detailed.output <- do.call(rbind.data.frame, lapply(results , function(x) x[[2]])) 

write.csv(diagnostics,file = paste0("Results/",species,"/diagnostics_",run,"_MSE_V2.csv"), row.names=FALSE)
write.csv(detailed.output,file = paste0("Results/",species,"/detailed.output_",run,"_MSE_V2.csv"), row.names=FALSE)




