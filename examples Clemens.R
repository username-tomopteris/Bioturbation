# package turbo examples (git hub version 2025)
# functions for fitting model profiles to tracer data
# Clemens Abraham

### Housekeeping ####
library(FME)
library(deSolve)
library(rootSolve)
library(coda)
library(tidyverse)
library(plotrix)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
graphics.off()

# load functions from "turbo"
load(file = "turbo functions.RData")

#### parameters, data, initialization ####

# parameters

slicenumber <- 200 # number of modeled slices
dx <- 0.1 # slice thickness <= cakethickness
cakethickness <- 0.5 # luminophore cake thickness >= dx
days <- 14 # incubation/ observation time

# read data (test data set "luminophore")

datafile <- "C:/Users/cleme/Desktop/Uni/HIWI Job/BenTools Maps/turbo R GitHub/pck turbo/data/luminophore.txt"
data <- read.table(datafile, header=TRUE)

# profile names

profilenames <- as.vector(unique(data$profile))

# loop through profiles

dblist <- NULL
nldblist <- NULL
sllist <- NULL
wtlist <- NULL

### automated modeling (png) ####
for(profilename in profilenames){
  
  cat(profilename,"\n")
  
  dataframe <- subset(data,profile==profilename)
  datalimits <- unique(c(dataframe$start,dataframe$end))
  datamidpoints <- midpoints(datalimits)
  dataprofile <- numtoconc(dataframe$lum,datalimits)
  dataprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=dataprofile))
  
  # create initial profile
  
  initprofile <- initialprofile(dataprofile,datalimits,slicenumber,dx,cakethickness)
  initlimits <- seq(0,slicenumber*dx,by=dx)
  initmidpoints <- midpoints(initlimits)
  
  # fit nonlocal model
  
  fit <- modFit(p=c(2,20),f=nlobjective,lower=c(0,0))
  fitsl <- fit$par[1]
  fitwt <- fit$par[2]
  fitnldb <- fitsl^2/(2*fitwt)
  
  times <- c(0,days)
  final <- nonlocal(times,initprofile,list(waitingtime=fitwt,steplength=fitsl,dx=dx))
  
  # plot data and nonlocal model
  
  png(file=paste(profilename,".png",sep=""))
  plot(dataprofile, datamidpoints, type="n", ylab="depth",
       xlab="luminophore concentration", ylim = c(max(datamidpoints), 0))
  lines(final[2,2:(slicenumber+1)], initmidpoints,col="red",lwd=2)
  
  # fit diffusive model
  
  fit <- modFit(p=c(0.01),f=diffobjective,lower=c(0))
  fitdb <- fit$par[1]
  
  times <- c(0,days)
  final <- diffusive(times,initprofile,list(db=fitdb,dx=dx))
  lines(final[2,2:(slicenumber+1)], initmidpoints, col="green",lwd=2)
  
  # plot data
  
  points(dataprofile, datamidpoints, pch=21, bg="white")
  
  # legend
  
  legend(10,max(dataprofile)/2,pch=c(NA,NA,21),lty=c("solid","solid","blank"),
         lwd=c(2,2,1),col=c("red","green","black"),
         legend=c("nonlocal model","diffusive model","data"),bty="n")
  
  # results
  
  # cat("\nfitted step length: ",round(fitsl,digits=3),"\nfitted waiting time: ",round(fitwt,digits=3),"\nfitted db: ",round(fitdb,digits=3),"\nfitted nldb: ",round(fitnldb,digits=3),"\n\n",sep="")
  
  dblist <- c(dblist,round(fitdb,digits=4))
  nldblist <- c(nldblist,round(fitnldb,digits=4))
  sllist <- c(sllist,round(fitsl,digits=4))
  wtlist <- c(wtlist,round(fitwt,digits=4))    
  
  # mtext(dblist, side = 3, line = 0)
  
  # close graphics device
  
  dev.off()
  
}

### automated modeling (R plot) ####
for(profilename in profilenames){
  
  cat(profilename,"\n")
  
  dataframe <- subset(data,profile==profilename)
  datalimits <- unique(c(dataframe$start,dataframe$end))
  datamidpoints <- midpoints(datalimits)
  dataprofile <- numtoconc(dataframe$lum,datalimits)
  dataprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=dataprofile))
  
  # create initial profile
  
  initprofile <- initialprofile(dataprofile,datalimits,slicenumber,dx,cakethickness)
  initlimits <- seq(0,slicenumber*dx,by=dx)
  initmidpoints <- midpoints(initlimits)
  
  # fit nonlocal model
  
  fit <- modFit(p=c(2,20),f=nlobjective,lower=c(0,0))
  fitsl <- fit$par[1]
  fitwt <- fit$par[2]
  fitnldb <- fitsl^2/(2*fitwt)
  
  times <- c(0,days)
  final_nl <- nonlocal(times,initprofile,list(waitingtime=fitwt,steplength=fitsl,dx=dx))
  
  # fit diffusive model
  
  fit <- modFit(p=c(0.01),f=diffobjective,lower=c(0))
  fitdb <- fit$par[1]
  
  times <- c(0,days)
  final_diff <- diffusive(times,initprofile,list(db=fitdb,dx=dx))
  
  # plot data and models
  
  plot(dataprofile, datamidpoints, type="n", ylab="depth",
       xlab="luminophore concentration", main = profilename, ylim = c(max(datamidpoints), 0),
       xlim = c(0, sort(dataprofile, decreasing = TRUE)[2]))
    
    # nonlocal
    final_nl <- final_nl[2,2:(slicenumber+1)]
    lines(final_nl, initmidpoints,col="orange2",lwd=2, lty = 2)
    
    # diffusive
    final_diff <- final_diff[2,2:(slicenumber+1)]
    lines(final_diff, initmidpoints, col="palegreen3",lwd=2)
    
    # data points
    
    points(dataprofile, datamidpoints, pch=21, bg="white")
    
    legend("right",max(dataprofile)/2,pch=c(NA,NA,21),lty=c(1, 2, 0),
    lwd=c(2,2,1),col=c("palegreen3","orange2","black"),
    legend=c("diffusive model","nonlocal model","data"),bty="n")
    
    # R square values (sum of squared errors)
    
    R2_diff <- 1 - sum((dataprofile - final_diff[seq(1, length(final_diff),
                by = length(final_diff) / length(datamidpoints))])^2) /
                sum((dataprofile - mean(dataprofile))^2)
    
    R2_nl <- 1 - sum((dataprofile - final_nl[seq(1, length(final_nl),
                by = length(final_nl) / length(datamidpoints))])^2) /
                sum((dataprofile - mean(dataprofile))^2)
    
    # results
    
    legend("bottomright", 
           legend = c(paste("Db (diff):", round(fitdb,digits=4)),
                      paste("Db (nl):", round(fitnldb,digits=4)),
                      paste("waiting time (nl):", round(fitwt,digits=4)),
                      paste("step length (nl):", round(fitsl,digits=4)),
                      paste("R2 (diff):", round(R2_diff, digits = 4)),
                      paste("R2 (nl):", round(R2_nl, digits = 4))),
           bty = "n") # no frame
  
  cat("\nfitted step length: ",round(fitsl,digits=3),"\nfitted waiting time: ",round(fitwt,digits=3),"\nfitted db: ",round(fitdb,digits=3),"\nfitted nldb: ",round(fitnldb,digits=3),"\n\n",sep="")
  
  dblist <- c(dblist,round(fitdb,digits=4))
  nldblist <- c(nldblist,round(fitnldb,digits=4))
  sllist <- c(sllist,round(fitsl,digits=4))
  wtlist <- c(wtlist,round(fitwt,digits=4))
}

### R square values ####


### all results ####

allresults <- data.frame(profilenames,dblist,nldblist,sllist,wtlist,stringsAsFactors=FALSE)
names(allresults) <- c("profile","db","nldb","steplength","waitingtime")

test_results <- write.csv(allresults,file="results.csv")