# turbo example for luminophore pulse experiments
# author: Clemens Abraham

################################################################################
### Housekeeping ###
################################################################################

# library
library(deSolve)
library(rootSolve)
library(coda)
library(FME)
library(tidyverse)
library(Matrix)
library(plotrix)

# load customized functions
load(file = "turbo pulse functions.RData")

################################################################################
### Example ###
################################################################################

# data
data <- read.table("luminophore.txt", header = TRUE)
data <- subset(data, profile == "a")

# global parameters
parms <- list(
  depth = max(data$end),
  db = 0.01,
  days = 14,
  dx = 0.5,
  cakethickness = 0.5,
  sl = 2,
  wt =20,
  k =0,
  flux = 0,
  fluxintroduction = 0
)

parms$slicenumber <- parms$depth/parms$dx
parms$times <- c(0,parms$days)

# data profile
datalimits <- unique(c(data$start,data$end))
datamidpoints <- midpoints(datalimits)
conc_profile <- numtoconc(data$lum,datalimits)
dataprofileframe <- data.frame(depth=datamidpoints, concentration=conc_profile)
initprof <- initialprofile(conc_profile,datalimits,parms$slicenumber,
                           parms$dx,parms$cakethickness)
initlimits <- seq(0,parms$slicenumber*parms$dx,by=parms$dx)
initmidpoints <- midpoints(initlimits)
initprofileframe <- data.frame(depth=initmidpoints, concentration=initprof)

### diffusive model ###

# fit diffusive model (ODE)
fit <- modFit(p=parms$db,f=diffobjective,lower=c(0))                            # minimizes cost of "diffobjective", p = init. value
                                                                                # lower = lower bounds on db
fitdb <- as.numeric(fit$par[1])

diff_modelprofileframe <- modelprofileframe_diff(fitdb)

### non local model ###

# fit non local model (CTRW)
fitnl <- modFit(p = c(sl = parms$sl, wt = parms$wt), f = nlobjective_analytic,  # minimizes cost of "nlobjective_analytic", p = init. values
                lower = c(0,0))                                                 # lower = lower bounds on sl and wt

fitsl <- as.numeric(fitnl$par[1])                                               # fitted step length (sl)
fitwt <- as.numeric(fitnl$par[2])                                               # fitted waiting time (wt)
fitnldb <- fitsl^2/(2*fitwt)                                                    # calculates db (sl^2/2*wt)

nl_modelprofileframe <- modelprofileframe_nl(sl = fitsl, wt = fitwt)

### results ###

# plot the results
par(lwd = 1, cex = 1, cex.lab = 1, cex.axis = 1, mar = c(5, 5, 2, 2))
revaxis(dataprofileframe$concentration, dataprofileframe$depth,
        xlab = "concentration of luminophores", ylab = "depth in cm", pch = 4)

lines(initprofileframe$concentration,      -initprofileframe$depth,      lty = 1)
lines(diff_modelprofileframe$concentration, -diff_modelprofileframe$depth, lty = 2)
lines(nl_modelprofileframe$concentration,  -nl_modelprofileframe$depth,   lty = 3)

legend("right",
       legend = c("initial profile", "observed data",
                  paste0("diffusive (Db = ", round(fitdb, 3), ")"),
                  paste0("non-local (Db = ", round(fitnldb, 3), ")")),
       lty = c(1, NA, 2, 3),
       pch = c(NA, 4, NA, NA),
       bty = "n")
