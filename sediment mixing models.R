# Implementation of different bioturbation/ transport models DRAFT
# C. Abraham/ Claude.AI

### Housekeeping ####
library(rootSolve)
library(deSolve)
library(shape)
library(ReacTran)
library(ggplot2)
library(gridExtra)

# graphics
par(cex = 1.5)
par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
par(mfrow = c(1, 3))


### Diffusive sediment mixing model (diffusion only!) ####

### analytical solution

diffusive <- function(C0 = 100, Db = 0.1, t = 30, depth)
{
  conc <- vector(length = length(depth))
  conc <- (C0/(2*sqrt(pi*Db*t)))*exp(-depth^2/(4*Db*t))
  return(conc)
}

depth <- seq(0, 15, by = 0.01)
lum <- diffusive(depth = depth)


plot(lum, depth, ylim = c(15, 0), xlab = "Luminophores",
     ylab = "cm", type = "l", lwd = 2, main = "Diffusion after 30d")

### differential equation
# graphics
par(mfrow = c(1, 2), cex = 1.5)

# 1D Diffusionsgleichung: ∂c/∂t = D * ∂²c/∂x²

diffusionsgleichung <- function(t, c, parms) {
  # Parameter
  D <- parms$D          # Diffusionskoeffizient
  dx <- parms$dx        # Räumliche Schrittweite
  n <- length(c)
  
  # Zweite räumliche Ableitung (zentrale Differenzen)
  d2cdx2 <- numeric(n)
  d2cdx2[2:(n-1)] <- (c[3:n] - 2*c[2:(n-1)] + c[1:(n-2)]) / dx^2
  
  # Randbedingungen (keine Diffusion über die Ränder)
  d2cdx2[1] <- (c[2] - c[1]) / parms$dx^2 # Oberer Rand: kein Fluss nach außen
  d2cdx2[n] <- (c[n-1] - c[n]) / parms$dx^2  # Unterer Rand: kein Fluss nach außen
  
  # Zeitableitung
  dcdt <- D * d2cdx2
  
  return(list(dcdt))
}

# Parameter definieren
x <- 0:15                   # Raumpunkte von x0 bis x15
nx <- length(x)
dx <- x[2] - x[1]
D <- 0.1                    # Diffusionskoeffizient

# Anfangsbedingung: Konzentration c0=100 bei x0, sonst 0
c0 <- numeric(nx)
c0[1] <- 100                # Bei x0 (erster Index)

# Zeitpunkte für die Lösung
zeiten <- c(0, 30, 100)

# Parameter-Liste
parms <- list(D = D, dx = dx)

# PDE lösen
loesung <- ode(y = c0, times = zeiten, func = diffusionsgleichung, 
               parms = parms, method = "bdf")
result_diffusion <- data.frame(c_t0 = loesung[1, -1], c_t30 = loesung[2, -1], c_t100 = loesung[3, -1], Tiefe = x)
plot(result_diffusion$c_t0, result_diffusion$Tiefe, ylim = c(15,0), xlim = c(0, 50), type = "l",
     main = "Diffusion", lwd = 2, xlab = "", ylab = "Tiefe in cm")
lines(result_diffusion$c_t30, result_diffusion$Tiefe, lwd = 2)
lines(result_diffusion$c_t100, result_diffusion$Tiefe, lwd = 2)


### Advection diffusion equation ####

diffusion_advektion <- function(t, c, parms) {
  n <- length(c)
  dx <- parms$dx
  Db <- parms$Db
  u <- parms$u
  
  # Diffusion: ∂²c/∂x²
  d2c <- numeric(n)
  d2c[2:(n-1)] <- (c[3:n] - 2*c[2:(n-1)] + c[1:(n-2)]) / dx^2
  
  # RANDBEDINGUNGEN für Diffusion:
  # Oberer Rand (x=0): Kein diffusiver Fluss nach oben
  d2c[1] <- -c[1] / dx^2  # Neumann BC: ∂c/∂x = 0 bei x=0
  
  # Unterer Rand (x=15): Kein diffusiver Fluss nach unten
  d2c[n] <- -c[n] / dx^2  # Neumann BC: ∂c/∂x = 0 bei x=15
  
  # Advektion: ∂c/∂x (upwind scheme für u > 0, Transport nach unten)
  dc <- numeric(n)
  dc[1] <- c[1] / dx  # Oberer Rand: Material wird advektiv eingegraben
  dc[2:n] <- (c[2:n] - c[1:(n-1)]) / dx
  
  # RANDBEDINGUNG für Advektion am unteren Rand:
  # Material verlässt das System am unteren Rand durch Advektion
  # Dies ist physikalisch korrekt: alter Tracer wird nach unten transportiert
  
  # ∂c/∂t = Db·∂²c/∂x² - u·∂c/∂x
  dcdt <- Db * d2c - u * dc
  
  list(dcdt)
}

x <- 0:15
dx <- 1
c0 <- c(100, rep(0, length(x) - 1))
parms <- list(Db = 0.1, u = 0.001, dx = dx)

# Lösung über Zeit
times <- seq(0, 100, by = 1)
loesung <- ode(c0, times, diffusion_advektion, parms, method = "bdf")

# Ergebnisse
result_advection_diffusion <- data.frame(
  c_t0 = loesung[1, -1], 
  c_t30 = loesung[31, -1], 
  c_t100 = loesung[101, -1], 
  Tiefe = x
)

loesung <- ode(c0, c(0, 30, 100), diffusion_advektion, parms, method = "bdf")
result_advection_diffusion <- data.frame(c_t0 = loesung[1, -1], c_t30 = loesung[2, -1], c_t100 = loesung[3, -1], Tiefe = x)
plot(result_advection_diffusion$c_t0, result_advection_diffusion$Tiefe, ylim = c(15, 0), xlim = c(0, 50), type = "l",
     main = "Diffusion-Advektion", lwd = 2, xlab = "", ylab = "Tiefe in cm")
lines(result_advection_diffusion$c_t30, result_advection_diffusion$Tiefe, lwd = 2)
lines(result_advection_diffusion$c_t100, result_advection_diffusion$Tiefe, lwd = 2)

### Non-local exchange sediment model ####
# example from Soetaert 2009

par(cex = 1.5)

nonlocal <- function(
    k=0.031, u=0.001, Db=0.1, depo=0.1,
    L=5, injectflux=0.3, sed)
{
  a <- (u/Db - sqrt ( (u/Db)^2 + 4*k / Db)) /2
  b <- (u/Db + sqrt ( (u/Db)^2 + 4*k / Db)) /2
  expaL <- exp(a*L)
  expbL <- exp(b*L)
  
  A <- matrix(nrow = 3, ncol = 3, byrow = TRUE, data = c(
    
    u-Db*a, u-Db*b, 0,
    expaL, expbL, -expaL,
    (u-Db*a)*expaL, (u-Db*b)*expbL, (Db*a-u)*expaL)
  )
  B <- c(depo, 0, -injectflux)
  X <- solve(A, B)
  
  s1 <- which(sed<L)
  s2 <- which(sed>=L)
  
  conc <- vector(length = length(sed))
  conc[s1] <- X[1]*exp(a*sed[s1])+X[2]*exp(b*sed[s1])
  conc[s2] <- X[3]*exp(a*sed[s2])
  return(conc)
}

depth <- seq(0, 15, by = 0.01) # cm
Pb <- nonlocal(sed=depth)

plot(Pb, depth, ylim = c(15, 0), xlab = "Pb, dpm/cm2/yr",
     ylab = "Tiefe in cm", type = "l", lwd = 2)

### Biodiffusion with regard to porosity gradients (Mulsow and Bourdreau 1998) ####

c0 <- c(100, rep(0, 30)) # initial concentration

Db <- 0.5


## empirical solid voulume fraction gradient equation (Mulsow and Bourdreau 1998) ##
solid_fraction_gradient <- function(x) {
  Phi_sInf-(Phi_sInf-Phi_s0)*exp(-alpha*x)
}

# deep sea sediment
depth <- seq(0, 15, by = 0.01)

Phi_s0 <- c(0.1)
Phi_sInf <- c(0.2)
alpha <- c(0.2)

results_solid_fraction_mud <- solid_fraction_gradient(depth)
plot(results_solid_fraction_mud, depth, ylim = c(15, 0), xlim = c(0, 1),
     type = "l", xlab = "Sedimentfraktion", ylab = "Tiefe in cm", lwd = 2)

# sandy sediment
Phi_s0 <- c(0.5)
Phi_sInf <- c(0.9)
alpha <- c(0.4)

results_solid_fraction_sand <- solid_fraction_gradient(depth)
lines(results_solid_fraction_sand, depth, lwd = 2, lty = 2)
legend("bottom", c("Tiefseeton", "Sand"), lty = c(1, 2), lwd = 2)

## constant porosity mixing: ∂c/∂t = Phi_s*Db*∂²c/∂x² ##

Phi_s <- mean(results_solid_fraction_mud) # deep sea sediment

diffusion <- function(t, c, parms) {
  n <- length(c)
  d2c <- numeric(n)
  d2c[2:(n-1)] <- (c[3:n] - 2*c[2:(n-1)] + c[1:(n-2)]) / parms$dx^2
  d2c[1] <- (c[2] - c[1]) / parms$dx^2 # Oberer Rand: kein Fluss nach außen
  d2c[n] <- (c[n-1] - c[n]) / parms$dx^2 # Unterer Rand: kein Fluss nach außen
  list(parms$D_eff * d2c)
}

x <- seq(0, 15, by = 0.5) # depth 0-15cm by 0.5cm
dx <- x[2] - x[1]
parms <- list(D_eff = Phi_s * Db, dx = dx)

loesung <- ode(y = c0, times = c(0, 30), func = diffusion, parms = parms, method = "bdf")
result_constant_porosity <- data.frame(c_t30 = loesung[nrow(loesung), -1], Tiefe = x)
plot(result_constant_porosity, type = "l", main = "Constant porosity mixing", ylim = c(15, 0))


## intraphase mixing: ∂Phi_s*c/∂t = Db * ∂/∂x(Phi_s*∂c/∂x) ##

x <- seq(0, 15, by = 0.5) # depth 0-15cm by 0.5cm
results_solid_fraction_mud <- solid_fraction_gradient(x)
Phi_s <- results_solid_fraction_mud

diffusion <- function(t, c, parms) {
  n <- length(c)
  Phi_s <- parms$Phi_s
  Db <- parms$Db
  dx <- parms$dx
  
  # ∂c/∂x an den Zellgrenzen (i+1/2)
  dcdx <- numeric(n + 1)
  dcdx[2:n] <- (c[2:n] - c[1:(n-1)]) / dx
  dcdx[1] <- 0  # Oberer Rand
  dcdx[n+1] <- 0  # Unterer Rand
  
  Phi_s_grenzen <- numeric(n + 1)
  Phi_s_grenzen[2:n] <- (Phi_s[2:n] + Phi_s[1:(n-1)]) / 2
  Phi_s_grenzen[1] <- Phi_s[1]
  Phi_s_grenzen[n+1] <- Phi_s[n]
  
  # ∂/∂x(Phi_s*∂c/∂x)
  d_Phi_dcdx <- (Phi_s_grenzen[2:(n+1)] * dcdx[2:(n+1)] - 
                   Phi_s_grenzen[1:n] * dcdx[1:n]) / dx
  
  # ∂c/∂t = (D/Phi_s) * ∂/∂x(Phi_s*∂c/∂x)
  dcdt <- (Db / Phi_s) * d_Phi_dcdx
  
  list(dcdt)
}

dx <- x[2] - x[1]
parms <- list(Phi_s = Phi_s, Db = Db, dx = dx)

loesung <- ode(y = c0, times = c(0, 30), func = diffusion, parms = parms, method = "bdf")
result_intraphase_mixing <- data.frame(c_t1 = loesung[nrow(loesung), -1], Tiefe = x)
plot(result_intraphase_mixing, type = "l", main = "Intraphase porosity mixing", ylim = c(15, 0))


## interphase mixing: ∂Phi_s*c/∂t = Db * ∂/∂x(∂Phi_s*c/∂x) ##

x <- seq(0, 15, by = 0.5)
results_solid_fraction_mud <- solid_fraction_gradient(x)
Phi_s <- results_solid_fraction_mud

diffusion <- function(t, c, parms) {
  n <- length(c)
  Phi_s <- parms$Phi_s
  Db <- parms$Db
  dx <- parms$dx
  
  # Phi_s * c
  Phi_c <- Phi_s * c
  
  # ∂(Phi_s*c)/∂x an den Zellgrenzen
  d_Phi_c <- numeric(n + 1)
  d_Phi_c[2:n] <- (Phi_c[2:n] - Phi_c[1:(n-1)]) / dx
  d_Phi_c[1] <- 0
  d_Phi_c[n+1] <- 0
  
  # ∂/∂x(∂(Phi_s*c)/∂x)
  d2_Phi_c <- (d_Phi_c[2:(n+1)] - d_Phi_c[1:n]) / dx
  
  # ∂(Phi_s*c)/∂t = Db * ∂²(Phi_s*c)/∂x²
  # → ∂c/∂t = (Db/Phi_s) * ∂²(Phi_s*c)/∂x²
  dcdt <- (Db / Phi_s) * d2_Phi_c
  
  list(dcdt)
}

dx <- x[2] - x[1]

c0 <- c(100, rep(0, length(x) - 1))
parms <- list(Phi_s = Phi_s, Db = Db, dx = dx)

loesung <- ode(y = c0, times = c(0, 30), diffusion, parms, method = "bdf")
result_interphase_mixing <- data.frame(c_t30 = loesung[nrow(loesung), -1], Tiefe = x)
plot(result_interphase_mixing, type = "l", main = "Interphase porosity mixing", ylim = c(15, 0))


# combinated results
par(cex = 1.5)
plot(result_constant_porosity, type = "l", ylim = c(15, 0), lwd = 2, xlab = "", ylab = "Tiefe in cm")
lines(result_intraphase_mixing, ylim = c(15, 0), lty = 2, lwd = 2)
lines(result_interphase_mixing, ylim = c(15, 0), lty = 3, lwd = 2)
legend("bottom", c("constant porosity mixing", "intraphase mixing", "interphase mixing"), lty = c(1, 2, 3), lwd = 2)



### non-local exchange with regard to porosity gradients (DIY) ####

# intraphase mixing

Phi_s <- results_solid_fraction_mud

nonlocal_intraphase <- function(
    k=0.031, u=0.001, Db=0.1, depo=0.1,
    L=5, injectflux=0.3, sed)
{
  a <- (u/Db_eff - sqrt ( (u/Db_eff)^2 + 4*k / Db_eff)) /2
  b <- (u/Db_eff + sqrt ( (u/Db_eff)^2 + 4*k / Db_eff)) /2
  expaL <- exp(a*L)
  expbL <- exp(b*L)
  
  A <- matrix(nrow = 3, ncol = 3, byrow = TRUE, data = c(
    
    u-Db_eff*a, u-Db_eff*b, 0,
    expaL, expbL, -expaL,
    (u-Db_eff*a)*expaL, (u-Db_eff*b)*expbL, (Db_eff*a-u)*expaL)
  )
  B <- c(depo, 0, -injectflux)
  X <- solve(A, B)
  
  s1 <- which(sed<L)
  s2 <- which(sed>=L)
  
  conc <- vector(length = length(sed))
  conc[s1] <- X[1]*exp(a*sed[s1])+X[2]*exp(b*sed[s1])
  conc[s2] <- X[3]*exp(a*sed[s2])
  return(conc)
}

depth <- seq(0, 15, by = 0.5) # cm
Pb_intraphase <- rep(0, 31)
for (i in depth) {
  Db_eff <- Db*Phi_s[i*2+1]
  Pb_new <- nonlocal_intraphase(sed=i)
  Pb_intraphase[i*2+1] <- Pb_new
}

# interphase mixing

Phi_s <- results_solid_fraction_mud

nonlocal_interphase <- function(
    k=0.031, u=0.001, Db=0.1, depo=0.1,
    L=5, injectflux=0.3, sed)
{
  a <- (u/Db_eff - sqrt ( (u/Db_eff)^2 + 4*k / Db_eff)) /2
  b <- (u/Db_eff + sqrt ( (u/Db_eff)^2 + 4*k / Db_eff)) /2
  expaL <- exp(a*L)
  expbL <- exp(b*L)
  
  A <- matrix(nrow = 3, ncol = 3, byrow = TRUE, data = c(
    
    u-Db_eff*a, u-Db_eff*b, 0,
    expaL, expbL, -expaL,
    (u-Db_eff*a)*expaL, (u-Db_eff*b)*expbL, (Db_eff*a-u)*expaL)
  )
  B <- c(depo, 0, -injectflux)
  X <- solve(A, B)
  
  s1 <- which(sed<L)
  s2 <- which(sed>=L)
  
  conc <- vector(length = length(sed))
  conc[s1] <- X[1]*exp(a*sed[s1])+X[2]*exp(b*sed[s1])
  conc[s2] <- X[3]*exp(a*sed[s2])
  return(conc)
}

depth <- seq(0, 15, by = 0.5) # cm
Pb_interphase <- rep(0, 31)
for (i in depth) {
  Db_eff <- Db*Phi_s[i*2+1]
  Pb_new <- nonlocal_interphase(sed=i)
  Pb_interphase[i*2+1] <- Pb_new
}

# multi plot non-local with porosity gradient
par(cex = 1.5)
depth <- seq(0, 15, by = 0.5)
plot(Pb_intraphase, depth, ylim = c(15, 0), xlab = "Pb, dpm/cm2/yr",
     ylab = "Tiefe in cm", type = "l", lwd = 2, lty = 2)
depth <- seq(0, 15, by = 0.01)
lines(Pb, depth, lty = 1, lwd = 2) # assumes constant porosity
# lines(Pb_interphase, depth, lty = 3)
legend("bottom", c("constant porosity mixing", "intraphase mixing"), lty = c(1, 2), lwd = 2)



