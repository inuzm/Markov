#---------------------------------------------------------------
# 0) SOURCE FILES
#---------------------------------------------------------------
# The directory to run the simulations must have the following 
# structure
#
# Directory structure:
#
# dir 
#  |
#  |- cpp 
#      |- Orntein_Uhlenbeck.cpp
#      |- Feller_Harris.cpp
#  |- R
#      |- PlotSimulations.R
#  |- plots 


#Clear all
rm(list = ls())

#Call libraries
library("ggplot2")
library("gridExtra")
library("Rcpp")
library("latex2exp")

#Source cpp functions
sourceCpp("cpp/Orntein_Uhlenbeck.cpp")
sourceCpp("cpp/Feller_Harris.cpp")

#Source R functions
source("R/PlotSimulations.R")


#---------------------------------------------------------------
# 1) ORNSTEIN-UHLENBECK MODEL
#---------------------------------------------------------------

#Run model
nsim  <- 250
tspan <- 0:1000/10
mu    <- 10
tau   <- 1
alpha <- 5
simmat_direct <- rornstein.uhlenbeck(tspan, mu, tau, alpha, nsim, "direct")
simmat_mix    <- rornstein.uhlenbeck(tspan, mu, tau, alpha, nsim, "mixture")
simmat_sde    <- rornstein.uhlenbeck(tspan, mu, tau, alpha, nsim, "sde")

#Plot
pdf("plots/Direct_OU.pdf", width = 10, height = 5)
plotsims(tspan, simmat_direct, mytitle = paste0("Ornstein-Uhlenbeck under direct method simulation (", nsim, " simulations)."))
dev.off()

pdf("plots/Mixture_OU.pdf", width = 10, height = 5)
plotsims(tspan, simmat_mix, mytitle = paste0("Ornstein-Uhlenbeck under mixture method simulation (", nsim, " simulations)."))
dev.off()

pdf("plots/SDE_OU.pdf", width = 10, height = 5)
plotsims(tspan, simmat_sde, mytitle = paste0("Ornstein-Uhlenbeck under sde discretization method simulation (", nsim, " simulations)."))
dev.off()

#Timing of functions
if (require(microbenchmark)){
  outime <- microbenchmark(
    rornstein.uhlenbeck(tspan, mu, tau, alpha, nsim, "direct"),
    rornstein.uhlenbeck(tspan, mu, tau, alpha, nsim, "mixture"),
    rornstein.uhlenbeck(tspan, mu, tau, alpha, nsim, "sde")
  )
}

#---------------------------------------------------------------
# 2) FELLER-HARRIS MODEL
#---------------------------------------------------------------

#Run model
nsim  <- 250
tspan <- 0:1000/10
lambda_a  <- 0.4
theta_a   <- 0.3
alpha_a   <- lambda_a*theta_a/(1 + theta_a)

#Gamma simulations
#---------------------------------
shape <- 2
rate  <- 1
x0        <- rgamma(nsim, shape = shape, rate = rate)
simmat_pp <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "poissonprocess", dist = rgamma, shape = shape, rate = rate)
simmat_ck <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "chapmankolmogorov", dist = rgamma, shape = shape, rate = rate)

#Plot
pdf("plots/Gamma_PP_Feller-Harris.pdf", width = 10, height = 5)
plotsims(tspan, simmat_pp, mytitle = paste0("Feller-Harris under direct method simulation with Q gamma (", nsim, " simulations)."))
dev.off()

pdf("plots/Gamma_CK_Feller-Harris.pdf", width = 10, height = 5)
plotsims(tspan, simmat_ck, mytitle = paste0("Feller-Harris under direct method simulation with Q gamma (", nsim, " simulations)."))
dev.off()

#Normal simulations
#---------------------------------
mu      <- 10
desvest <- 0.25
x0        <- rnorm(nsim, mean = mu, sd = desvest)
simmat_pp <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "poissonprocess", dist = rnorm, mean = mu, sd = desvest)
simmat_ck <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "chapmankolmogorov", dist = rnorm, mean = mu, sd = desvest)

#Plot
pdf("plots/Normal_PP_Feller-Harris.pdf", width = 10, height = 5)
plotsims(tspan, simmat_pp, mytitle = paste0("Feller-Harris under direct method simulation with Q normal (", nsim, " simulations)."))
dev.off()

pdf("plots/Normal_CK_Feller-Harris.pdf", width = 10, height = 5)
plotsims(tspan, simmat_ck, mytitle = paste0("Feller-Harris under direct method simulation with Q normal (", nsim, " simulations)."))
dev.off()

#Poisson simulations
#---------------------------------
lambda    <- 1.5
x0        <- rpois(nsim, lambda = lambda)
#Falla al poner 
#simmat_pp <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "poissonprocess", dist = rpois, lambda = lambda)
simmat_pp <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "poissonprocess", dist = function(n){rpois(n, lambda)})
simmat_ck <- rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "chapmankolmogorov", dist = function(n){rpois(n, lambda)})

#Plot
pdf("plots/Poisson_PP_Feller-Harris.pdf", width = 10, height = 5)
plotsims(tspan, simmat_pp, mytitle = paste0("Feller-Harris under direct method simulation with Q gamma (", nsim, " simulations)."))
dev.off()

pdf("plots/Poisson_CK_Feller-Harris.pdf", width = 10, height = 5)
plotsims(tspan, simmat_ck, mytitle = paste0("Feller-Harris under direct method simulation with Q gamma (", nsim, " simulations)."))
dev.off()

#Timing of functions
if (require(microbenchmark)){
  fellertime <- microbenchmark(
    rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "poissonprocess", dist = rnorm, mean = mu, sd = desvest),
    rfeller.harris(tspan, lambda_a, theta_a, alpha_a, nsim, x0, type = "chapmankolmogorov", dist = rnorm, mean = mu, sd = desvest)
  )
}
