##################### pseudo-MLE estimation procedure: nested clayton ##################### 
rm(list=ls())
setwd("./clayton_example/")
library(dplyr)
library(survival)
library(copula)
library(plyr)
library(muhaz)
library(mvtnorm)
library(cubature)

source("SourceFunctions_clayton.R")

# Setting with Paramter - Clayton copula
tau = 0.5  # value of tau12 in paper
tau0 = 0.4 # value of tau in paper
n = 2000 # sample size
theta.c = 2*tau/(1-tau) # clayton inner theta
theta0.c = 2*tau0/(1-tau0) # clayton outer theta

# Read Generated Data
mydat = readRDS("mydat.RDS")


# Return estimated SD and S1star, S1star and store in mydat1
out = npest.star_12(mydat)
mydat1 = out$X

#Estimation Procedure
#### outer: clatyon, inner: clayton ####
iter = 1; K = 50
procC.iter.clayton.c.theta0 <- matrix(NA, nrow = K, ncol = 1)
procC.iter.clayton.c.theta  <- matrix(NA, nrow = K, ncol = 1)

procC.iter.clayton.c.tau0 <- matrix(NA, nrow = K, ncol = 1)
procC.iter.clayton.c.tau  <- matrix(NA, nrow = K, ncol = 1)

procC.iter.clayton.c.theta0[iter,1] <- 1
procC.iter.clayton.c.theta[iter,1] <-  1

procC.iter.clayton.c.tau0[iter,1] <-  procC.iter.clayton.c.theta0[iter,1]/(procC.iter.clayton.c.theta0[iter,1]+2)
procC.iter.clayton.c.tau[iter,1] <-   procC.iter.clayton.c.theta[iter,1]/(procC.iter.clayton.c.theta[iter,1]+2)


print(paste(iter, "run:", "theta0:", procC.iter.clayton.c.theta0[iter,1], ";theta:", procC.iter.clayton.c.theta[iter,1] ))

repeat{
  
  iter = iter + 1
  if(iter == K){break}
  
  procC.iter.clayton.c.theta0[iter,1]  = procC.iter.clayton.c.theta0[(iter-1),1]
  procC.iter.clayton.c.theta[iter,1]   = procC.iter.clayton.c.theta[(iter-1),1]
  
  procC.iter.clayton.c.tau0[iter,1]  = procC.iter.clayton.c.tau0[(iter-1),1]
  procC.iter.clayton.c.tau[iter,1]   = procC.iter.clayton.c.tau[(iter-1),1]
  
  c.S1.hat = g.fn(mydat1$S1.star.hat, mydat1$SD.U1.hat, procC.iter.clayton.c.theta0[iter,1])
  c.S2.hat = g.fn(mydat1$S2.star.hat, mydat1$SD.U2.hat, procC.iter.clayton.c.theta0[iter,1])
  
  keep_ind = which(between(c.S1.hat,0.001,0.999) & between(c.S2.hat,0.001,0.999) & !is.na(c.S1.hat) & !is.na(c.S2.hat))
  c.S1.hat.keep = c.S1.hat[keep_ind]
  c.S2.hat.keep = c.S2.hat[keep_ind]
  mydat1.keep = mydat1[keep_ind,]
  
  out = optim(c(procC.iter.clayton.c.theta[iter,1],procC.iter.clayton.c.theta0[iter,1]), 
              clayton.cc.stage.T1T2.D.lik, X0 = mydat1.keep,
              S1.hat = c.S1.hat.keep, S2.hat = c.S2.hat.keep,
              SD.hat = mydat1.keep$SD.hat
  )

  procC.iter.clayton.c.theta0[iter,1] <- out$par[2]
  procC.iter.clayton.c.theta[iter,1] <- out$par[1]
  
  
  print(paste(iter, "run:", "theta:", procC.iter.clayton.c.theta0[iter,1], ";theta12:", procC.iter.clayton.c.theta[iter,1] ))
  print(paste(iter, "run:", "tau:", procC.iter.clayton.c.tau0[iter,1], ";tau12:", procC.iter.clayton.c.tau[iter,1] ))
  
  procC.iter.clayton.c.tau0[iter,1] = procC.iter.clayton.c.theta0[iter,1]/(procC.iter.clayton.c.theta0[iter,1]+2)
  procC.iter.clayton.c.tau[iter,1] = procC.iter.clayton.c.theta[iter,1]/(procC.iter.clayton.c.theta[iter,1]+2)
 
  if(abs(procC.iter.clayton.c.theta0[iter,1] - procC.iter.clayton.c.theta0[(iter-1),1] ) < 1e-5 &
     abs(procC.iter.clayton.c.theta[iter,1] - procC.iter.clayton.c.theta[(iter-1),1] )< 1e-5 ) break
  
}

