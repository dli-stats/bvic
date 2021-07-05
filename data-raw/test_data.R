## code to prepare `test_data` dataset goes here

##################### Data generation: nested clayton ##################### 

rm(list=ls())
# setwd("./clayton_example")
library(dplyr)
library(survival)
library(copula)
library(plyr)
library(muhaz)
library(mvtnorm)
library(cubature)


# Nested Clayton copula
# type = "cc" # "cc" or "ff"

# tau and tau12 parameter
tau = 0.5  # 0.5 or 0.8
tau0 = 0.4 # 0.4 or 0.3

#sample size 
n = 2000

theta.c = 2*tau/(1-tau) 
theta0.c = 2*tau0/(1-tau0)


sam = rnacopula(n, onacopula("Clayton", C(theta0.c, 3,C(theta.c, c(1,2)))))


T1 = qweibull(sam[,1], shape = 2, scale = 70, lower.tail = F)
T2 = qweibull(sam[,2], shape = 2, scale = 60, lower.tail = F)
D  = qweibull(sam[,3], shape = 2, scale = 85, lower.tail = F)
Ca = rexp(nrow(sam), rate = 1/350)

C = pmin(D, Ca)
U1 = pmin(T1, C)
U2 = pmin(T2, C)

T1star = pmin(T1, D)
T2star = pmin(T2, D)

eta_1 = ifelse(T1star <= Ca, 1, 0)
eta_2 = ifelse(T2star <= Ca, 1, 0)

delta_D = ifelse(D  <= Ca, 1, 0)
delta_1 = ifelse(T1 <= C, 1, 0)
delta_2 = ifelse(T2 <= C, 1, 0)


# S1 <- sam[,1]
# S2 <- sam[,2]
# SD <- sam[,3]
# 
# S1.U1 <- pweibull(U1, shape = 2, scale = 70, lower.tail = F)
# S2.U2 <- pweibull(U2, shape = 2, scale = 60, lower.tail = F)
# SD.C  <- pweibull(C,  shape = 2, scale = 85, lower.tail = F)
# SD.U1 <- pweibull(U1, shape = 2, scale = 85, lower.tail = F)
# SD.U2 <- pweibull(U2, shape = 2, scale = 85, lower.tail = F)

theta = theta.c

test_data = data.frame(U1, U2, C, 
                       # T1, T2,
                   # T1star, T2star, D, 
                   # Ca, 
                   delta_1, delta_2, delta_D
                   # eta_1, eta_2
                   # S1, S2, SD, S1.U1, S2.U2, SD.C, SD.U1, SD.U2
                   )


usethis::use_data(test_data, overwrite = TRUE)
