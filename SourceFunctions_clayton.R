##################### Source Function for pMLE: clayton ##################### 


# Function to return estimated SD, S1star, S2star
npest.star_12<-function (X){
  
  k <- ncol(X)
  
  fit.SD <- survival::survfit(survival::Surv(X$C, X$delta_D)~1)
  SD.hat.fn  <- approxfun(x = fit.SD$time, y = fit.SD$surv, yleft = 1, yright = 0, method = "constant")
  SD.hat <- SD.hat.fn(X$C)
  
  fit.S1.star <- survival::survfit(survival::Surv(X$U1, X$eta_1)~1)
  fit.S2.star <- survival::survfit(survival::Surv(X$U2, X$eta_2)~1)
  S1.star.fn  <- approxfun(x = fit.S1.star$time, y = fit.S1.star$surv, yleft = 1, yright = 0,method = "constant")
  S2.star.fn  <- approxfun(x = fit.S2.star$time, y = fit.S2.star$surv, yleft = 1, yright = 0,method = "constant")
  
  S1.star.hat <- S1.star.fn(X$U1)
  S2.star.hat <- S2.star.fn(X$U2)
  
  #True S1/S1 valued at observed U1/U1 for MLE
  S1.fn <- approxfun(x = X$T1, y = X$S1, method = "constant",yleft = 1, yright = 0)
  S2.fn <- approxfun(x = X$T2, y = X$S2, method = "constant",yleft = 1, yright = 0)
  
  S1.U1 <- S1.fn(X$U1)
  S2.U2 <- S2.fn(X$U2)
  
  #S1.U1 <- pweibull(X$U1, shape = 2, scale = 70, lower.tail = F)
  
  #True SD valued at observed C
  SD.fn <- approxfun(x = X$D, y = X$SD, method = "constant", yleft = 1, yright = 0)
  SD.C <- SD.fn(X$C)
  #SD.C  <- pweibull(X$C, shape = 2, scale = 85, lower.tail = F)
  
  #Calculate SD.hat at values of U1/U2
  SD.U1.hat <- SD.hat.fn(X$U1)
  SD.U2.hat <- SD.hat.fn(X$U2)
  
  #"True" SD at values of U1/U2
  SD.U1 <- SD.fn(X$U1)
  SD.U2 <- SD.fn(X$U2)
  
  
  #store S and f for likelihood calc later
  X[,k+1] <- S1.star.hat; X[,k+2] <- SD.hat
  X[,k+3] <- S1.U1; X[,k+4] <- SD.C
  X[,k+5] <- SD.U1.hat; X[,k+6] <- SD.U1
  X[,k+7] <- S2.star.hat; X[,k+8] <- S2.U2
  X[,k+9] <- SD.U2.hat; X[,k+10] <- SD.U2

  
  colnames(X)[(k+1):(k+10)] <- c("S1.star.hat", "SD.hat", "S1.U1", "SD.C", "SD.U1.hat", "SD.U1",
                                 "S2.star.hat", "S2.U2", "SD.U2.hat",
                                 "SD.U2")
  
  return( list(X= X, SD.hat.fn = SD.hat.fn,
                S1.star.fn= S1.star.fn, S2.star.fn= S2.star.fn))
}



####################################################################################################
#                                    Clayton Copula Likelihood Function                            #
####################################################################################################

clayton.cc.stage.T1T2.D.lik <- function(X0, thetas, S1.hat, S2.hat, SD.hat)
{
  #  theta = exp(b)-1
  theta = thetas[1]
  theta0 = thetas[2]
  
  X <- data.frame(X0, theta, theta0)
  X$S1 = S1.hat
  X$S2 = S2.hat
  X$SD = SD.hat
  
  case1.1 <- X[which(X$delta_1 == 1 & X$delta_2  == 1 & X$delta_D == 1),]
  if(nrow(case1.1)!=0){L.part1.1 <- with(case1.1, clayton.C21(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C.1(S1, S2, theta)*
                                           clayton.C1.(S1, S2, theta)+
                                           clayton.C11(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C11(S1, S2, theta))}else{L.part1.1 = 0}
  
  case1.2 <- X[which(X$delta_1 == 1 & X$delta_2 == 0 & X$delta_D == 1),]
  if(nrow(case1.2)!=0){L.part1.2 <- with(case1.2, clayton.C11(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C1.(S1, S2, theta))}else{L.part1.2 = 0}
  
  case1.3 <- X[which(X$delta_1 == 0 & X$delta_2 == 1 & X$delta_D == 1),]
  if(nrow(case1.3)!=0){L.part1.3 <- with(case1.3, clayton.C11(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C.1(S1, S2, theta))}else{L.part1.3 = 0}
  
  case1.4 <- X[which(X$delta_1 == 0 & X$delta_2 == 0 & X$delta_D == 1),]
  if(nrow(case1.4)!=0){L.part1.4 <- with(case1.4, clayton.C.1(clayton.C..(S1, S2, theta), SD, theta0))}else{L.part1.4 = 0}
  
  case2.1 <- X[which(X$delta_1 == 1 & X$delta_2  == 1 & X$delta_D == 0),]
  if(nrow(case2.1)!=0){L.part2.1 <- with(case2.1, clayton.C2.(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C.1(S1, S2, theta)*
                                           clayton.C1.(S1, S2, theta)+
                                           clayton.C1.(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C11(S1, S2, theta))}else{L.part2.2 = 0}
  
  case2.2 <- X[which(X$delta_1 == 1 & X$delta_2 == 0 & X$delta_D == 0),]
  if(nrow(case2.2)!=0){L.part2.2 <-  with(case2.2, clayton.C1.(clayton.C..(S1, S2, theta), SD, theta0)*
                                            clayton.C1.(S1, S2, theta))}else{L.part2.2 = 0}
  
  case2.3 <- X[which(X$delta_1 == 0 & X$delta_2 == 1 & X$delta_D == 0),]
  if(nrow(case2.3)!=0){L.part2.3 <- with(case2.3, clayton.C1.(clayton.C..(S1, S2, theta), SD, theta0)*
                                           clayton.C.1(S1, S2, theta))}else{L.part2.3 = 0}
  
  case2.4 <- X[which(X$delta_1 == 0 & X$delta_2 == 0 & X$delta_D == 0),]
  if(nrow(case2.4)!=0){L.part2.4 <-with(case2.4, clayton.C..(clayton.C..(S1, S2, theta), SD, theta0))}else{L.part2.4 = 0}
  
  L.result <- c(L.part1.1, L.part1.2, L.part1.3, L.part1.4, 
                L.part2.1, L.part2.2, L.part2.3, L.part2.4)
  
  suppressWarnings(-sum(log(L.result[L.result != 0]), na.rm = TRUE))  
  
}

####################################################################################################
#                                    Clayton Copula SourceFunctions                                #
####################################################################################################

cdfExpr.Clayton <- function(n) {
  expr <- "u1^(-theta) - 1"
  for (i in 2:n) {
    cur <- paste0("u", i, "^(-theta) - 1")
    expr <- paste(expr, cur, sep = " + ")
  }
  expr <- paste("(1 + (", expr, "))^ (-1/theta)")
  parse(text = expr)
}
pdfExpr <- function(cdf, n) {
  val <- cdf
  for (i in 1:n) {
    val <- D(val, paste0("u", i))
  }
  val
}

clayton.f <- function(theta)
{
  paste(clayton.C11(u1, u2, theta))
  
}


C.clayton.expr <- cdfExpr.Clayton(2)
#C.clayton <- expression((u^(-theta) + v^(-theta) -1)^(-1/theta))
clayton.C1. <- function(u1,u2,theta){}; body(clayton.C1.) = D(C.clayton.expr, "u1")
clayton.C2. <- function(u1,u2,theta){}; body(clayton.C2.) = D(D(C.clayton.expr, "u1"),"u1")
clayton.C.1 <- function(u1,u2,theta){}; body(clayton.C.1) = D(C.clayton.expr, "u2")
clayton.C11 <- function(u1,u2,theta){}; body(clayton.C11) = D(D(C.clayton.expr, "u1"), "u2")
clayton.C21 <- function(u1,u2,theta){}; body(clayton.C21) = D(D(D(C.clayton.expr, "u1"), "u1"), "u2")
#clayton.C11 <- pdfExpr(C.clayton.expr, 2)
clayton.C.. <- function(u1,u2,theta){}; body(clayton.C..) = C.clayton.expr

clayton.dtheta.C11 <- function(u1,u2,theta){}; body(clayton.dtheta.C11) = D(D(D(C.clayton.expr, "u1"), "u2"), "theta")
clayton.dtheta.C.1 <- function(u1, u2, theta){}; body(clayton.dtheta.C.1) = D(D(C.clayton.expr, "u2"), "theta")
clayton.dtheta.C1. <- function(u1, u2, theta){}; body(clayton.dtheta.C1.) = D(D(C.clayton.expr, "u1"), "theta")
clayton.dtheta.C.. <- function(u1, u2, theta){}; body(clayton.dtheta.C..) = D(C.clayton.expr, "theta")

clayton.dtheta2.C11 <- function(u1,u2,theta){}; body(clayton.dtheta2.C11) = D(D(D(D(C.clayton.expr, "u1"), "u2"), "theta"), "theta")
clayton.dtheta2.C1. <- function(u1,u2,theta){}; body(clayton.dtheta2.C1.) = D(D(D(C.clayton.expr, "u1"), "theta"),"theta")
clayton.dtheta2.C.1 <- function(u1,u2,theta){}; body(clayton.dtheta2.C.1) = D(D(D(C.clayton.expr, "u2"), "theta"),"theta")
clayton.dtheta2.C.. <- function(u1,u2,theta){}; body(clayton.dtheta2.C..) = D(D(C.clayton.expr, "theta"),"theta")


C.clayton.expr.3d <- cdfExpr.Clayton(3)
#C.clayton <- expression((u^(-theta) + v^(-theta) -1)^(-1/theta))
clayton.C111 <- function(u1,u2,u3,theta){}; body(clayton.C111) = D(D(D(C.clayton.expr.3d, "u1"),"u2"),"u3")
clayton.C101 <- function(u1,u2,u3,theta){}; body(clayton.C101) = D(D(C.clayton.expr.3d, "u1"),"u3")
clayton.C011 <- function(u1,u2,u3,theta){}; body(clayton.C011) = D(D(C.clayton.expr.3d, "u2"), "u3")
clayton.C001 <- function(u1,u2,u3,theta){}; body(clayton.C001) = D(C.clayton.expr.3d,"u3")

clayton.C110 <- function(u1,u2,u3,theta){}; body(clayton.C110) = D(D(C.clayton.expr.3d, "u1"),"u2")
clayton.C100 <- function(u1,u2,u3,theta){}; body(clayton.C100) = D(C.clayton.expr.3d, "u1")
clayton.C010 <- function(u1,u2,u3,theta){}; body(clayton.C010) = D(C.clayton.expr.3d, "u2")
clayton.C000 <- function(u1,u2,u3,theta){}; body(clayton.C000) = C.clayton.expr.3d



g.fn <- function(w1,w2,theta)
{
  (w1^(-theta) - w2^(-theta) + 1)^(-1/theta)
  
}

g.expression <- expression( (w1^(-theta) - w2^(-theta) + 1)^(-1/theta))

g1. <- function(w1,w2,theta){}; body(g1.) = D(g.expression,"w1")
g.1 <- function(w1,w2,theta){}; body(g.1) = D(g.expression,"w2")
dtheta.g <- function(w1,w2,theta){}; body(dtheta.g) = D(g.expression,"theta")
dw1.g <- function(w1,w2,theta){}; body(dw1.g) = D(g.expression,"w1")
dw2.g<- function(w1,w2,theta){}; body(dw2.g) = D(g.expression,"w2")

dtheta.g1. <-function(w1,w2,theta){}; body(dtheta.g1.) = D((D(g.expression,"w1")), "theta")
dtheta2.g <- function(w1,w2,theta){}; body(dtheta2.g) = D(D(g.expression,"theta"), "theta")
dtheta2.g1. <- function(w1,w2,theta){}; body(dtheta2.g1.) = D(D((D(g.expression,"w1")), "theta"),"theta")

algebra.X <- function(X, theta)
{
  k <- ncol(X)
  X[,k+1] <- g.fn(w1 = X$S1.star.hat, w2 =  X$SD.U1.hat, theta = theta)
  colnames(X)[k+1] = "S1.hat"
  
  return (X)
}

algebra.X_12 <- function(X, theta1, theta2)
{
  k <- ncol(X)
  X[,k+1] <- g.fn(w1 = X$S1.star.hat, w2 =  X$SD.U1.hat, theta = theta1)
  colnames(X)[k+1] = "S1.hat"
  
  X[,k+2] <- g.fn(w1 = X$S2.star.hat, w2 =  X$SD.U2.hat, theta = theta2)
  colnames(X)[k+2] = "S2.hat"
  return (X)
}

algebra.X_2 <- function(X, theta2)
{
  k <- ncol(X)
  
  X[,k+1] <- g.fn(w1 = X$S2.star.hat, w2 =  X$SD.U2.hat, theta = theta2)
  colnames(X)[k+1] = "S2.hat"
  return (X)
}


clayton.logL.expr <- expression (d1*dD * log((-1)*((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^(((-1/theta) - 1) - 1) * (((-1/theta) - 1) * (u2^((-theta) - 1) * (-theta))) * ((-1/theta) * (u1^((-theta) - 1) * (-theta))))*((w1^(-theta) - w2^(-theta) + 1)^((-1/theta) - 1) * ((-1/theta) * (w1^((-theta) - 1) * (-theta)))))
                                 + d1*(1-dD)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^((-1/theta) - 1) *  ((-1/theta) * (u1^((-theta) - 1) * (-theta)))*((w1^(-theta) - w2^(-theta) + 1)^((-1/theta) - 1) * ((-1/theta) * (w1^((-theta) - 1) * (-theta)))))
                                 + (1-d1)*dD*log((-1)*(1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^((-1/theta) - 1) * ((-1/theta) * (u2^((-theta) - 1) * (-theta))))
                                 + (1-d1)*(1-dD)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^(-1/theta))
)

clayton.logL.expr_12 <- expression (d1*d2 *((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^(((-1/theta) - 1) - 
                                                                                         1) * (((-1/theta) - 1) * (u2^((-theta) - 1) * (-theta))) * 
                                              ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
                                    +d1*(1-d2)*((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^((-1/theta) - 1) * 
                                                  ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
                                    +(1-d1)*d2*((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^((-1/theta) - 1) * 
                                                  ((-1/theta) * (u2^((-theta) - 1) * (-theta))))
                                    +(1-d1)*(1-d2)*((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^(-1/theta))
                                    
)

clayton.logL <- function(d1, dD, u1, u2, w1, w2, theta){};body(clayton.logL) = clayton.logL.expr
clayton.dtheta.logL <- function(d1, dD, u1, u2, w1, w2, theta){}; body(clayton.dtheta.logL) = D(clayton.logL.expr, "theta")
clayton.dtheta2.logL <- function(d1, dD, u1, u2, w1, w2, theta){}; body(clayton.dtheta2.logL) = D(D(clayton.logL.expr, "theta"),"theta")

clayton.logL_12 <- function(d1, d2, u1, u2,theta){};body(clayton.logL_12) = clayton.logL.expr_12
clayton.dtheta.logL_12 <- function(d1, d2, u1, u2,  theta){}; body(clayton.dtheta.logL_12) = D(clayton.logL.expr_12, "theta")
clayton.dtheta2.logL_12 <- function(d1, d2, u1, u2, theta){}; body(clayton.dtheta2.logL_12) = D(D(clayton.logL.expr_12, "theta"),"theta")



clayton.logL.expr.naive <- expression(d1*dD*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^(((-1/theta) - 1) - 
                                                                                             1) * (((-1/theta) - 1) * (u2^((-theta) - 1) * (-theta))) * 
                                                  ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
                                      +d1*(1-dD)*log((-1)*((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^((-1/theta) - 1) * 
                                                             ((-1/theta) * (u1^((-theta) - 1) * (-theta)))))
                                      +(1-d1)*dD*log((-1)*((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^((-1/theta) - 1) * 
                                                             ((-1/theta) * (u2^((-theta) - 1) * (-theta)))))
                                      +(1-d1)*(1-dD)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1))^(-1/theta))
)

clayton.dtheta.logL.naive <- function(d1, dD, u1, u2, theta){}; body(clayton.dtheta.logL.naive) = D(clayton.logL.expr.naive, "theta")
clayton.dtheta2.logL.naive <- function(d1, dD, u1, u2, theta){}; body(clayton.dtheta2.logL.naive) = D(D(clayton.logL.expr.naive, "theta"),"theta")

#for calculation of sandwich variance 
clayton.logL.expr.T1T2D <- expression(
  d1*d2*d3*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^((((-1/theta) - 
                                                                                1) - 1) - 1) * ((((-1/theta) - 1) - 1) * (u3^((-theta) - 
                                                                                                                                1) * (-theta))) * (((-1/theta) - 1) * (u2^((-theta) - 1) * 
                                                                                                                                                                         (-theta))) * ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
  +d1*(1-d2)*d3*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^(((-1/theta) - 
                                                                                    1) - 1) * (((-1/theta) - 1) * (u3^((-theta) - 1) * (-theta))) * 
                      ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
  +(1-d1)*d2*d3*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^(((-1/theta) - 
                                                                                    1) - 1) * (((-1/theta) - 1) * (u3^((-theta) - 1) * (-theta))) * 
                      ((-1/theta) * (u2^((-theta) - 1) * (-theta))))
  +(1-d1)*(1-d2)*d3*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^((-1/theta) - 
                                                                                       1) * ((-1/theta) * (u3^((-theta) - 1) * (-theta))))
  +d1*d2*(1-d3)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^(((-1/theta) - 
                                                                                    1) - 1) * (((-1/theta) - 1) * (u2^((-theta) - 1) * (-theta))) * 
                      ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
  +d1*(1-d2)*(1-d3)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^((-1/theta) - 
                                                                                       1) * ((-1/theta) * (u1^((-theta) - 1) * (-theta))))
  +(1-d1)*d2*(1-d3)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^((-1/theta) - 
                                                                                       1) * ((-1/theta) * (u2^((-theta) - 1) * (-theta))))
  +(1-d1)*(1-d2)*(1-d3)*log((1 + (u1^(-theta) - 1 + u2^(-theta) - 1 + u3^(-theta) - 1))^(-1/theta))
  
)
clayton.dtheta.logL.T1T2D <- function(d1,d2,d3, u1,u2, u3, theta){};body(clayton.dtheta.logL.T1T2D)=D(clayton.logL.expr.T1T2D,"theta")
clayton.dtheta2.logL.T1T2D <- function(d1,d2,d3, u1,u2,u3, theta){};body(clayton.dtheta2.logL.T1T2D)=D(D(clayton.logL.expr.T1T2D ,"theta"), "theta")


