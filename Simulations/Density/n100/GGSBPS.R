#------------------------------------------------------------------------------#
#               Griddy-Gibbs routine for Poisson model                         #
#         Copyright, Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

GriddyGibbs <- function(x, y, niter = 100, penorder = 2, 
                        burn = "half", xl = min(x), xr = max(x), K = 10,
                        progressbar = TRUE){

n <- length(y)
B <- Rcpp_KercubicBspline(x, lower = xl, upper = xr, K = K) # C++ call
D <- diag(K)
for(k in 1:penorder){D <- diff(D)}
P <- t(D) %*% D           
rankP <- K-penorder

# Gamma prior on penalty parameter lambda > 0.
a_lambda <- 1e-04
b_lambda <- 1e-04
lambda_shape <- 0.5 * rankP + a_lambda

# Define quantities for compact notation
psik <- colSums(y*B)

if(penorder == 2){
Sk <- function(k, theta, lambda){
  if(k == 1){
    lambfactor <- (2 * theta[k + 1] - theta[k + 2]) 
  } else if(k == 2){
    lambfactor <- (2 * theta[k - 1] + 4 * theta[k + 1] - theta[k + 2])
  } else if (k >= 3 & k <= (K-2)){
    lambfactor <- 4 * (theta[k-1] + theta[k+1])-(theta[k-2]+theta[k+2])
  } else if (k == (K-1)){
    lambfactor <- 4*theta[k-1]-theta[k-2]+2*theta[k+1]
  } else if (k == K){
    lambfactor <- 2*theta[k-1]-theta[k-2]
  }
  val <- psik[k] + lambda * lambfactor
  return(val)
}

Ilamb <- function(k, lambda){
  if(k == 1 | k == K){
    val <- lambda
  } else if(k == 2 | k == (K-1)){
    val <- 5 * lambda
  } else if (k >= 3 & k <= (K-2)){
    val <- 6 * lambda
  }
  return(val)
}
} else if(penorder == 3){
  Sk <- function(k, theta, lambda){
    if(k == 1){
      lambfactor <- (theta[k + 3] + 3 * theta[k+1] - 3 * theta[k+2]) 
    } else if(k == 2){
      lambfactor <- (-6 * theta[k + 2] + theta[k + 3] + 3 * theta[k - 1] + 12 * theta[k + 1])
    } else if (k == 3){
      lambfactor <- (15 * theta[k+1] - 6 * theta[k+2] + theta[k+3] -
                     3 * theta[k-2] + 12 * theta[k-1])
    } else if (k >= 4 & k <= (K-3)){
      lambfactor <- 15 * (theta[k-1] + theta[k+1]) - 
        6 * (theta[k-2] + theta[k+2]) + (theta[k-3] + theta[k+3])
    } else if (k == (K-2)){
      lambfactor <- 15 * (theta[k-1] + theta[k+1]) - 
        6 * (theta[k-2] + theta[k+2]) - 3 * theta[k+1] + 3 * theta[k+2]
    } else if (k == (K-1)){
      lambfactor <- 15 * (theta[k-1] + theta[k+1]) - 12 * theta[k+1] -
        3 * theta[k-1]
    } else if (k == K){
      lambfactor <- (-12 * theta[k-1] + 3 * theta[k-2])
    }
    val <- psik[k] + lambda * lambfactor
    return(val)
  }
  
  Ilamb <- function(k, lambda){
    if(k == 1 | k == K){
      val <- lambda
    } else if(k == 2 | k == (K-1)){
      val <- 10 * lambda
    } else if (k == 3 | k == (K-2)){
      val <- 19 * lambda 
    } else if (k >= 4 & k <= (K-3)){
      val <- 20 * lambda
    }
    return(val)
  }  
}

logp <- function(x, k, theta, lambda){
  xlen <- length(x)
  Xdes <- matrix(rep(x, each = n), ncol = xlen, byrow = FALSE)
  val <- (-0.5) * Ilamb(k, lambda) * x^2 + x * Sk(k, theta, lambda) -
    colSums((exp(B[,k]*Xdes) * as.numeric(exp(B[,-k]%*%theta[-k]))))
  return(val)
}

dlogp <- function(x, k, theta, lambda){
  val <- (-Ilamb(k, lambda) * x)+Sk(k, theta, lambda)-
    sum(exp(B[,k]*x) * B[,k] * exp(B[,-k]%*%theta[-k]))
  return(val)
}

d2logp <- function(x, k, theta, lambda){
  val <- (-1) * (Ilamb(k, lambda) +
    sum(exp(B[,k]*x) * (B[,k])^2 * exp(B[,-k]%*%theta[-k])))
  return(val)
}

# Starting values for B-spline coefficients
theta <- rep(0, K) # maximum a priori of B-splines amplitudes
lambda <- c()
thetamat <- matrix(0, nrow = niter, ncol = K)

if(isTRUE(progressbar)){
cat(paste0("Griddy-Gibbs algorithm in Bayesian P-splines Poisson model running for ",niter, " iterations \n"))
progbar <- utils::txtProgressBar(min = 1, max = niter, initial = 1,
                                 style = 3, char =">")
}
eps <- sqrt(.Machine$double.eps)

for(m in 1:niter){
  
  # Draw lambda
  lambda_rate <- as.numeric(0.5 * (t(theta)%*%P%*%theta)+b_lambda)
  lambda[m] <- rgamma(n = 1, shape = lambda_shape, rate = lambda_rate)
  
  # Draw thetas
  for(k in 1:K){
    
    # 1: Find maximum of conditional posterior
    deriv0 <- dlogp(x=0, k = k, theta = theta, lambda = lambda[m])
    if(deriv0 > 0){
      rbound <- (deriv0/Ilamb(k, lambda[m])) + eps
      xmax <- suppressWarnings(uniroot(f = dlogp, lower = 0, upper = rbound,
                                       k = k, theta = theta, 
                      lambda = lambda[m])$root)
    } else if(deriv0 < 0){
      lbound <- (deriv0/Ilamb(k, lambda[m]))-eps
      xmax <- suppressWarnings(uniroot(f = dlogp, lower = lbound, 
                                       upper = 0, k = k, theta = theta, 
                                       lambda = lambda[m])$root)
    } else{
      xmax <- 0
    }
    
    # 2: Determine grid bounds with grid-grower
    logfmax <- logp(xmax, k = k, theta = theta, lambda = lambda[m])
    sdmax <- sqrt(((-1) * d2logp(xmax, k = k, theta = theta, 
                                   lambda = lambda[m]))^(-1))
    b <- 0.01
    logb <- log(b)
    diff <- 2 * abs(logb)
    j <- 0
    while(diff >= logb){# left grid-grower
      xleft <- xmax - (2^j) * sdmax
      diff <- logp(xleft, k = k, theta = theta, lambda = lambda[m])-logfmax
      j <- j+1
    }
    diff <- 2 * abs(logb)
    j <- 0
    while(diff >= logb){# right grid-grower
      xright <- xmax + (2^j) * sdmax
      diff <- logp(xright, k = k, theta = theta, lambda = lambda[m])-logfmax
      j <- j+1
    }
    
    # 3: Construct grid and sample with inverse cdf method
    xgrid <- seq(xleft, xright, length = 100)
    dxgrid <- xgrid[2] - xgrid[1]
    logpxgrid <- logp(x = xgrid, k = k, theta = theta, lambda = lambda[m])
    logpxgrid <- logpxgrid - max(logpxgrid)
    cnorm <- 1/sum(exp(logpxgrid)*dxgrid)
    densval <- cnorm * exp(logpxgrid)
    cdf <- cumsum(densval* dxgrid)
    u <- runif(n = 1, min = 0, max = 1)
    while(u<cdf[1]){
      u <- runif(n = 1, min = 0, max = 1)
    }
    theta_draw <- xgrid[sum(cdf<=u)]
    theta[k] <- theta_draw
  }
  
  thetamat[m, ] <- theta
  if(isTRUE(progressbar)) utils::setTxtProgressBar(progbar, m)
}

if(burn == "half"){
  burnin <- niter * 0.5
  thetachain <- thetamat[-(1:burnin),]
  lambdachain <- lambda[-(1:burnin)]
} else{
  burnin <- burn
  thetachain <- thetamat[-(1:burn),]
  lambdachain <- lambda[-(1:burn)]
}

outlist <- list(niter = niter, K = K, thetachain = thetachain, 
                lambdachain = lambdachain)
return(outlist)
}





