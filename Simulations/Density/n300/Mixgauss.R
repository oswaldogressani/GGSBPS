#------------------------------------------------------------------------------#
#           Routine to simulate samples from a mixture of 3 Gaussians          #
#               Copyright Oswaldo Gressani. All rights reserved.               #
#------------------------------------------------------------------------------#


mixgauss <-function(n, alpha = c(0.25,0.50,0.25), mu = c(0.1,0.50,0.9),
                    range = c(0,1)){
  
  # Gaussian sd
  Nsd <- c(0.03, 0.06, 0.03)
  
  # Sampling with replacement with weights alpha.
  w <- sample(seq_len(3), n, replace = TRUE, prob = c(alpha[1], alpha[2], alpha[3]))
  
  # Sample of size n from the mixture of 3 Gaussian densities
  sample <- (w == 1) * rnorm(n, mu[1], Nsd[1]) +
              (w == 2) * rnorm(n, mu[2], Nsd[2]) + 
                (w == 3) * rnorm(n, mu[3], Nsd[3])
  
  # Keep sampled values in range
  outbounds <- which(sample < range[1] | sample > range[2])
  outboundslen <- length(outbounds)
  if(outboundslen>0){
    for(j in 1:outboundslen){
      jidx <- outbounds[j]
      sample_temp <- (w[jidx] == 1) * rnorm(1, mu[1], Nsd[1]) + 
                      (w[jidx] == 2) * rnorm(1, mu[2], Nsd[2]) + 
                        (w[jidx] == 3) * rnorm(1, mu[3], Nsd[3])
      while(sample_temp < 0 | sample_temp > 1){
        sample_temp <- (w[jidx] == 1) * rnorm(1, mu[1], Nsd[1]) + 
                        (w[jidx] == 2) * rnorm(1, mu[2], Nsd[2]) + 
                          (w[jidx] == 3) * rnorm(1, mu[3], Nsd[3])
      }
      sample[jidx] <- sample_temp
    }
  }
  btwrange <- as.logical(1-(any(sample < range[1]) | any(sample > range[2])))
  
  densfun <- function(x){
    val <- alpha[1] * dnorm(x = x, mean = mu[1], sd = Nsd[1]) +
             alpha[2] * dnorm(x = x, mean = mu[2], sd = Nsd[2]) +
               alpha[3] * dnorm(x = x, mean = mu[3], sd = Nsd[3])
    return(val)
  }
  x <- seq(range[1], range[2], length = 500)
  denseval <- sapply(x, densfun)
  
  outlist <- list(n = n, alpha = alpha, sample = sample,
                  densfun = densfun, x = x, fx = denseval,
                  range = range,
                  btwrange = btwrange)
  
  return(outlist)
}




