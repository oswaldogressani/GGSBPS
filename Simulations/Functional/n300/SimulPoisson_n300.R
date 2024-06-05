#------------------------------------------------------------------------------#
#       Simulation study for GGSBPS algorithm (Poisson data)                   #
#          Copyright Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#               

#-- Load libraries and source code
library("Rcpp")
library("ggplot2")
rm(list = ls())
sourceCpp("KercubicBspline.cpp")
source("GGSBPS.R")

S <- 100
n <- 300

set.seed(2024)

#-- Hosting matrices/vectors
xl <- (-0.1)
xr <- (1.1)
xp <- seq(xl, xr, by = 0.05)  # coarse grid for performance measurement
xf <- seq(xl, xr, by = 0.001) # fine grid for plot
fxp <- matrix(0, nrow = S, ncol = length(xp))

fxf <- matrix(0, nrow = S, ncol = length(xf))
penorder <- 2

# Target function to be estimated
mu <- function(x, j = 3){
  val <- sqrt(x * (1 - x) * 1.7) * 6 * cos(2 * pi * (1 + 2 ^ ((9 - 4 * j) / 5)) / 
                                             (x + 2 ^ ((9 - 4 * j) / 5))) + 5
  return(val)
}

# Progress bar showing simulation status
cat(paste0("GGSBPS simulation running for ",S, " iterations \n"))
progbar <- utils::txtProgressBar(min = 1, max = S, initial = 1,
                                 style = 3, char =">")

tic <- proc.time()

for(s in 1:S){

#-- Generate data
x <- runif(n = n, min = 0, max = 1)     # Random grid of length n
y <- rpois(n = n, lambda = sapply(x, mu)) # Generate Poisson data

#-- Fit model with GGSBPS algorithm
GGSfit <- GriddyGibbs(x = x, y = y, niter = 1000, penorder = penorder,
                      xl = xl, xr = xr, progressbar = FALSE)
Bxf <- Rcpp_KercubicBspline(x = xf, lower = xl, upper = xr, K = GGSfit$K)
thetahat <- colMeans(GGSfit$thetachain)
fhat <- as.numeric(exp(Bxf%*%thetahat))
Bxp <- Rcpp_KercubicBspline(x = xp, lower = xl, upper = xr, K = GGSfit$K)
fhatxp <- as.numeric(exp(Bxp%*%thetahat))

#-- Record fit (point estimate)
fxf[s, ] <- fhat
fxp[s, ] <- fhatxp

utils::setTxtProgressBar(progbar, s)
}

toc <- proc.time() - tic

#-- Computation of performance measures
ftrue <- mu(xp)
Bias <- colMeans(fxp - matrix(rep(ftrue, S), nrow = S, byrow = TRUE))
fbar <- colMeans(fxp)
ESE <- sqrt((1/(S-1)) * 
              colSums((fxp - matrix(rep(fbar, S), nrow = S, byrow = TRUE))^2))
RMSE <- sqrt(colMeans((fxp - matrix(rep(ftrue, S), nrow = S, byrow = TRUE))^2))

perf <- matrix(0, nrow = length(xp), ncol = 4)
colnames(perf) <- c("x", "Bias", "ESE", "RMSE")
perf[,1] <- xp
perf[,2] <- Bias
perf[,3] <- ESE
perf[,4] <- RMSE

round(perf, 3)

subperf <- round(perf[c(5,7,9,11,13,15,17,19,21),],3)
write.csv(file = "ScenarioD.csv", subperf, row.names = FALSE)


#-- Plot with ggplot2
xdom <- seq(0,1, length = 500)
fxdom <- sapply(xdom, mu)

datatar <- data.frame(x = xdom, y = fxdom)

medcol <- grDevices::rgb(240, 101, 41, maxColorValue = 255)
colors <- c("Target" = "black", "GGSBPS" = "gray65", "Median" = medcol)
linetypes <- c("Target" = 1, "GGSBPS" = 1, "Median" = 4)

tarfun <- ggplot(data = datatar, aes(x = x, y = y)) +
  labs(x = expression(x), y = expression(mu(x)), color = "Legend",
       linetype = "Legend")  +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) + 
  theme_classic() + xlim(0,1)

for(s in 1:S){
  datafxf <- data.frame(x = xf, y = fxf[s,])
  tarfun <- tarfun + geom_line(data = datafxf, aes(y=y, color = "GGSBPS", 
                                                   linetype = "GGSBPS"))
}
tarfun <- tarfun +  geom_line(data = datatar, 
                              aes(y=y, color = "Target", linetype = "Target"),
                              linewidth = 1.1) 
datamedian <- data.frame(x = xf, y = apply(fxf, 2, "median"))
tarfun <- tarfun + geom_line(data = datamedian, aes(y=y, color = "Median", 
                                                    linetype = "Median"),
                             linewidth = 1.1) +
  theme(legend.position = "top",
        legend.title = ggplot2::element_blank()) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) 

png(file = "ScenarioD.png", width = 1000, height = 500)
tarfun
dev.off()

pdf(file = "ScenarioD.pdf", width = 12.3, height = 6)
tarfun
dev.off()

subperf

toc



















