#------------------------------------------------------------------------------#
#                     Analysis of Stamps data                                  #
#        PART I of this algorithm is code adapted from the JOPS package        #
#------------------------------------------------------------------------------#

library("ggplot2")
library("Rcpp")
library("JOPS")
library("locfit")

rm(list = ls())
sourceCpp("KercubicBspline.cpp")
source("GGSBPS.R")

#----- PART I Application to Stamps data
data(stamp)
raw = rep(stamp$thick, times = stamp$count)
set.seed(123)
m = length(raw)
raw = 1000 * sample(raw, m / 3)

xl = 55
xr = 135
brk <- seq(xl, xr, by = 1)
h = hist(raw, breaks = brk, plot = F)
x = h$mids
y = h$counts
m = length(y)
h <- x[2] - x[1]

# B-spline basis matrices, for fitting and plotting
nseg <- 27
B <- bbase(x, xl, xr, nseg = nseg)
xg <- seq(xl, xr, by = 0.1)
Bg <- bbase(xg, xl, xr, nseg = nseg)
n <- ncol(B)

# Penalty matrix and parameters for gamma prior
dk = (xr - xl) / nseg;
per = 10;
phi = 2 * pi * dk / per;
dp = c(1, -2 * cos(phi), 1)
Dp = matrix(0, n - 2, n)
for (j in 1:(n - 2)) Dp[j, j:(j + 2)] = dp;
D = diff(Dp)
P = t(D) %*% D
D2 = diff(diag(n), diff = 2)
P2 = t(D2) %*% D2

# Fit the model
llas <- seq(-1, 2, by = 0.05)
laps <- LAPS_dens(B, P, y, llas, mon = F)
fhat <- exp(Bg %*% laps$alpha)
laps2 <- LAPS_dens(B, P2, y, llas, mon = F)
fhat2 <- exp(Bg %*% laps2$alpha)

# Make plots
DF = data.frame(x = x, y = y)
DF1 = data.frame(x = xg, fhat = fhat, fhat2 = fhat2)
DF2 = data.frame(x = llas, w = laps$weights, w2 = laps2$weights)

plt1 = ggplot(data = DF) +
  geom_bar(aes(x = x, y = y), stat = "identity",
           fill= "wheat3", col = "white", size = 0.1) +
  geom_hline(yintercept = 0) +
  xlab("Thickness (micrometer)") + ylab("Frequency") +
  ggtitle("Hidalgo stamps, harmonic penalty") +
  geom_line(data = DF1, aes(x = x, y = fhat2), col = "red", size = 0.8) +
  geom_line(data = DF1, aes(x = x, y = fhat), col = "blue", size = 0.8) +
  scale_x_continuous(breaks = seq(50, 140, by = 10)) +
  JOPS_theme()

titl =  bquote(lambda ~ "weights")
plt2 = ggplot(data = DF2) +
  geom_line(aes(x = x, y = w2), col = 'red') +
  geom_point(aes(x = x, y = w2), col = 'red', pch = 15) +
  geom_line(aes(x = x, y = w), col = 'blue') +
  geom_point(aes(x = x, y = w), col = 'blue') +
  xlab(expression(log10(lambda))) + ylab('Weight') +
  #ggtitle('lambda weights') +
  ggtitle(titl) +
  JOPS_theme()

KLps <- dim(B)[2]

#------------------------ PART II Griddy-Gibbs routine
set.seed(1234)
xl <- 55
xr <- 135
modfit <- GriddyGibbs(x = x, y = y, niter = 10000, xl = xl, xr = xr, K = KLps,
                      penorder = 2)
thetahat <- colMeans(modfit$thetachain)

xx <- seq(xl, xr, length = 200)
Bx <- Rcpp_KercubicBspline(xx, lower = xl, upper = xr, K = modfit$K)
muhat <- as.numeric(exp(Bx%*%thetahat))

# Compare psPoisson fit with Griddy-Gibbs fit
pdf(file = "Stamps.pdf", width = 12.3, height = 6)
par(mfrow = c(1,1))
plot(hist(raw, breaks = brk, plot = F), 
     main = "Hidalgo stamps", xlab = "Thickness (micrometer)",
     xlim = c(60,140))
lines(DF1$x, DF1$fhat2, type = "l", col = "brown", lwd = 2)
lines(xx, muhat, type = "l", col = "cornflowerblue", lwd = 2)
legend("top",col = c("cornflowerblue","brown"),
       c("GGSBPS","Laplacian-P-splines"), lty=1, lwd = 2,
       bty = "n")
dev.off()

# Lambda value with LPS
round(as.numeric(laps2$lambda),3)

# Lambda value with GGS-BPSP
round(mean(modfit$lambdachain),3)




