#------------------------------------------------------------------------------#
#               Application of GGS-BPSP to Geyser data                         #
#       PART I of this algorithm is code adapted from the JOPS package         #
#------------------------------------------------------------------------------#

library("ggplot2")
library("Rcpp")
library("JOPS")

rm(list = ls())
sourceCpp("KercubicBspline.cpp")
source("GGSBPS.R")

#----- PART I Application to Old faithful Geyser data from JOPS package
data(faithful)
v = faithful[, 1]  # Eruption length

# Histogram with narrow bin width
bw = 0.05
hst = hist(v, breaks = seq(1, 6, by = bw), plot = F)
x = hst$mids
y = hst$counts
Data = data.frame(x = x, y = y)

# Poisson smoothing
lambdas = 10^seq(-3, -0, by = 0.1)
aics = NULL
for (lambda in lambdas) {
  fit = psPoisson(x, y, nseg = 17, pord = 2, lambda = lambda, show = F)
  aics = c(aics, fit$aic)
}

# Finding min AIC
ka = which.min(aics)
lambda = lambdas[ka]

# Gridded data for plotting using opt tuning
fit = psPoisson(x, y, nseg = 17, pord = 2, lambda = lambda, show = F)
F1 = data.frame(x = log10(lambdas), y = aics)
Fit = data.frame(x = fit$xgrid, y = fit$mugrid)
plt1 = ggplot(aes(x = x, y = y), data = F1) +
  geom_point() +  ggtitle(paste("Bin width", bw)) +
  xlab(expression(log10(lambda))) + ylab('AIC') +
  JOPS_theme()

plt2 = ggplot(aes(x = x, y = y,  fill = I("wheat3")), data = Data) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  xlab("Eruption length (min)") + ylab("Frequency") +
  ggtitle("Old Faithtful geyser") +
  geom_line(data = Fit, col = I("steelblue"), size = 1.0) +
  JOPS_theme()

# Number of B-spline basis functions in the penalized likelihood approach
KpL <- dim(fit$B)[2]

#------------------------ PART II Griddy-Gibbs routine
set.seed(1234)
modfit <- GriddyGibbs(x = x, y = y, niter = 10000, penorder = 2, K = KpL)
thetahat <- colMeans(modfit$thetachain)

xx <- seq(min(x), max(x), length = 200)
Bx <- Rcpp_KercubicBspline(xx, lower = min(x), upper = max(x), K = modfit$K)
muhat <- as.numeric(exp(Bx%*%thetahat))

# Compare psPoisson fit with Griddy-Gibbs fit
pdf(file = "Geyser.pdf", width = 12.3, height = 6)
plot(hst, main = "Old Faithful Geyser", xlab = "Eruption length (min)" )
lines(Fit$x, Fit$y, type = "l", col = "brown", lwd = 2)
lines(xx, muhat, type = "l", col = "cornflowerblue", lwd = 2)
legend("top",col = c("cornflowerblue","brown"),
       c("GGSBPS","Penalized likelihood"), lty=1, lwd = 2,
       bty = "n", cex = 0.90)
dev.off()

# Lambda value with AIC from penalized likelihood method
round(lambda,3)

# Lambda value with GGS-BPSP
round(mean(modfit$lambdachain),3)



