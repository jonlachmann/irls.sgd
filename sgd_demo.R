# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-04-23

sgd.demo.data <- function (nobs, means=c(0.5,2)) {
  x1 <- rnorm(nobs, means[1],2)
  x2 <- rnorm(nobs, means[2],2)
  return(matrix(c(x1,x2), ncol=2))
}

sgd.demo.fun <- function (x, data) {
  sum(t(data) * c(sin(x[1]), cos(x[2]))) #t(data[,1, drop=F])%*%sin(x[1]) + t(data[,2, drop=F])*cos(x[2])
}

sgd.demo.grad <- function (x, data) {
  rowMeans(t(data) * c(cos(x[1]), -sin(x[2]))) #c(t(data[,1, drop=F])%*%cos(x[1]), -t(data[,2, drop=F])%*%sin(x[2]))
}

sgd.demo.grid <- function (data, fun=sgd.demo.fun, lims=c(0.9*pi,pi*2.1,0,pi*2), resolution=pi/20) {
  xgrid <- seq(lims[1], lims[2], resolution)
  ygrid <- seq(lims[3], lims[4], resolution)
  fvals <- matrix(NA, length(xgrid), length(ygrid))
  for (i in 1:length(xgrid)) {
    for (j in 1:length(ygrid)) {
      x <- c(xgrid[i], ygrid[j])
      fvals[i,j] <- fun(x, data)
    }
  }
  return(list(x=xgrid, y=ygrid, z=fvals))
}

set.seed(1234)
data <- sgd.demo.data(1000)
grid <- sgd.demo.grid(data)
gd_res <- sgd(sgd.demo.grad, data, start=c(pi*1.5-1.5,pi-2.8), alpha=0.25, decay=0.999, maxit=50)
sgd_res <- sgd(sgd.demo.grad, data, subs=0.001, start=c(pi*1.5-1.5,pi+2.8), alpha=0.1, decay=0.999, maxit=300)
bsgd_res <- sgd(sgd.demo.grad, data, subs=0.05, start=c(pi*1.5+1.5,pi+2.8), alpha=0.1, decay=0.999, maxit=100)

{
par(mfrow=c(1,1))
contour(grid, main="Gradient Descent Variations")


  lines(gd_res$xhist[,1], gd_res$xhist[,2], col="blue")
points(gd_res$xhist[,1], gd_res$xhist[,2], col="blue", pch=15)

  lines(sgd_res$xhist[,1], sgd_res$xhist[,2], col="red")
points(sgd_res$xhist[,1], sgd_res$xhist[,2], col="red", pch=16)

  lines(bsgd_res$xhist[,1], bsgd_res$xhist[,2], col="darkgreen")
points(bsgd_res$xhist[,1], bsgd_res$xhist[,2], col="darkgreen", pch=17)

legend(5.5,1, legend=c("Gradient Descent", "Stochastic GD", "Batch Stochastic GD"),col=c("blue", "red", "darkgreen"), pch=c(15,16,17))
}

filled.contour(grid, main="Gradient Descent Variations")