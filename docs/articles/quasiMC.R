## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
library(randtoolbox)
library(geoR)
library(VGAM)

#Parameters of the model
phi <- 0.25; kappa <- 0.5; sigma2 <- 1; nugget <- 0.2; delta <- phi*0.3; dij <- 10; yij <- 1

n <- 2

#First solution
halt1 <- halton(n)
u.star <- qrice(halt1, sigma = sqrt(2)*delta, vee = dij)
mean(dnorm(yij, mean = 0, 
           sd = sqrt(2*(nugget + sigma2*(1 - matern(u.star, phi = phi, kappa = kappa))))))

#Second solution
xi <- 1; yi <- 2; xj <- 11; yj <- 2
halt2 <- halton(n, dim = 4)
xi.star <- qnorm(halt2[,1], mean = xi, sd = delta)
xj.star <- qnorm(halt2[,2], mean = xj, sd = delta)
yi.star <- qnorm(halt2[,3], mean = yi, sd = delta)
yj.star <- qnorm(halt2[,4], mean = yj, sd = delta)
u.star <- sqrt((xi.star - xj.star)^2 + (yi.star - yj.star)^2)
mean(dnorm(yij, mean = 0, 
          sd = sqrt(2*(nugget + sigma2*(1 - matern(u.star, phi = phi, kappa = kappa))))))


