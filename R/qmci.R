#' Likelihood calculation with Quasi Monter Carlo Integration
#'
#' This function calculates the liklelihood of the geostatistical model in
#' presence of positional error using Quasi Monte Carlo Integration. It is used
#' internally in the optimisation routine.
qmci <- function(data, sequence, sigma2, phi, nugget, mu, delta, kappa) {
  ave <- apply(data, MARGIN = 1, FUN = function(x, delta) {
    xi.star <- qnorm(sequence[,1], mean = x[3], sd = delta)
    yi.star <- qnorm(sequence[,2], mean = x[4], sd = delta)
    xj.star <- qnorm(sequence[,3], mean = x[5], sd = delta)
    yj.star <- qnorm(sequence[,4], mean = x[6], sd = delta)
    u.star <- sqrt((xi.star - xj.star)^2 + (yi.star - yj.star)^2)
    sigma <- sigma2 + nugget
    rho <- sigma2*matern(u.star, phi = phi, kappa = kappa)/sigma
    z <- ((x[1] - mu)^2 + (x[2] - mu)^2 - 2*rho*(x[1] - mu)*(x[2] - mu))/sigma
    den <- 2*pi*sigma*sqrt(1 - rho^2)
    mean((1/den)*exp(-z/(2*(1 - rho^2))), na.rm = T)
  }, delta = delta)
  return(-sum(log(ave)))
}
