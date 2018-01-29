#' Variogram fit adjusting for positional error
#'
#' This function fit the parameters of a variogram through N-weighetd least
#' squares taking also into account the presence of positional error.
#' @param data A numeric vector of spatial data.
#' @param locations A two column matrix containing the locations (coordinates x
#'   and y).
#' @param displacement The type of geomasking to be applied: either "gaussian"
#'   or "uniform".
#' @param delta A number that specify the standard deviation of the positional
#'   error in the case of Gaussian geomasking or the maximum displacement
#'   distance in the case of Uniform geomasking.
#' @param kappa Numerical value for the additional smoothness parameter of the
#'   matern correlation function.
#' @param ini Initial values for the parameters to be passed to the optimisation
#'   algorithm.
#' @param ... Control argumenets for the optimiser.
#'
#' @return A numeric vector containg the estimated parameters.
#'
#' @examples
varioadj <- function(data, locations, displacement = "gaussian", delta, kappa, ini, ...) {

  if(displacement == "uniform") delta <- delta/sqrt(6)

  #Adjusting for positional error (Matern correlation)
  integrand <- function(u.star, u, phi, delta) {
    exp(drice(u.star, sigma = sqrt(2)*delta, vee = u, log = T) +
          log(matern(u.star, phi = phi, kappa = kappa)))
  }

  rho <- function(u, phi, delta) {
    integrate(function(x) integrand(x, u, phi, delta),
              lower=0, upper = Inf)$value
  }

  rho <- Vectorize(rho)

  v.loss.adj <- function(parms, u, v, n) {
    sigmasq <- exp(parms[1])
    phi <- exp(parms[2])
    tausq <- exp(parms[3])
    v.mod <- (tausq + sigmasq) - sigmasq*rho(u, phi = phi, delta = delta)
    loss <- sum(n*(v - v.mod)^2)
    return(loss)
  }

  #Calculates the sample variogram
  df <- data.frame(x1 = locations[,1], x2 = locations[,2], y = data)
  coordinates(df) <- ~ x1 + x2
  vario <- variogram(y ~ 1, df, cutoff = 3)

  #Optimisation
  res <- optim(par = log(ini), fn = v.loss.adj,
               u = vario$dist, v = vario$gamma, n = vario$np)
  params <- c(exp(res$par[1]), exp(res$par[2]), exp(res$par[3]))
  names(params) <- c("sigma2", "phi", "nugget")
  return(params)
}
