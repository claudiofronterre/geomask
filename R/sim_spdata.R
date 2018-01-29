#' Simulation of gaussian spatial data using Cholesky decomposition
#'
#' Given a set of spatial coordinates this function generates gaussian distributed data whose spatial structure is defined by the user through the specification of a spatial correlation function.
#'
#' @param locations A two column matrix containing the locations (coordinates x
#'   and y).
#' @param cov_model A string indicating the type of the correlation function.
#'  For the available choices see \code{\link[geoR]{cov.spatial}}.
#' @param cov_pars A numeric vector with 2 elements with the covariance
#' parameters. The first element corresponds to the variance parameter
#'  \eqn{\sigma^2}. The second element or corresponds to the range parameter
#'   \eqn{\phi} of the correlation function.
#' @param nugget Value of the nugget parameter \eqn{\tau^2}.
#' @param kappa Numerical value for the additional smoothness parameter of
#'  the correlation function. Only required by the following correlation
#'  functions: "matern", "powered.exponential", "cauchy", "gencauchy" and
#'  "gneiting.matern".
#'
#' @return A numeric vector containing the simualted spatial data.
#'
#'
#' @examples
sim_spdata <- function(locations, cov_model = "matern", cov_pars = c(1, 0.16), nugget = 0.5, kappa = 0.5) {

  npts <- nrow(locations)
  kappa <- kappa
  nugget <- nugget

  #Prepare cholesky matrix to simulate data
  d_true <- as.matrix(dist(locs.true))
  sigma <- diag(nugget, nrow = npts) +
    cov.spatial(d_true, cov.pars = cov_pars, cov.model = cov_model, kappa = kappa)
  L <- t(chol(sigma))
  Z <- matrix(rnorm(npts), nrow = npts, ncol = 1)
  data <- L %*% Z
  return(as.numeric(data))
}
