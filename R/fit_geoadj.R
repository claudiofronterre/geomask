#' Fit a geostatistical model to geomasked data
#'
#' This function fit a geostatistical model using composite likelihood to
#' spatial data that have positional error due to geoamsking.
#'
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
#' @param n_sequence A numeric value. It will define the lenght of the halton
#'   sequence for the quasi monte carlo integration. A longer sequence requires
#'   more computational time but provides more accurate results. Defaul to 10.
#' @param approx If TRUE (defautl is FALSE) it will use an approximation to
#'   calculate the composite likelihood. If set to TRUE a threshold value in the
#'   argument thresh needs to be provided.
#' @param thresh If approx is TRUE this defines the level of the approximation.
#'   By default is 0.000005 and garuantess a good comprise between speed and
#'   accuracy. Bigger values will make the computation faster but less accurate.
#' @param ini Initial values for the parameters to be passed to the optimisation
#'   algorithm.
#' @param method The optimisation method to be used. Default is "BFGS".
#'
#' @return A list containing the set of estimated parameters. The likelihood
#'   evaluated at the estimated parameters and a code to asses convergence of
#'   the algorithm.
#' @export
#'
#' @examples
fit_geoadj <- function(data, locations, displacement = "gaussian", delta, kappa, n_sequence = 10, approx = F, thresh = 0.000005, ini, method = "BFGS") {

  if(displacement == "uniform") delta <- delta/sqrt(6)
  npts <- nrow(locations)
  locs.obs <- locations
  r <- thresh

  # Prepare data
  id <- t(combn(npts, 2))
  res <- t(apply(id, MARGIN = 1, FUN = function(x) {
    cbind(yi =data[x[1]] , yj = data[x[2]],  x1 = locs.obs[x[1], 1],
          y1 = locs.obs[x[1], 2],  x2 = locs.obs[x[2], 1],  y2 = locs.obs[x[2], 2])
  }))
  dij <- sapply(1:nrow(id), FUN = function(x) dist(rbind(locs.obs[id[x,1],], locs.obs[id[x,2], ])))
  res <- cbind(res, dij)
  colnames(res) <- c("yi", "yj", "x1", "y1", "x2", "y2", "dij")

  halt2 <- halton(n_sequence, dim = 4)

  if (approx == F) {
    f.loss <- function(params, res, delta){
      sigma2 <- exp(params[1])
      phi <- exp(params[2])
      nugget <- exp(params[3])
      loss <- qmci(data = res, sequence = halt2, sigma2 = sigma2, phi = phi, nugget = nugget, mu = 0,
                   delta = delta, kappa = kappa)
    }
    CL <- optim(par = log(ini), method = method, fn = f.loss, res = res, delta = delta)
  } else {
    f.loss <- function(params, res, r, delta){
      sigma2 <- exp(params[1])
      phi <- exp(params[2])
      nugget <- exp(params[3])

      cutoff <- uniroot(fcut, interval = c(0, 50 * phi + 1),
                        phi = phi, kappa = kappa, r = r, tol = 0.000001)$root
      res2 <- res[dij > cutoff, ]
      res <- res[dij < cutoff, ]

      l1 <- qmci(data = res, sequence = halt2, sigma2 = sigma2, phi = phi, nugget = nugget, mu = 0,
                 delta = delta, kappa = kappa)
      sigma <- sigma2 + nugget
      z <- ((res2[, "yi"] - 0)^2 + (res2[, "yj"] - 0)^2)/sigma
      c <- (1/(2*pi*sigma))*exp(-z/2)
      l2 <- - sum(log(c))
      loss <- l1 + l2
    }
    CL <- optim(par = log(ini), method = method, fn = f.loss, res = res, delta = delta, r = r)
  }

  params = exp(CL$par)
  names(params) <- c("sigma2", "phi", "nugget")
  results <- list(params = params, lik = CL$value, convergence = CL$convergence)
  return(results)
}
