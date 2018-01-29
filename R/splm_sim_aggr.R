#' Spatial model for spatially aggregated data
#'
#' This function provide  model based parameters estiamtion for a geostatical model when spatial data are avialable at a corarser scale than their natural resolution. It is mainly for simulation purposing and testing.
#' @param data A numeric vector of spatial data.
#' @param mc_points Number of Monte Carlo points to calculate the correlation matrix.
#' @param tau.sq Numeric value for the nugger parameter \eqn{\nugget^2}.
#' @param ncovariates Number of covariates to simulate.
#' @param beta Values for the beta parameters.
#' @param aggr.cell To how many cells the points should be aggregated.
#' @param fix.nugget A logical value. Should the nugget be fixed?
#' @param fix.nug If fix.nugget = TRUE then a numeric value should be provided.
#' @param ini Initial parameters values for the optimisation algorithm.
#' @param message If FALSE suppress all the messages from the fitting algorithm.
#'
#' @return
#' @export
#'
#' @examples
splm_sim_aggr <- function(data,
                          mc_points = 100,
                          tau.sq = 0,
                          ncovariates = 2,
                          beta = c(3.5, -1.2),
                          aggr.cell = 9,
                          fix.nugget = F,
                          fix.nug = 0.5,
                          ini = c(0.5, 0.5),
                          message = F) {
  # gres = grid resolution
  # grlim = limits of the grid
  # covmodel = true correlation function
  # kappa = true kappa
  # ncovariates = number of simulated covariates
  #
  #

  S <- data

  # Simulate covariates
  D <- matrix(data = rnorm(n*ncovariates), nrow = n, ncol = ncovariates)

  # Beta parameters
  beta <-  beta <- as.matrix(beta)

  # Nugget
  tau <- tau.sq

  #Generate the data
  Y <- rnorm(n, mean = D%*%beta + S, sd = sqrt(tau))

  # Convert to raster
  xyz = data.frame(x = coords[,1], y = coords[,2], z = Y, D)
  rs <- rasterFromXYZ(xyz)

  # Aggregate everything to a coarse resolution
  rs.coarse <- aggregate(rs, fact = ceiling(sqrt(ncell(rs)/aggr.cell)))

  # Convert back to points
  df <- rasterToPoints(rs.coarse)


  #Run the model
  nblocks <- dim(rs.coarse$z)[1]^2
  xy <- spsample(as(extent(rs.coarse$z), "SpatialPolygons"), n = mc_points*nblocks, type = "regular")
  pp <- data.frame(x = coordinates(xy)[,1], y = coordinates(xy)[ , 2], block = cellFromXY(rs.coarse$z, xy))

  mc_corr <- function(i, j, phi) {
    p1 <- pp[pp[,"block"] == i, c(1,2)]
    p2 <- pp[pp[,"block"] == j, c(1,2)]
    blockdist <- dist2(p1, p2)[upper.tri((dist2(p1, p2)), diag = T)]
    sum(matern(blockdist, phi = phi, kappa))/(nrow(p1)*nrow(p2))
  }


  Dm <- df[,c(4,5)]
  Y <- as.matrix(df[,3])

  f.loss <- function(pars, ff, y) {
    nu <- exp(pars[1])
    phi <- exp(pars[2])

    comb <- t(combn(nblocks,2))
    res <- apply(comb, 1, function(x) mc_corr(i = x[1], j = x[2], phi = phi))
    corr <- diag(nblocks)
    corr[upper.tri(corr)] <- res
    corr[lower.tri(corr)] <- res
    corr

    diag(corr) <- 1 + nu

    corr <- chol(corr)
    ldet.5 <- sum(log(diag(corr)))
    corr <- chol2inv(corr)

    beta_hat <- solve(crossprod(ff, corr)%*%ff,
                      crossprod(ff, corr)%*%y)

    z <- y - ff%*%beta_hat
    sigma2_hat <- mean(crossprod(corr,z)*z)

    log_lik <- ldet.5 + nblocks*(1 + log(2*pi*sigma2_hat))/2

    attr(log_lik, "param") <- c(beta = beta_hat, sigma2 = sigma2_hat, phi = phi, nugget = nu*sigma2_hat)
    return(log_lik)
  }


  res <- optim(par = log(ini), fn = f.loss, ff = cbind(1, Dm), y = Y)
  lkhat <- attr(f.loss(res$par, ff = cbind(1, Dm), y = Y), "param")
  ls <- list(lkhat, xyz, df)
  names(ls) <- c("params", "original_data", "aggregated_data")
  return(ls)
}

