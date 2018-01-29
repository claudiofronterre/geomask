#' Function needed for internal usage
fcut <- function(range, phi, kappa, r) {
  matern(range, phi = phi, kappa = kappa) - r
}
