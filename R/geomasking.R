#' Apply Gaussian or Uniform geomasking to a set of spatial coordinates.
#'
#' @param locations A two column matrix containing the locations (coordinates x
#'   and y) to which apply geomasking.
#' @param displacement The type of geomasking to be applied: either "gaussian"
#'   or "uniform".
#' @param delta A number that specify the standard deviation of the positional
#'   error in the case of Gaussian geomasking or the maximum displacement
#'   distance in the case of Uniform geomasking.
#' @return A matrix with the displaced locations.
#' @examples
geomasking <- function(locations, displacement, delta) {
  if(ncol(locations) < 2) warning(paste("The number of columns of the", class(locations), "object provided must be equal to 2 (coordinates x and y)"))
  n <- nrow(locations)
  if (displacement == "gaussian") {
    new_locations <- locations + matrix(rnorm(2 * n, mean = 0, sd = delta), ncol = 2)
  } else {
    rd <- runif(n, min = 0, max = delta)
    theta <- runif(n, min = 0, max = 2 * pi)
    new_locations <- cbind(locations[,1] + rd*cos(theta), locations[,2] + rd*sin(theta))
  }
  return(new_locations)
}
