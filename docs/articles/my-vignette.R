## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 5
)
library(geomask)

## ---- fig.align='center'-------------------------------------------------
n <- 10000
max_dist <- 2
delta <- max_dist/sqrt(6) 

# Generate the original points
true_locations <- matrix(0, nrow = n, ncol = 2)

# DHS displacement (Uniform geomasking)
unif_locs <- geomasking(locations = true_locations, displacement = "uniform", delta = max_dist)

# Gaussian geomasking
gauss_locs <- geomasking(locations = true_locations, displacement = "gaussian", delta = delta)

# Plot the displaced locations
par(mfrow = c(1, 2))
plot(unif_locs, col = rgb(50, 205, 50, alpha = 35, maxColorValue = 255), xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), xlab = "", ylab = "", main = "Unifrom")
plot(gauss_locs, col = rgb(50, 205, 50, alpha = 35, maxColorValue = 255), xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), xlab = "", ylab = "", main = "Gaussian")

## ---- fig.align='center'-------------------------------------------------
# Compare displaced distances
dist1 <- sqrt(unif_locs[, 1]^2 + unif_locs[, 2]^2)
dist2 <- sqrt(gauss_locs[, 1]^2 + gauss_locs[, 2]^2)
par(mfrow = c(1, 2))
hist(dist1, col = "lightgreen", main = "Uniform")
hist(dist2, col = "lightgreen", main = "Gaussian")

# Compare within distances of displaced points
d1 <- dist(gauss_locs)
d2 <- dist(unif_locs)
hist(d1, col = "lightgreen", xlim = c(0, 5), main = "Uniform")
hist(d2, col = "lightgreen", xlim = c(0, 5), main = "Gaussian")

