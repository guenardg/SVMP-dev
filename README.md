# R Package SVMP-dev: Spatial Variance Modeling and Prediction -- Development

Documents herein contains the development files and R implementations of a
framework for spatial variance modeling and prediction. It involves a set of
distance-based weighting functions that are used to calculate spatial
covariances from sets of pairwise distance between locations. These spatial
covariances operators are then used to obtain an eigensystem whose eigenvectors
represent a suite of spatial variation patterns (each eigenvalue represents the
size of the structure represented by its associated eigenvectors).

The method is flexible. The space can simply be defined in terms of a Euclidean
space, such as a series (one-dimensional; e.g., a transect) or a cartesian plane
(two-dimensional; e.g., a plot), or any other relevant manifolds, for instances
geodesics using great circle distances on spheroids or ellipsoids as a distance
metric, or curves surface using shortest path lengths as a metric.

The R language package called SVMP that is included here is for research
purposes and should be expected to change over time. It is an unpublished and
unrevised work that has been put here for informational, public commenting, and
revision purposes. Therefore, it must not be relied upon for application until
the method and implementation it be duly revised and published.

## Authored by

Guillaume Guénard

Université de Montréal, June 2019 -- Feb. 2020
