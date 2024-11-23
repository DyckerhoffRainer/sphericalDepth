#' Angular halfspace depth
#'
#' This function computes the angular halfspace depth (it is not assumed that the points are in general position).
#'
#'@param x a matrix of format n x d that contains the data points with respect to which the depth has to be computed
#'@param mass a vector of length n that contains either the probabilities or the (integer) multiplicities of the n data points
#'@param z a matrix of format m x d containing the points whose depth should be computed
#'@param ind a vector of length l containing the indices of the data points x whose depth should be computed
#'@param target the dimension to which the data is projected, possible values are 1, 2, 3
#'@param method string, possible values are "standard", "adaptive", "single", "multiple", "nGP", "GP".
#'Denotes the variant of the algorithm that is used to compute the signed halfspace depth in the target dimension.
#'
#'For target=2 only "standard", "adaptive", "single", "multiple" are valid. For target=2 "standard" is the same as "adaptive", meaning that depending on the number of points for which the depth has to be computed, either "single" or "multiple" is selected.
#'
#'For target=3 only "standard", "nGP", "GP" are valid. For target=3 "standard" is the same as "nGP", meaning that the points do not have to be in general position. If method = "GP", then is is assumed that all points are in general position, so that a faster algorithm can be used.
#'
#'@param nThreads the number of threads to use in a multiprocessor environment. For the default 0 the number of threads is automatically chosen depending on the used hardware.
#'@return a vector of length m + l containing the angular halfspace depths of the points z and the data points x whose indices are given in ind.
#'If the argument 'mass' is missing or if an integer vector is supplied for 'mass', then the depths are given as integer values, i.e., multiplied by n.
#'If a double vector is supplied for the argument 'mass', then the depths are given as doubles.
#'
#'@details
#'The routine uses the recursive algorithm of Dyckerhoff, Nagy (2024). Unless method='GP', it does not assume that the points are in general position. The points in z may coincide with points in x. z or ind may be missing (but not both).
#'
#'Regarding the choice of parameters 'target' and 'method' the default values usually will give the best results.
#'In most cases this routine will be slower than 'ahD_Comb'.

#'@references Dyckerhoff, R., Nagy, S. (2024). Exact computation of angular halfspace depth.
#'
#'@author Rainer Dyckerhoff
#'
#'@export
#'@useDynLib sphericalDepth, .registration=TRUE
#'
#'@examples
#'d <- 4
#'n <- 100
#'m <- 100
#'l <- 100
#'# simulate data uniformly distributed on the sphere
#'x <- matrix(rnorm(n*d),nrow=n)
#'x <- x / sqrt(rowSums(x*x))
#'z <- matrix(rnorm(m*d),nrow=m)
#'z <- z / sqrt(rowSums(z*z))
#'# vector of probabilities
#'prob <- rep(1/n, n)
#'# vector of point multiplicities
#'count <- as.integer(rep(1,n))
#'
#'res <- matrix(nrow=20, ncol=m+l)
#'
#'# combinatorial algorithm, w/o argument 'mass', different target dimensions and methods
#'system.time(res[ 1,] <- ahD_Rec(x, z = z, ind = 1:l, target = 1))
#'system.time(res[ 2,] <- ahD_Rec(x, z = z, ind = 1:l, target = 2, method = "single"))
#'system.time(res[ 3,] <- ahD_Rec(x, z = z, ind = 1:l, target = 2, method = "multiple"))
#'system.time(res[ 4,] <- ahD_Rec(x, z = z, ind = 1:l, target = 3, method = "nGP"))
#'system.time(res[ 5,] <- ahD_Rec(x, z = z, ind = 1:l, target = 3, method = "GP"))
#'
#'# combinatorial algorithm, pass a vector of probabilitiues to mass, different target dimensions
#'#and methods
#'system.time(res[ 6,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 1))
#'system.time(res[ 7,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 2, method = "single"))
#'system.time(res[ 8,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 2, method = "multiple"))
#'system.time(res[ 9,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 3, method = "nGP"))
#'system.time(res[10,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 3, method = "GP"))
#'# multiply by n since the depths are since we want to compare with the integer version of the depth
#'res[6:10,] <- res[6:10,] * n
#'
#'# combinatorial algorithm, pass a vector of multiplicieties to mass, different target dimensions
#'# and methods
#'system.time(res[11,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 1))
#'system.time(res[12,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 2, method = "single"))
#'system.time(res[13,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 2, method = "multiple"))
#'system.time(res[14,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 3, method = "nGP"))
#'system.time(res[15,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 3, method = "GP"))
#'
#'# recursive algorithm, different target dimensions and methods
#'system.time(res[16,] <- ahD_Comb(x, z = z, ind = 1:l, target = 1))
#'system.time(res[17,] <- ahD_Comb(x, z = z, ind = 1:l, target = 2, method = "single"))
#'system.time(res[18,] <- ahD_Comb(x, z = z, ind = 1:l, target = 2, method = "multiple"))
#'system.time(res[19,] <- ahD_Comb(x, z = z, ind = 1:l, target = 3, method = "nGP"))
#'system.time(res[20,] <- ahD_Comb(x, z = z, ind = 1:l, target = 3, method = "GP"))
#'
#'print(paste("Number of different results:  ",sum(apply(res, 2, max) - apply(res, 2, min) > 1e-13)))


ahD_Rec <- function(x, mass = NULL, z = NULL, ind = NULL, target = 2, method = c("standard", "adaptive", "single", "multiple", "nGP", "GP"), nThreads = 0) {
  method <- match.arg(method)
  method <- switch(method,
                   standard = 0,
                   adaptive = 1,
                   single = 2,
                   multiple = 3,
                   nGP = 4,
                   GP = 5)
  if ((target == 2) && !(method %in% c(0,1,2,3))) {
    stop("For target dimension target = 2 method must be one of 'standard', 'adaptive', 'single', 'multiple'!")
  }
  if ((target == 3) && !(method %in% c(0,4,5))) {
    stop("For target dimension target = 3 method must be one of 'standard', 'GP', 'nGP'!")
  }
  if (missing(mass)) {
    .C(
      "aHD_R_int",
      as.double(t(x)),
      as.integer(nrow(x)),
      as.integer(ncol(x)),
      if (!is.null(z)) as.double(t(z)) else NA,
      as.integer(length(z) / as.integer(ncol(x))),
      if (!is.null(ind)) as.integer(ind - 1) else NA,
      as.integer(length(ind)),
      as.integer(target),
      result = integer(length(z) / as.integer(ncol(x)) + length(ind)),
      as.integer(method),
      as.integer(nThreads),
      NAOK = TRUE,
      PACKAGE = "sphericalDepth"
    )$result
  }
  else
    if (is.integer(mass)) {
      .C(
        "aHD_R_int_val",
        as.double(t(x)),
        as.integer(mass),
        as.integer(nrow(x)),
        as.integer(ncol(x)),
        if (!is.null(z)) as.double(t(z)) else NA,
        as.integer(length(z) / as.integer(ncol(x))),
        if (!is.null(ind)) as.integer(ind - 1) else NA,
        as.integer(length(ind)),
        as.integer(target),
        result = integer(length(z) / as.integer(ncol(x)) + length(ind)),
        as.integer(method),
        as.integer(nThreads),
        NAOK = TRUE,
        PACKAGE = "sphericalDepth"
      )$result
    }
  else {
    .C(
      "aHD_R_double_val",
      as.double(t(x)),
      as.double(mass),
      as.integer(nrow(x)),
      as.integer(ncol(x)),
      if (!is.null(z)) as.double(t(z)) else NA,
      as.integer(length(z) / as.integer(ncol(x))),
      if (!is.null(ind)) as.integer(ind - 1) else NA,
      as.integer(length(ind)),
      as.integer(target),
      result = double(length(z) / as.integer(ncol(x)) + length(ind)),
      as.integer(method),
      as.integer(nThreads),
      NAOK = TRUE,
      PACKAGE = "sphericalDepth"
    )$result
  }
}

