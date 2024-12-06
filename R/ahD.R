#' Angular halfspace depth
#'
#' This function computes the angular halfspace depth (it is not assumed that the points are in general position).
#'
#'@param x a matrix of format \code{n}-times-\code{d} that contains the data points with respect to which the depth has to be computed
#'@param mass a vector of length \code{n} that contains either the probabilities or the (integer) multiplicities of the \code{n} data points
#'@param z a matrix of format \code{m}-times-\code{d} containing the points whose depth should be computed
#'@param ind a vector of length \code{l} containing the indices of the data points \code{x} whose depth should be computed
#'@param alg a string, possible values are \code{"comb"} and \code{"rec"}.
#'Denotes the algorithm that is used to compute the signed halfspace depth in the target dimension.
#'
#'If \code{alg="comb"}, the combinatorial algorithm is used to compute the angular halfspace depth. If \code{alg="rec"} the recursive algorithm is used. For details, see Dyckerhoff and Nagy (2024).
#'@param target the dimension to which the data is projected, possible values are \code{1}, \code{2}, or \code{3}
#'@param method a string, possible values are \code{"standard"}, \code{"adaptive"}, \code{"single"}, \code{"multiple"}, \code{"nGP"}, \code{"GP"}.
#'Denotes the variant of the algorithm that is used to compute the signed halfspace depth in the target dimension.
#'
#'For \code{target=2}, only \code{"standard"}, \code{"adaptive"}, \code{"single"}, \code{"multiple"} are valid.
#'For \code{target=2}, \code{method="standard"} is the same as \code{method="adaptive"}, meaning that depending on the number of points for which the depth has to be computed, either \code{"single"} or \code{"multiple"} is selected.
#'
#'For \code{target=3}, only \code{"standard"}, \code{"nGP"}, \code{"GP"} are valid.
#'For \code{target=3}, \code{method=} \code{"standard"} is the same as \code{method="nGP"}, meaning that the points do not have to be in general position. If \code{method="GP"}, then it is assumed that all points are in general position, so that a faster algorithm can be used.
#'
#'@param nThreads the number of threads to use in a multiprocessor environment. For the default \code{0}, the number of threads is automatically chosen depending on the used hardware.
#'@return a vector of length \code{m+l} containing the angular halfspace depths of the points \code{z} and the data points \code{x} whose indices are given in \code{ind}.
#'If the argument \code{mass} is missing or if an integer vector is supplied for \code{mass}, then the depths are given as integer values, i.e., multiplied by \code{n}.
#'If a double vector is supplied for the argument \code{mass}, then the depths are given as doubles.
#'
#'@details
#'Unless \code{target=3} and \code{method="GP"}, the algorithm does not assume that the points are in general position. The points in \code{z} may coincide with points in \code{x}. \code{z} or \code{ind} may be missing. If both \code{z} and \code{ind} are missing, thenthe depths of all points in \code{x} are calculated.
#'
#'Regarding the choice of parameters \code{alg}, \code{target}, and \code{method}, the default values usually will give the best results.
#'In most cases the combinatorial algorithm (\code{alg="comb"}) called with the default values for \code{target} and \code{method} will give the best results.
#'
#'
#'
#'@references Dyckerhoff, R., and Nagy, S. (2024). Exact computation of angular halfspace depth.
#'
#'@author Rainer Dyckerhoff
#'
#'@export
#'@useDynLib sphericalDepth, .registration=TRUE
#'
#'@examples
#'d <- 4
#'n <- 50
#'m <- 50
#'l <- 50
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
#'res <- matrix(nrow=30, ncol=m+l)
#'
#'# combinatorial algorithm, w/o argument 'mass', different target dimensions and methods
#'res[ 1,] <- ahD(x, z = z, ind = 1:l, alg = "comb", target = 1)
#'res[ 2,] <- ahD(x, z = z, ind = 1:l, alg = "comb", target = 2, method = "single")
#'res[ 3,] <- ahD(x, z = z, ind = 1:l, alg = "comb", target = 2, method = "multiple")
#'res[ 4,] <- ahD(x, z = z, ind = 1:l, alg = "comb", target = 3, method = "nGP")
#'res[ 5,] <- ahD(x, z = z, ind = 1:l, alg = "comb", target = 3, method = "GP")
#'
#'# combinatorial algorithm, pass a vector of probabilities to mass, different target dimensions
#'# and methods
#'res[ 6,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "comb", target = 1)
#'res[ 7,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "comb", target = 2, method = "single")
#'res[ 8,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "comb", target = 2, method = "multiple")
#'res[ 9,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "comb", target = 3, method = "nGP")
#'res[10,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "comb", target = 3, method = "GP")
#'# multiply by n since the depths are since we want to compare with the integer version of the depth
#'res[6:10,] <- res[6:10,] * n
#'
#'# combinatorial algorithm, pass a vector of multiplicities to mass, different target dimensions
#'# and methods
#'res[11,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "comb", target = 1)
#'res[12,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "comb", target = 2, method = "single")
#'res[13,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "comb", target = 2, method = "multiple")
#'res[14,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "comb", target = 3, method = "nGP")
#'res[15,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "comb", target = 3, method = "GP")
#'
#'# recursive algorithm, w/o argument 'mass', different target dimensions and methods
#'res[16,] <- ahD(x, z = z, ind = 1:l, alg = "rec", target = 1)
#'res[17,] <- ahD(x, z = z, ind = 1:l, alg = "rec", target = 2, method = "single")
#'res[18,] <- ahD(x, z = z, ind = 1:l, alg = "rec", target = 2, method = "multiple")
#'res[19,] <- ahD(x, z = z, ind = 1:l, alg = "rec", target = 3, method = "nGP")
#'res[20,] <- ahD(x, z = z, ind = 1:l, alg = "rec", target = 3, method = "GP")
#'
#'# recursive algorithm, pass a vector of probabilities to mass, different target dimensions
#'# and methods
#'res[21,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "rec", target = 1)
#'res[22,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "rec", target = 2, method = "single")
#'res[23,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "rec", target = 2, method = "multiple")
#'res[24,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "rec", target = 3, method = "nGP")
#'res[25,] <- ahD(x, mass = prob, z = z, ind = 1:l, alg = "rec", target = 3, method = "GP")
#'# multiply by n since the depths are since we want to compare with the integer version of the depth
#'res[21:25,] <- res[21:25,] * n
#'
#'# recursive algorithm, pass a vector of multiplicities to mass, different target dimensions
#'# and methods
#'res[26,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "rec", target = 1)
#'res[27,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "rec", target = 2, method = "single")
#'res[28,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "rec", target = 2, method = "multiple")
#'res[29,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "rec", target = 3, method = "nGP")
#'res[30,] <- ahD(x, mass = count, z = z, ind = 1:l, alg = "rec", target = 3, method = "GP")
#'
#'print(paste("Number of different results:  ",sum(apply(res, 2, max) - apply(res, 2, min) > 1e-13)))

ahD <- function(x, mass = NULL, z = NULL, ind = if (missing(z)) 1:nrow(as.matrix(x)) else NULL, alg = c("comb", "rec"), target = 2, method = c("standard", "adaptive", "single", "multiple", "nGP", "GP"), nThreads = 0) {
  alg <- match.arg(alg)
  alg <- switch(alg,
                comb = 0,
                rec = 1)
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
      "aHD_int",
      as.double(t(x)),
      as.integer(nrow(x)),
      as.integer(ncol(x)),
      if (!is.null(z)) as.double(t(z)) else NA,
      as.integer(length(z) / as.integer(ncol(x))),
      if (!is.null(ind)) as.integer(ind - 1) else NA,
      as.integer(length(ind)),
      result = integer(length(z) / as.integer(ncol(x)) + length(ind)),
      as.integer(alg),
      as.integer(target),
      as.integer(method),
      as.integer(nThreads),
      NAOK = TRUE,
      PACKAGE = "sphericalDepth"
    )$result
  }
  else
    if (is.integer(mass)) {
      .C(
        "aHD_int_val",
        as.double(t(x)),
        as.integer(mass),
        as.integer(nrow(x)),
        as.integer(ncol(x)),
        if (!is.null(z)) as.double(t(z)) else NA,
        as.integer(length(z) / as.integer(ncol(x))),
        if (!is.null(ind)) as.integer(ind - 1) else NA,
        as.integer(length(ind)),
        result = integer(length(z) / as.integer(ncol(x)) + length(ind)),
        as.integer(alg),
        as.integer(target),
        as.integer(method),
        as.integer(nThreads),
        NAOK = TRUE,
        PACKAGE = "sphericalDepth"
      )$result
    }
  else {
    .C(
      "aHD_double_val",
      as.double(t(x)),
      as.double(mass),
      as.integer(nrow(x)),
      as.integer(ncol(x)),
      if (!is.null(z)) as.double(t(z)) else NA,
      as.integer(length(z) / as.integer(ncol(x))),
      if (!is.null(ind)) as.integer(ind - 1) else NA,
      as.integer(length(ind)),
      result = double(length(z) / as.integer(ncol(x)) + length(ind)),
      as.integer(alg),
      as.integer(target),
      as.integer(method),
      as.integer(nThreads),
      NAOK = TRUE,
      PACKAGE = "sphericalDepth"
    )$result
  }
}


