% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ahD_Rec.R
\name{ahD_Rec}
\alias{ahD_Rec}
\title{Angular halfspace depth}
\usage{
ahD_Rec(
  x,
  mass = NULL,
  z = NULL,
  ind = NULL,
  target = 2,
  method = c("standard", "adaptive", "single", "multiple", "nGP", "GP"),
  nThreads = 0
)
}
\arguments{
\item{x}{a matrix of format \code{n}-times-\code{d} that contains the data points with respect to which the depth has to be computed}

\item{mass}{a vector of length \code{n} that contains either the probabilities or the (integer) multiplicities of the \code{n} data points}

\item{z}{a matrix of format \code{m}-times-\code{d} containing the points whose depth should be computed}

\item{ind}{a vector of length \code{l} containing the indices of the data points \code{x} whose depth should be computed}

\item{target}{the dimension to which the data is projected, possible values are \code{1}, \code{2}, or \code{3}}

\item{method}{a string, possible values are \code{"standard"}, \code{"adaptive"}, \code{"single"}, \code{"multiple"}, \code{"nGP"}, \code{"GP"}.
Denotes the variant of the algorithm that is used to compute the signed halfspace depth in the target dimension.

For \code{target=2}, only \code{"standard"}, \code{"adaptive"}, \code{"single"}, \code{"multiple"} are valid.
For \code{target=2}, \code{method="standard"} is the same as \code{method="adaptive"}, meaning that depending on the number of points for which the depth has to be computed, either \code{"single"} or \code{"multiple"} is selected.

For \code{target=3}, only \code{"standard"}, \code{"nGP"}, \code{"GP"} are valid.
For \code{target=3}, \code{method=} \code{"standard"} is the same as \code{method="nGP"}, meaning that the points do not have to be in general position. If \code{method="GP"}, then it is assumed that all points are in general position, so that a faster algorithm can be used.}

\item{nThreads}{the number of threads to use in a multiprocessor environment. For the default \code{0}, the number of threads is automatically chosen depending on the used hardware.}
}
\value{
a vector of length \code{m+l} containing the angular halfspace depths of the points \code{z} and the data points \code{x} whose indices are given in \code{ind}.
If the argument \code{mass} is missing or if an integer vector is supplied for \code{mass}, then the depths are given as integer values, i.e., multiplied by \code{n}.
If a double vector is supplied for the argument \code{mass}, then the depths are given as doubles.
}
\description{
This function computes the angular halfspace depth using the recursive algorithm (it is not assumed that the points are in general position). \code{ahD_Rec} is deprecated and should no longer be used. Applications should use \code{ahD} instead.
}
\details{
The routine uses the recursive algorithm of Dyckerhoff and Nagy (2024). Unless \code{target=3} and \code{method="GP"}, it does not assume that the points are in general position. The points in z may coincide with points in x. z or ind may be missing (but not both).

Regarding the choice of parameters \code{target}, and \code{method}, the default values usually will give the best results.
In most cases this routine will be slower than 'ahD_Comb'.
}
\examples{
d <- 4
n <- 50
m <- 50
l <- 50
# simulate data uniformly distributed on the sphere
x <- matrix(rnorm(n*d),nrow=n)
x <- x / sqrt(rowSums(x*x))
z <- matrix(rnorm(m*d),nrow=m)
z <- z / sqrt(rowSums(z*z))
# vector of probabilities
prob <- rep(1/n, n)
# vector of point multiplicities
count <- as.integer(rep(1,n))

res <- matrix(nrow=20, ncol=m+l)

# recursive algorithm, w/o argument 'mass', different target dimensions and methods
res[ 1,] <- ahD_Rec(x, z = z, ind = 1:l, target = 1)
res[ 2,] <- ahD_Rec(x, z = z, ind = 1:l, target = 2, method = "single")
res[ 3,] <- ahD_Rec(x, z = z, ind = 1:l, target = 2, method = "multiple")
res[ 4,] <- ahD_Rec(x, z = z, ind = 1:l, target = 3, method = "nGP")
res[ 5,] <- ahD_Rec(x, z = z, ind = 1:l, target = 3, method = "GP")

# recursive algorithm, pass a vector of probabilitiues to mass, different target dimensions
#and methods
res[ 6,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 1)
res[ 7,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 2, method = "single")
res[ 8,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 2, method = "multiple")
res[ 9,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 3, method = "nGP")
res[10,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 3, method = "GP")
# multiply by n since the depths are since we want to compare with the integer version of the depth
res[6:10,] <- res[6:10,] * n

# recursive algorithm, pass a vector of multiplicieties to mass, different target dimensions
# and methods
res[11,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 1)
res[12,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 2, method = "single")
res[13,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 2, method = "multiple")
res[14,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 3, method = "nGP")
res[15,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 3, method = "GP")

# combinatorial algorithm, different target dimensions and methods
res[16,] <- ahD_Comb(x, z = z, ind = 1:l, target = 1)
res[17,] <- ahD_Comb(x, z = z, ind = 1:l, target = 2, method = "single")
res[18,] <- ahD_Comb(x, z = z, ind = 1:l, target = 2, method = "multiple")
res[19,] <- ahD_Comb(x, z = z, ind = 1:l, target = 3, method = "nGP")
res[20,] <- ahD_Comb(x, z = z, ind = 1:l, target = 3, method = "GP")

print(paste("Number of different results:  ",sum(apply(res, 2, max) - apply(res, 2, min) > 1e-13)))
}
\references{
Dyckerhoff, R., and Nagy, S. (2024). Exact computation of angular halfspace depth.
}
\author{
Rainer Dyckerhoff
}
