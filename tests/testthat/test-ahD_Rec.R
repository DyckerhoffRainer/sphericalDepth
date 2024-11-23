test_that("ahD_Rec works", {

  d <- 4
  n <- 200
  m <- 200
  l <- 200
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
  system.time(res[ 1,] <- ahD_Rec(x, z = z, ind = 1:l, target = 1))
  system.time(res[ 2,] <- ahD_Rec(x, z = z, ind = 1:l, target = 2, method = "single"))
  system.time(res[ 3,] <- ahD_Rec(x, z = z, ind = 1:l, target = 2, method = "multiple"))
  system.time(res[ 4,] <- ahD_Rec(x, z = z, ind = 1:l, target = 3, method = "nGP"))
  system.time(res[ 5,] <- ahD_Rec(x, z = z, ind = 1:l, target = 3, method = "GP"))
  # recursive algorithm, pass a vector of probabilitiues to mass, different target dimensions and methods
  system.time(res[ 6,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 1))
  system.time(res[ 7,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 2, method = "single"))
  system.time(res[ 8,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 2, method = "multiple"))
  system.time(res[ 9,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 3, method = "nGP"))
  system.time(res[10,] <- ahD_Rec(x, mass = prob, z = z, ind = 1:l, target = 3, method = "GP"))
  # multiply by n since the depths are since we want to compare with the integer version of the depth
  res[6:10,] <- res[6:10,] * n
  # recursive algorithm, pass a vector of multiplicieties to mass, different target dimensions and methods
  system.time(res[11,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 1))
  system.time(res[12,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 2, method = "single"))
  system.time(res[13,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 2, method = "multiple"))
  system.time(res[14,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 3, method = "nGP"))
  system.time(res[15,] <- ahD_Rec(x, mass = count, z = z, ind = 1:l, target = 3, method = "GP"))
  # combinatorial algorithm, w/o argument 'mass', different target dimensions and methods
  system.time(res[16,] <- ahD_Comb(x, z = z, ind = 1:l, target = 1))
  system.time(res[17,] <- ahD_Comb(x, z = z, ind = 1:l, target = 2, method = "single"))
  system.time(res[18,] <- ahD_Comb(x, z = z, ind = 1:l, target = 2, method = "multiple"))
  system.time(res[19,] <- ahD_Comb(x, z = z, ind = 1:l, target = 3, method = "nGP"))
  system.time(res[20,] <- ahD_Comb(x, z = z, ind = 1:l, target = 3, method = "GP"))
  print(paste("Number of different results: ",sum(apply(res, 2, max) - apply(res, 2, min) > 1e-13)))

  expect_equal(sum(apply(res, 2, max) - apply(res, 2, min) > 1e-13), 0)
})
