test_that("reduction works", {

  testLowerDim <- function(n, m1, m2, dim, dimLow) {

    m <- m1 + m2
    x <- matrix(rnorm(n * dimLow), nrow = n)
    x <- x / sqrt(rowSums(x * x))

    z1 <- matrix(rnorm(m1 * dimLow), nrow = m1)
    z1 <- z1 / sqrt(rowSums(z1 * z1))

    Z2 <- matrix(rnorm(m2 * dim), nrow = m2)
    Z2 <- Z2 / sqrt(rowSums(Z2 * Z2))

    X <- cbind(x, matrix(0, ncol = dim - dimLow, nrow = n))
    Z1 <- cbind(z1, matrix(0,ncol = dim - dimLow, nrow = m1))

    v <- rnorm(dim)
    v <- v / sqrt(sum(v*v))

    X <- X - 2 * X %*% v %*% t(v)      # Householdertransformation für X
    Z1 <- Z1 - 2 * Z1 %*% v %*% t(v)   # Householdertransformation für Z1

    Z <- matrix(nrow = m1+m2, ncol = dim)
    r <- sample(1:(m1+m2), m1)
    Z[r,] <- Z1
    Z[-r,] <- Z2

    res1 <- integer(m1+m2+n)
    res1 <- ahD_Comb(x=X,z=Z,ind=1:n)

    res2 <- integer(m1 + m2)
    res2[r] <- ahD_Comb(x=x,z=z1)
    res2[-r] <- ahD_Comb(x=cbind(x,0), z = c(rep(0, dimLow), 1))
    res2 <- c(res2, ahD_Comb(x=x,ind=1:n))

    return(abs(sum(res1 - res2)))
  }

  result <- replicate(1000,
                testLowerDim(n = 200, m1 = 1000, m2 = 100, dim = 4, dimLow = 2))

  expect_equal(sum(result), 0)
})
