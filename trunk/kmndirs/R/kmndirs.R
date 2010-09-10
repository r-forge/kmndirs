kmndirs <-
function(x, k, nrandom = 100L, maxiter = 10L)
{
    ## Currently only dense matrix support.
    x <- as.matrix(x)
    ## Normalize.
    x <- x / sqrt(rowSums(x ^ 2))

    n <- nrow(x)

    .C(R_kmndirs, x, n, ncol(x), as.integer(k), as.integer(nrandom),
       as.integer(maxiter), ids = integer(n))$ids
}
