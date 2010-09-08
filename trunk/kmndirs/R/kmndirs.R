kmndirs <-
function(x, k, nruns = 100L, maxiter = 10L)
{
    ## Currently only dense matrix support.
    x <- as.matrix(x)
    ## Normalize.
    x <- x / sqrt(rowSums(x ^ 2))

    n <- nrow(x)

    .C(R_kmndirs, x, n, ncol(x), as.integer(k), as.integer(nruns),
       as.integer(maxiter), ids = integer(n))$ids
}
