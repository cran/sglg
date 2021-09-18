
# @importFrom Rcpp sourceCpp
# sourceCpp('src/cpp_code.cpp')

blockgs <- function(A, b, x0, ps, iter, tol) {
    if (missingArg(iter))
        iter <- 5000
    if (missingArg(tol))
        tol <- 0.01

    k <- length(ps)
    cond <- 1
    m <- 1
    x <- x0

    mchol <- function(A, b) {
        R <- chol(A)
        x <- backsolve(R, b, transpose = TRUE)
        x <- backsolve(R, x)
        return(x)
    }

    while (m <= iter & cond > tol) {
        aps <- c(0, ps)
        for (j in 1:k) {
            sp <- sum(aps[1:j])
            l <- 1 + sp
            r <- sp + ps[j]
            bb <- (b[l:r] - A[l:r, -(l:r)] %*% x[-(l:r)])
            x[l:r] <- mchol(A[l:r, l:r], bb)
        }
        m <- m + 1
        cond <- sqrt(sum((x - x0)^2))
        x0 <- x
    }
    if (m > iter) {
        stop("Sorry, convergence was not successful.")
    }
    return(list(x = x, iter = m, cond = cond))
}
