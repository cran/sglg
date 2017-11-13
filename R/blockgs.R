
blockgs <- function(A, b, x0, ps, iter) {
    if (missingArg(iter)) 
        iter <- 100
    
    tol <- 1e-05
    k <- length(ps)
    cond <- 1
    l <- 1
    x <- x0
    
    while (l <= iter & cond > tol) {
        aps <- c(0, ps)
        for (j in 1:k) {
            sp <- sum(aps[1:j])
            l <- 1 + sp
            r <- sp + ps[j]
            x[l:r] <- solve(A[l:r, l:r]) %*% (b[l:r] - A[l:r, 
                -(l:r)] %*% x[-(l:r)])
        }
        l <- l + 1
        cond <- sqrt(sum((x - x0)^2))
        x0 <- x
    }
    if (l > iter) {
        stop("Sorry, convergence was not successful.")
    }
    return(list(x = round(x, digits = 5), iter = l, cond = cond))
}
