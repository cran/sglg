
#set.seed(8142031)
#n1 <- 2
#n2 <- 10
#n3 <- 2
#n <- n1 + n2 + n3

#sqrtA <- matrix(runif(n * n), n, n)
#A <- t(sqrtA) %*% sqrtA
#b <- runif(n)
#true <- solve(A) %*% b
#x0 <- true + runif(n, -1.5, 1.5)
#ps <- c(n1, n2, n3)
#app_sol <- blockgs(A, b, x0, ps)
#(app_sol$x - true)/true

# Algunas simulaciones

#R <- 1000
#results<- rep(0,R)

#for(i in 1:R){
#   x0 <- true + runif(n, -1.5, 1.5)
#   results[i] <- blockgs(A, b, x0, ps)$iter
#}

#results
#summary(results)
#hist(results)
#quantile(results,0.99)
# > 4000

blockgs <- function(A, b, x0, ps, iter, tol) {
    if (missingArg(iter))
        iter <- 5000
    if (missingArg(tol))
        tol <- 0.0001

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
    return(list(x = x, iter = m, cond = round(cond, 5)))
}
