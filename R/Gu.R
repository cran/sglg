#'Tool to build the basis matrix and the penalty matrix of natural cubic splines.
#'
#'\code{Gu} builds the basis matrix and penalty matrix to approximate a smooth function using
#'natural cubic splines based on the Gu basis form.
#'
#' @param t the covariate.
#' @param knot a integer value that represent the number of knots of the natural cubic spline.

#' @return nknot number of knots.
#' @return knots set of knots.
#' @return N  basis matrix.
#' @return K penalty matrix.

#' @references Wood, S. (2006) Generalized additive models: An R introduction. Chapman and Hall.
#' @references Carlos Alberto Cardozo Delgado. Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' t <- runif(1000)
#' knot <- 6
#' N_gu <- Gu(t,knot)
#' @export Gu
Gu <- function(t, knot) {
    r3 <- function(t, z) {
        output <- ((z - 0.5)^2 - 1/12) * ((t - 0.5)^2 - 1/12)/4 - ((abs(t -
            z) - 0.5)^4 - 0.5 * (abs(t - z) - 0.5)^2 + 7/240)/24
        return(output)
    }

    spl.N <- function(t, tk) {
        N <- cbind(t,outer(t, tk, FUN = r3))
        return(N)
    }

    spl.S <- function(tk) {
        S <- matrix(0, q, q)
        S[2:q, 2:q] <- outer(tk, tk, FUN = r3)
        return(S)
    }

    tks <- function(t, knot) {
        ts <- t - min(t)
        ts <- ts/max(ts)
        le <- as.numeric(knot)
        ts <- seq(0, 1, length = le)
        knt <- le - 1
        output <- ts[1:knt]
        return(output)
    }

    knots <- tks(t, knot)
    q <- length(knots) + 1
    N <- spl.N(t, knots)
    K <- spl.S(knots)
    return(list(nknot = knot, knots = knots, N = N, K = K))
}

