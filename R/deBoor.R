#' Build the values of cubic B-spline basis function.

#' @param t a vector of values.
#' @param knots a set of knots distributed along the range of t.
#' @param i i-th basis element in which the t values are going to be evaluated.
#' @param m order of the basis elements. In the case of cubic B-spline basis, m=2.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models.  Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.

bspline = function(t, knots, i, m) {
    if (m == -1) {
        final <- length(knots) - 1
        if (i == final) {
            output = as.numeric(t <= knots[i + 1] & t >= 
                knots[i])
        } else output = as.numeric(t < knots[i + 1] & t >= 
            knots[i])
    } else {
        z0 = (t - knots[i])/(knots[i + m + 1] - knots[i])
        z1 = (knots[i + m + 2] - t)/(knots[i + m + 2] - knots[i + 
            1])
        output = z0 * bspline(t, knots, i, m - 1) + z1 * 
            bspline(t, knots, i + 1, m - 1)
    }
    return(output)
}

#' A model matrix for the B-spline

#' @param t a vector of values.
#' @param knots a set of knots distributed along the range of t.

Bspl.N = function(t, knots) {
    N = matrix(0, length(t), length(knots) - 4)
    J <- length(knots) - 4
    for (j in 1:J) {
        N[, j] = bspline(t, knots, j, 2)
    }
    return(N)
}

#' Penalized B-spline matrix, given the number of internal knots.
#' @param knots a set of internal knot.

spl.S = function(knots) {
    output <- diff(diag(knots), difference = 2)
    output <- t(output) %*% output
    return(output)
}

#' Build a set of knots for a given number of internal knots.
#' @param knots a set of internal knot.
#' @param x a vector of values.

vknt = function(x, knots) {
    int <- range(x)
    grid <- seq(int[1], int[2], length = knots)
    h <- grid[2] - grid[1]
    ts <- seq(int[1] - 2 * h, int[2] + 2 * h, by = h)
    return(ts)
}

#' Build the basis matrix of cubic B-spline basis.
#' @param knots a set of internal knot.
#' @param x a vector of values.

Nbas = function(x, knots) {
    valueknt <- vknt(as.matrix(x), knots)
    Nb <- Bspl.N(as.matrix(x), valueknt)
    return(Nb)
}

#' Build the basis matrix and the penalty matrix  of cubic B-spline basis.
#' @param knots a set of internal knot.
#' @param x a vector of values.

deBoor <- function(x, knots) {
    N <- Nbas(x, knots)
    K <- spl.S(knots)
    return(list(N = N, K = K))
}
