#' Build the basis matrix and the penalty matrix  of cubic B-spline basis.
#'
#'\code{deBoor} builds the basis matrix and penalty matrix to approximate a smooth function using
#' cubic B-spline cubic.
#' @param knots a set of internal knot.
#' @param t a vector of values.
#'
#' @return nknot number of knots.
#' @return knots set of knots.
#' @return N  basis matrix.
#' @return K penalty matrix.
#'
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models.  Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' set.seed(1)
#' t_1 <- runif(120)
#' range(t_1)
#' t_2 <- t_1 + 2  #runif(120,2,3)
#' range(t_2)
#' knot <- 10
#' dB1 <- deBoor2(t_1,knot)
#' dB2 <- deBoor2(t_2,knot)
#' dB1$knots
#' dB2$knots
#' plot(0,0,xlim=c(-0.5,3.5))
#' points(dB1$knots,rep(0,length(dB1$knots)),pch=20)
#' delta <- dB2$knots[1] - dB1$knots[1]
#' points(dB2$knots-delta,rep(0,length(dB2$knots)),pch=2,col= 'blue')
#' dB1$K
#' dB2$K
#' zeros <- vector()
#' plot(t_1,dB1$N[,1],pch=20)
#' for(j in 1:knot){
#' points(t_1,dB1$N[,j],pch=20,col=j)
#' zeros[j] <- sum(dB1$N[,j]==0)
#' }
#' zeros/120
#' cond_tNN <- vector()
#' KnotS <- 3:50
#' for(j in KnotS){
#' dB1 <- deBoor2(t_1,j)
#' print(dB1$knots[2]- dB1$knots[1])
#' min_max <- range(eigen(t(dB1$N)%*%dB1$N)$values)
#' cond_tNN[j-2] <- min_max[1]/min_max[2]
#' }
#' cond_tNN
#' plot(KnotS,cond_tNN,pch=20,ylim=c(0,0.07))
#' @export deBoor2
deBoor2 <- function(t, knots) {

    bspline <- function(t, knots, i, m) {
        if (m == -1) {
            final <- length(knots) - 1
            if (i == final) {
                output <- as.numeric(t <= knots[i + 1] & t >= knots[i])
            } else output <- as.numeric(t < knots[i + 1] & t >= knots[i])
        } else {
            z0 <- (t - knots[i])/(knots[i + m + 1] - knots[i])
            z1 <- (knots[i + m + 2] - t)/(knots[i + m + 2] - knots[i + 1])
            output <- z0 * bspline(t, knots, i, m - 1) + z1 * bspline(t,
                knots, i + 1, m - 1)
        }
        return(output)
    }

    Bspl.N <- function(t, knots) {
        N <- matrix(0, length(t), length(knots) - 4)
        J <- length(knots) - 4
        for (j in 1:J) {
            N[, j] <- bspline(t, knots, j, 2)
        }
        return(N)
    }

    spl.S <- function(knots) {
        output <- diff(diag(knots), difference = 2)
        output <- t(output) %*% output
        return(output)
    }


    vknt <- function(x, knots) {
        int <- range(x)
        grid <- seq(int[1], int[2], length = knots)
        h <- grid[2] - grid[1]
        ts <- seq(int[1] - 2 * h, int[2] + 2 * h, by = h)
        return(ts)
    }
    valueknt <- vknt(as.matrix(t), knots)

    Nbas <- function(x, knots) {
        Nb <- Bspl.N(as.matrix(x), valueknt)
        return(Nb)
    }

    N <- Nbas(t, knots)
    K <- spl.S(knots)
    return(list(knot = knots, knots = valueknt, N = N, K = K))
}
