#'Plotting a natural cubic splines or P-splines.
#'
#'\code{plotnpc} displays a graph of a fitted nonparametric effect, either natural cubic spline or P-spline, from an object of class sglg.
#'
#' @param fit an object of the class sglg. This object is returned from the call to sglg(), ssurvglg().

#' @references Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.
#' @references Wood, S. (2006). Additive generalized models: An R introduction. Chapman and Hall.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @import graphics
#' @examples
#' rows <- 175 # Number of observations
#' columns <- 2 # Number of parametric components
#' library(ssym)
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#' t_lambda <- 1
#' t_knot1 <- 7
#' ts1 <- seq(0, 1, length = t_knot1)
#' t_g1 <- 0.4 * sin(pi * ts1)
#'
#' BasisN <- function(n, knot) {
#'           N <- matrix(0, n, knot)
#'           m <- n/knot
#'           block <- matrix(1, m, 1)
#'           for (i in 1:knot) {
#'           l <- (i - 1) * m + 1
#'           r <- i * m
#'           N[l:r, i] <- block }
#'           return(N)
#'           }
#' s_N1 <- BasisN(rows, length(ts1))
#' x3 <- s_N1 %*% ts1
#' colnames(x3) <- 'x3'
#' set.seed(8142031)
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(rows, 0, 1)
#' X <- cbind(x1,x2)
#' error <- rglg(rows, 0, 1, t_lambda)
#' y1 <- X %*%t_beta + + s_N1 %*% t_g1 + t_sigma * error
#' data.example <- data.frame(y1,X,x3)
#' fit1 <- sglg(y1 ~ x1 + x2 - 1, npc=x3, method='FS',data=data.example)
#' plotnpc(fit1)
#' @export plotnpc

plotnpc <- function(fit) {
  if (fit$semi == FALSE) {
    stop("Sorry, for this kind of model it is not available this option.")
  }
  npc <- fit$npc
  y <- fit$y
  X <- fit$X
  p <- fit$p
  Knot <- fit$Knot
  betas <- fit$mu[1:p]
  as <- fit$mu[(p + 1):(p + Knot)]
  lambda <- fit$lambda
  rord <- fit$rord
  rdev <- fit$rdev
  y_est <- fit$y_est2
  scovar <- fit$scovar
  t_npc <- as.numeric(levels(factor(as.matrix(npc))))
  N_t <- deBoor2(t_npc, Knot)$N
  g_t <- N_t %*% as
  scovarred <- scovar[(p + 1):(p + Knot), (p + 1):(p + Knot)]
  Var <- N_t %*% (scovarred %*% t(N_t))
  nval <- diag(Var)
  nste <- sqrt(nval)
  n1 <- length(nste)
  percentil <- 0.025/n1
  quantil <- abs(qnorm(percentil))
  upper_g_t <- g_t + quantil * nste
  lower_g_t <- g_t - quantil * nste

  cloud <- y - X %*% betas
  xrange <- range(t_npc)
  yrange <- range(c(lower_g_t, upper_g_t))
  plot(as.matrix(npc), cloud, ylim = yrange, col = 3, pch = 20, xlab = colnames(npc),
       ylab = "g(x)", main = "Simultaneous 95% confidence intervals")

  f <- splinefun(t_npc, g_t, method = "natural")
  ls(envir = environment(f))
  splinecoef <- get("z", envir = environment(f))
  values <- curve(f, min(t_npc), max(t_npc), col = "black", add = TRUE)

  uf <- splinefun(t_npc, upper_g_t, method = "natural")
  ls(envir <- environment(uf))
  splinecoef <- get("z", envir = environment(uf))
  values <- curve(uf, min(t_npc), max(t_npc), col = "red", add = TRUE)

  lf <- splinefun(t_npc, lower_g_t, method = "natural")
  ls(envir <- environment(lf))
  splinecoef <- get("z", envir = environment(lf))
  values <- curve(lf, min(t_npc), max(t_npc), col = "red", add = TRUE)
}
