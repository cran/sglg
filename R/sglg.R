#'Fitting semi-parametric generalized log-gamma regression models
#'
#'\code{sglg} is used to fit a semi-parametric regression model suitable for analysis of data sets in which the response variable is continuous, strictly positive, and asymmetric.
#'In this setup, the location parameter of the response variable is explicitly modeled by semi-parametric functions, whose nonparametric components may be approximated by
#'natural cubic splines or cubic P-splines.
#'
#' @param formula a symbolic description of the systematic component of the model to be fitted. See details for further information.
#' @param npc a matrix with the nonparametric variables of the systematic part of the model to be fitted. Must be included the names of each variables.
#' @param basis a name of the cubic spline basis to be used in the model. Supported basis include deBoor and Gu basis
#'  which are a B-spline basis and a natural cubic spline basis, respectively.
#' @param data an optional data frame, list containing the variables in the model.
#' @param shape an optional value for the shape parameter of the error distribution of a generalized log-gamma distribution. Default value is 0.2.
#' @param method There are two possibles algorithms to estimate the parameters. The default algorithm is 'FS' Fisher-Scoring,
#' the other option is 'GSFS' an adequate combination between the block matrix version of non-linear Gauss-Seidel algorithm and Fisher-Scoring algorithm.
#' @param alpha0 is a vector of positive values for the smoothing parameters alpha. Default vector with 1 in each entry.
#' @param nknts is a vector of the number of knots in each non-linear component of the model.
#' @param Tolerance an optional positive value, which represents the convergence criterion. Default value is 5e-05.
#' @param Maxiter an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1e03.
#' @param format an optional string value that indicates if you want a simple or a complete report of the estimating process. Default value is 'complete'.
#'
#' @return mu a vector of parameter estimates associated with the location parameter.
#' @return sigma estimate of the scale parameter associated with the model.
#' @return lambda estimate of the shape parameter associated with the model.
#' @return interval estimate of a 95\% confidence interval for each estimate parameters associated with the model.
#' @return Deviance the deviance associated with the model.

#' @references Carlos A. Cardozo, G. Paula and L. Vanegas. Semi-parametric generalized log-gamma regression models. In preparation.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' set.seed(1)
#' rows<- 200
#' t_beta <- c(0.5,2)
#' t_sigma <- 0.5
#' t_lambda <- 1
#' x1 <- runif(rows,-3,3)
#' x2 <- rbinom(rows,1,0.5)
#' X <- cbind(x1,x2)
#' t <- as.matrix((2*1:rows - 1)/(2*rows))
#' colnames(t) <- "t"
#' f_t <- cos(4*pi*t)
#' error <- rglg(rows,0,1,t_lambda)
#' y <- X %*%t_beta + f_t + t_sigma*error
#' colnames(y) <- "y"
#' data <- data.frame(y,X,t)
#' fit1 <- sglg(y ~ x1 + x2 - 1,npc=t,data=data,basis = "deBoor",alpha0 = 0.1)
#' logLik(fit1) # -195.4538
#' quantile_residuals(fit1)
#' fit2 <- sglg(y ~ x1 + x2 - 1,npc=t,data=data,basis = "Gu",alpha0=0.5)
#' logLik(fit2)
#' #################################################
#' # An example with two non-parametric components #
#' #################################################
#' set.seed(2)
#' t_2 <- as.matrix(rnorm(rows,sd=0.5))
#' colnames(t_2) <- 't_2'
#' f_t_2 <- exp(t_2)
#' error <- rglg(rows,0,1,t_lambda)
#' y_2 <- X %*%t_beta + f_t + f_t_2 + t_sigma*error
#' colnames(y_2) <- 'y_2'
#' data2 <- data.frame(y_2,X,t,t_2)
#' npcs <- cbind(t,t_2)
#' fit3 <- sglg(y_2 ~ x1 + x2 - 1, npc=npcs, data=data2, alpha0 = c(0.45,0.65))
#' logLik(fit3)
#' #############################################################################
#' @import methods
#' @export sglg
#'
sglg = function(formula, npc, basis, data, shape=0.2, method, alpha0, nknts, Tolerance=5e-05, Maxiter=1000,format='complete') {
  if (missingArg(formula)) {
    stop("The formula argument is missing.")
  }
  if (missingArg(npc)) {
    stop("This kind of model need at least one non-parametric component.")
  }
  if (missingArg(data)) {
    stop("The data argument is missing.")
  }

  if (missingArg(method))
    method <- "FS"

  k <- dim(npc)[2]
  if (missingArg(basis))
    basis <- rep("deBoor", k)

  if (class(data) == "list")
    data <- as.data.frame(data)

  ######################################################################################################################################
  ######################################################################################################################################

  data1 <- model.frame(formula, data = data)
  X <- model.matrix(formula, data = data1)
  y <- model.response(data1)

  p <- ncol(X)
  n <- nrow(X)

  Knot <- vector()
  XX <- cbind(X, npc)

  if(missingArg(nknts)){
    op1 <- floor(n^(1/3))
    intknt <- function(x) {
      op2 <- length(as.numeric(levels(factor(x))))
      knt <- min(op1, op2)
      if (knt < 3) {
        stop("This covariate has not at least three different values.")
      }
      return(knt)
    }

    for (i in 1:k) {
      Knot <- append(Knot,intknt(XX[, (p + i)]))
    }
  }
  else{
    if(min(nknts)>2)
      Knot <- nknts
    else
      stop("Each covariate must have at least three knots.")
  }

  Tknot <- sum(Knot)

  ############################################################################################################
  formula2 <- formula
  for (j in 1:k) {
    formul <- paste(".~. + ", colnames(npc)[j])
    formul <- as.formula(formul)
    formula2 <- update(formula2, formul)
  }
  ############################################################################################################################################################

  # Initial values

  fit0    <- glg(formula2, shape = shape, data = data, format='simple')
  beta0   <- fit0$mu[1:p]
  g0s     <- rep(0,Tknot)
  sigma0  <- fit0$sigma
  lambda0 <- fit0$lambda

  if(missingArg(alpha0)) {
    alpha0 <- rep(1,k)
  }

  # Some fixed matrizes
  One <- matrix(1, n, 1)
  Ident <- diag(1, n)

  ## THE FISHER INFORMATION MATRIX

  N <- X
  K <- matrix(0, p + Tknot, p + Tknot)

  for (j in 1:k) {
    if(basis[j]=="deBoor")
      output <- deBoor2(npc[, j], Knot[j])
    else
      output <- Gu(npc[, j],Knot[j])

    N <- cbind(N, output$N)
    Knot1 <- c(0, Knot)
    l <- p + 1 + sum(Knot1[1:j])
    r <- p + sum(Knot[1:j])
    K[l:r, l:r] <- alpha0[j]*output$K
  }

  t_N <- t(N)

  I_gammas <- function(sigm) {
    return((1/(sigm^2)) * t_N %*% N + K)
  }

  I_gammassigma <- function(sigm, lambd) {
    return((1/sigm^2) * u_lambda(lambd) * t_N %*% One)
  }

  I_gammaslambda <- function(sigm, lambd) {
    return(-(sigm/lambd) * I_gammassigma(sigm, lambd))
  }

  ### Defining mu function

  mu <- function(bet, g) {
    return(N %*% c(bet, g))
  }

  ### First step: The residuals

  eps <- function(bet, g, sigm) {
    return((y - mu(bet, g))/sigm)
  }

  ### Second step: The matrix D

  D <- function(epsil, lambd) {
    return(diag(as.vector(exp(lambd * epsil)), n, n))
  }

  ### Third step: The matrix W

  W <- function(bet, g, sigm, lambd) {
    return(Ident - D(eps(bet, g, sigm), lambd))
  }

  ## The score functions

  t_One <- t(One)
  U_sigma <- function(bet, g, sigm, lambd) {
    return(-(1/sigm) * n - (1/(lambd * sigm)) * t_One %*% W(bet,g, sigm, lambd) %*% eps(bet, g, sigm))
  }

  U_lambda <- function(bet, g, sigm, lambd) {
    invlamb <- 1/lambd
    invlambtwo <- invlamb^2
    epsilons <- eps(bet, g, sigm)
    eta_lambd <- (invlamb) * (1 + 2 * (invlambtwo) * (digamma(invlambtwo) + 2 * log(abs(lambd)) - 1))
    Ds <- D(eps(bet, g, sigm), lambd)
    epsilons <- eps(bet, g, sigm)
    output <- n * eta_lambd - (invlambtwo) * t_One %*% epsilons + (2 * invlamb^3) * t_One %*% Ds %*% One - (invlambtwo) * t_One %*% Ds %*% epsilons
    return(output)
  }

  U_gammas <- function(bet, g, sigm, lambd){
    return((-1/(sigm * lambd)) * t_N %*% W(bet, g, sigm, lambd) %*% One - K %*% c(bet, g))
  }

  U_theta <- function(bet, g, sigm, lambd) {
    return(c(U_gammas(bet, g, sigm, lambd),U_sigma(bet, g, sigm, lambd),U_lambda(bet, g, sigm,lambd)))
  }

  ## Defining the components of the FIM

  I_33 <- function(sigm, lambd) {
    return((n/(lambd^2)) * K_1(lambd))
  }

  I_sl <- function(sigm, lambd) {
    return(matrix(c(I_22(n,sigm, lambd), I_23(n,sigm, lambd), I_23(n,sigm,lambd), I_33(sigm, lambd)), 2, 2))
  }

  I_theta <- function(sigm, lambd) {
    I_gs <- I_gammassigma(sigm, lambd)
    I_gl <- I_gammaslambda(sigm, lambd)
    output <- cbind(I_gammas(sigm), I_gs, I_gl)
    output1 <- cbind(I_gs, I_gl)
    output1 <- cbind(t(output1), I_sl(sigm, lambd))
    output <- rbind(output, output1)
    return(output)
  }

  ## LOG-LIKELIHOOD

  loglikglg <- function(bet, g, sigm, lambd) {
    epsilon <- eps(bet, g, sigm)
    invlambd <- 1/lambd
    output <- n * log(c_l(lambd)/sigm) + invlambd * t_One %*% epsilon - invlambd^2 * t_One %*% D(epsilon, lambd) %*% One -0.5 * t(c(bet, g)) %*% K %*% c(bet, g)
    return(output)
  }

  newpar <- function(bet, g, sigm, lambd) {
    output <- matrix(0, p + Tknot + 2, Maxiter)
    output[, 1] <- c(bet, g, sigm, lambd)
    new <- output[, 1]
    llglg <- loglikglg(new[1:p], new[(p + 1):(p + Tknot)], new[p + Tknot + 1], new[p + Tknot + 2])

    if (method == "FS") {
      output[, 2] <- output[, 1] + solve(I_theta(sigm, lambd)) %*% U_theta(bet, g, sigm, lambd)
      condition <- loglikglg(output[1:p, 2], output[(p + 1):(p + Tknot), 2], output[p + Tknot + 1, 2], output[p + Tknot + 2, 2]) - llglg
      if (condition > 0) {
        new <- output[, 2]
      }
      condition <- 1
      l <- 2
      while (abs(condition) > Tolerance & l < Maxiter) {
        l <- l + 1
        output[, l] <- output[, (l - 1)] + solve(I_theta(output[(p + Tknot + 1), (l - 1)], output[(p + Tknot + 2), (l - 1)])) %*% U_theta(output[1:p, (l - 1)], output[(p + 1):(p + Tknot), (l - 1)], output[(p + Tknot + 1), (l - 1)], output[(p + Tknot + 2), (l - 1)])
        llglg <- loglikglg(new[1:p], new[(p + 1):(p + Tknot)], new[p + Tknot + 1], new[p + Tknot + 2])
        condition <- loglikglg(output[1:p, l], output[(p + 1):(p + Tknot), l], output[p + Tknot + 1, l], output[p + Tknot + 2, l]) - llglg
        if (abs(condition) > 0) {
          new <- output[, l]
          llglg <- loglikglg(new[1:p], new[(p + 1):(p + Tknot)], new[p + Tknot + 1], new[p + Tknot + 2])
        }
      }
    }

    if (method == "GSFS") {
      ps <- c(p, Knot, 2)
      b <- I_theta(sigm, lambd) %*% output[, 1] + U_theta(bet,
                                                          g, sigm, lambd)
      output[, 2] <- blockgs(I_theta(sigm, lambd), b, output[,
                                                             1], ps)$x
      condition <- loglikglg(output[1:p, 2], output[(p + 1):(p + Tknot),2], output[p + Tknot + 1, 2], output[p + Tknot + 2, 2]) - llglg
      if (condition > 0) {
        new <- output[, 2]
      }
      condition <- 1
      l <- 2
      while (condition > Tolerance & l < Maxiter) {
        l <- l + 1
        b <- I_theta(output[(p + Tknot + 1), (l - 1)], output[(p +
                                                                 Tknot + 2), (l - 1)]) %*% output[, (l - 1)] + U_theta(output[1:p,
                                                                                                                              (l - 1)], output[(p + 1):(p + Tknot), (l - 1)], output[(p +
                                                                                                                                                                                        Tknot + 1), (l - 1)], output[(p + Tknot + 2), (l - 1)])
        output[, l] <- blockgs(I_theta(output[(p + Tknot + 1), (l -
                                                                  1)], output[(p + Tknot + 2), (l - 1)]), b, output[,
                                                                                                                    (l - 1)], ps)$x
        llglg <- loglikglg(new[1:p], new[(p + 1):(p + Tknot)], new[p +
                                                                     Tknot + 1], new[p + Tknot + 2])
        condition <- loglikglg(output[1:p, l], output[(p + 1):(p +
                                                                 Tknot), l], output[p + Tknot + 1, l], output[p + Tknot +
                                                                                                                2, l]) - llglg
        if (condition > 0) {
          new <- output[, l]
          llglg <- loglikglg(new[1:p], new[(p + 1):(p + Tknot)],
                             new[p + Tknot + 1], new[p + Tknot + 2])
        }
      }
    }
    return(list(est = new, ll = llglg, cond = condition, iter = l))
  }

  ## Effective degree freedom - EDF

  inv_root_A <- function(A) {
    e <- eigen(A)
    V <- e$vectors
    output <- solve(V %*% sqrt(diag(e$values)) %*% t(V))
    return(output)
  }

  edf <- function(sigm) {
    inv_root_tNN <- inv_root_A(t_N %*% N)
    output <- diag(1, p + sum(Knot)) + (sigm^2) * inv_root_tNN %*% K %*% inv_root_tNN
    output <- sum(diag(solve(output)))
    return(output)
  }

  Conv <- FALSE
  num.iter <- 1
  masterf <- function() {
    news <- newpar(beta0, g0s, sigma0, lambda0)
    if (news$iter >= Maxiter)
        stop('The algorithm did not converge!')
    else {
      Conv <- TRUE
      num.iter <- news$iter
      cond <- news$cond
      llglg <- news$ll
      df <- edf(news$est[(p + Tknot + 1)])
      aic <- -2 *(llglg - df - 2)
      bic <- -2 *llglg - log(n) * (df + 2)
      return(list(est = news$est, df = df, llglg = llglg, AIC = aic, BIC = bic,
                  Conv = Conv, iter = num.iter, cond = cond))
    }
  }

  if (format == 'simple')
  {
    output <- masterf()
    output_est <- output$est
    y_est <- N %*% output_est[1:(p + Tknot)]
    part2 <- ((output_est[p + Tknot + 1])/(output_est[p + Tknot + 2])) * (digamma((1/output_est[p + Tknot + 2])^2) - log((1/output_est[p + Tknot + 2])^2))
    y_est <- y_est + part2
    output <- list(alpha=alpha0,AIC=output$AIC,BIC=output$BIC,y=y,df=output$df,y_est=y_est)
    class(output) <- "sglg"
    return(output)
  }

  total_optimum <- function(start) {
    output0 <- c(start,masterf()$AIC)
    output1 <- output0[1:k]
    output3 <- masterf()
    df <- output3$df
    dfnpc <- df - p
    output <- output3$est
    llglg <- output3$llglg
    scores <- U_theta(output[1:p], output[(p + 1):(p + Tknot)], output[p + sum(Knot) + 1], output[p + Tknot + 2])
    covar <- I_theta(output[p + Tknot + 1], output[p + Tknot + 2])
    inter <- matrix(0, p + Tknot + 2, 2)
    scovar <- solve(covar)
    val <- diag(scovar)
    if (min(val) > 0) {
      ste <- sqrt(val)
      inter[, 1] <- as.matrix(output - 1.96 * ste)
      inter[, 2] <- as.matrix(output + 1.96 * ste)
      zs <- abs(output/ste)
      pval <- 1 - (pnorm(zs) - pnorm(-zs))
      pval2 <- pval[-((p + 1):(p + Tknot))]
      as <- output[(p + 1):(p + Tknot)]
    }

    y_est <- N %*% output[1:(p + Tknot)]
    ordresidual <- eps(output[1:p], output[(p + 1):(p + Tknot)], output[p +
                                                                          Tknot + 1])
    sgn <- sign(y - y_est)
    dev <- sgn * sqrt(2) * sqrt((1/output[p + Tknot + 2]^2) * exp(output[p +
                                                                           Tknot + 2] * ordresidual) - (1/output[p + Tknot + 2]) * ordresidual -
                                  (1/output[p + Tknot + 2])^2)
    Devian <- sum(dev^2)

    #good_fit <- gfit(ordresidual, output[p + Tknot + 2])
    #goodnessoffit = good_fit

    part2 <- ((output[p + Tknot + 1])/(output[p + Tknot + 2])) * (digamma((1/output[p +
                                                                                      Tknot + 2])^2) - log((1/output[p + Tknot + 2])^2))
    y_est2 <- y_est + part2

    return(list(formula = formula, npc = npc, basis =basis, size = n, mu = output[1:(p +
                                                                                       Tknot)], sigma = output[p + Tknot + 1], lambda = output[p + Tknot +
                                                                                                                                                 2], y = y, X = X, p = p, N = N, Knot = Knot, rord = ordresidual,
                rdev = dev, interval = inter, llglg = llglg, AIC = output3$AIC, BIC = output3$BIC,
                scores = scores, Itheta = covar, scovar = scovar, st_error = ste,
                Z_values = zs, p.values = pval, alpha = output1, d.f.model = df,
                d.f.npc = dfnpc, deviance = Devian,
                convergence = output3$Conv, condition = output3$cond, iterations = output3$iter,
                basis = basis, semi = TRUE, censored = FALSE, mu_est = y_est, y_est = y_est2))
  }
  output <- total_optimum(alpha0)
  class(output) <- "sglg"
  return(output)
}
