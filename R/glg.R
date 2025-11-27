#'Fitting multiple linear Generalized Log-gamma Regression Models
#'
#'\code{glg} is used to fit a multiple linear regression model suitable for analysis of data sets in which the response variable is continuous, strictly positive, and asymmetric.
#'In this setup, the location parameter of the response variable is explicitly modeled by a linear function of the parameters.
#'
#' @param formula a symbolic description of the systematic component of the model to be fitted.
#' @param data a data frame with the variables in the model.
#' @param shape an optional value for the shape parameter of the error distribution of a generalized log-gamma distribution. Default value is 0.2.
#' @param Tolerance an optional positive value, which represents the convergence criterion. Default value is 1e-04.
#' @param Maxiter an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1e03.
#' @param format an optional string value that indicates if you want a simple or a complete report of the estimating process. Default value is 'complete'.
#' @param envelope an optional and internal logical value that indicates if the glg function will be employed for build an envelope plot. Default value is 'FALSE'.
#'
#' @return mu a vector of parameter estimates associated with the location parameter.
#' @return sigma estimate of the scale parameter associated with the model.
#' @return lambda estimate of the shape parameter associated with the model.
#' @return interval estimate of a 95\% confidence interval for each estimate parameters associated with the model.
#' @return Deviance the deviance associated with the model.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @references Cardozo C.A.,  Paula G., and Vanegas L. (2022). Generalized log-gamma additive partial linear models with P-spline smoothing. Statistical Papers.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' set.seed(22)
#' rows <- 200
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(rows, 0, 1)
#' X <- cbind(x1,x2)
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#'
#' ######################
#' #                    #
#' # Extreme value case #
#' #                    #
#' ######################
#'
#' t_lambda <- -1
#' error <- rglg(rows, 0, 1, t_lambda)
#' y1 <- error
#' y1 <- X %*%t_beta + t_sigma*error
#' data.example <- data.frame(y1,X)
#' fit <- glg(y1 ~ x1 + x2 - 1, data=data.example)
#' fit$condition
#' logLik(fit)
#' summary(fit)
#' deviance_residuals(fit)
#' #############################
#' #                           #
#' # Normal case: A limit case #
#' #                           #
#' #############################
#' # When the parameter lambda goes to zero the GLG tends to a standard normal distribution.
#' set.seed(8142031)
#' y1 <- X %*%t_beta + t_sigma * rnorm(rows)
#' data.example <- data.frame(y1, X)
#' fit0 <- glg(y1 ~ x1 + x2 - 1,data=data.example)
#' logLik(fit0)
#' fit0$AIC
#' fit0$mu
#'
#' ############################################
#' #                                          #
#' #  A comparison with a normal linear model #
#' #                                          #
#' ############################################
#'
#' fit2 <- lm(y1 ~ x1 + x2 - 1,data=data.example)
#' logLik(fit2)
#' AIC(fit2)
#' coefficients(fit2)
#' @import methods
#' @importFrom Rcpp sourceCpp
#' @export glg

glg = function(formula, data, shape=0.75, Tolerance=5e-05, Maxiter=1000, format='complete', envelope= FALSE) {
    if (missingArg(formula)) {
        stop("The formula argument is missing.")
    }
    if (missingArg(data)) {
        stop("The data argument is missing.")
    }

    if (!is.data.frame(data))
      stop("The data argument must be a data frame.")

    data <- model.frame(formula, data = data)

    X <- model.matrix(formula, data = data)
    y <- model.response(data)
    p <- ncol(X)
    p1 <- p + 1
    p2 <- p + 2
    n <- nrow(X)
    # Initial values
    fit0 <- .lm.fit(X,y)
    beta0 <- coef(fit0)
    ########################################################
    s_fit0 <-  sort(fit0[[3]])
    breaks <- round(n*c(0.05,0.95))
    res_sam <- s_fit0[ seq(breaks[1],breaks[2],length=10)]
    sigma0 <- sqrt(sum((n*0.1)*(res_sam^2))/(n-p))
    ########################################################
    med <- median(y)
    mean <- mean(y)
    skew_val <- (mean - med)/med
    lambda0 <- shape

    if(skew_val > 0.05){
       lambda0 <- shape
    }else if(skew_val > 0.025){
       lambda0 <- -0.25
    }else if(skew_val > - 0.055){
       lambda0 <- 0.25
    }else{
      lambda0 <- -shape
    }
    ########################################################
    # Some fixed matrices
    One <- matrix(1, n, 1)
    I <- diag(1, n)
    t_X <- t(X)
    t_XX <- t_X %*% X
    t_XOne <- t_X %*% One
    # Fisher Matrix
    I_theta <- function(sigm,lambd) {
      inv_sigm <- 1/sigm
      inv_sigm2 <- inv_sigm^2
      inv_lamb <- 1/lambd
      inv_lamb2 <- inv_lamb^2
      output <- cbind(inv_sigm2*t_XX, inv_sigm2*u_lambda(lambd)*t_XOne, -inv_sigm * inv_lamb*u_lambda(lambd)*t_XOne)
      output <- rbind(output, c(output[1:p, p1], I_22(n,sigm,lambd), I_23(n, sigm, lambd)))
      output <- rbind(output, c(output[1:p, p2], output[p1, p2], n*inv_lamb2*K_1(lambd)))
      return(output)
    }
    ### Defining mu function
    mu <- function(bet) {
        return(X %*% bet)
    }
    ### First step: The residuals
    eps <- function(bet, sigm) {
        return( (y - mu(bet))/sigm)
    }
    ### Second step: The matrix D
    D <- function(epsilon, lambd) {
        return(diag(as.vector(exp(lambd * epsilon))))
    }
    ### Third step: The matrix W
    W <- function(bet, sigm, lambd) {
        return(I - D(eps(bet, sigm), lambd))
    }
    t_One <- t(One)
    U_theta <- function(bet, sigm, lambd) {
      inv_sigm <- 1/sigm
      inv_sigm2 <- inv_sigm^2
      inv_lamb <- 1/lambd
      inv_lamb2 <- inv_lamb^2
      W <- W(bet, sigm, lambd)
      epsilons <- eps(bet, sigm)
      eta_lambd <- inv_lamb * (1 + 2 * inv_lamb2 * (digamma(inv_lamb2) + 2 * log(abs(lambd)) - 1))
      Ds <- D(eps(bet, sigm), lambd)
      return(
             c(-inv_lamb*inv_sigm* t_X %*% (W %*% One),
               -inv_sigm*n - inv_lamb*inv_sigm * t_One %*% (W %*% epsilons),
              n*eta_lambd - inv_lamb2 * t_One %*% ( epsilons  + Ds %*% (epsilons - 2*inv_lamb))))
    }
################################################################################################################################
    ## LOG-LIKELIHOOD

    loglikglg <- function(bet, sigm, lambd) {
        epsilon <- eps(bet, sigm)
        invlamb <- 1/lambd
        return(n * log(c_l(lambd)/sigm) + invlamb * t_One %*% (epsilon - invlamb * (D(epsilon, lambd) %*% One)))
    }
################################################################################################################################
    newpar <- function(bet, sigm, lambd) {
      scores <- U_theta(bet, sigm, lambd)
      I <- I_theta(sigm, lambd)
      ini <- c(bet,sigm,lambd)
      dir <-  solve(I) %*% scores
      llglg_ini <- loglikglg(bet, sigm, lambd)
      condition <- -1
      M <- 0
      while (condition < 0 & M <= 10) {
          new <- ini + (0.8**M)*dir
          llglg_new <- loglikglg(new[1:p], new[p1], new[p2])
          condition <- llglg_new - llglg_ini
          M <- M + 1

      }
      return(new)
    }
    optimum <- function(bet, sigm, lambd){
      fin <- c(bet,sigm,lambd)
      condition <- 1
      l <- 1
      while (condition > Tolerance & l <= Maxiter){
        temp <- newpar(fin[1:p], fin[p1], fin[p2])
        llglg_temp <- loglikglg(temp[1:p], temp[p1], temp[p2])
        llglg_fin <- loglikglg(fin[1:p], fin[p1], fin[p2])
        condition <- llglg_temp - llglg_fin
        if (condition > 0){
          fin <- temp
        }
        l <- l + 1
      }
      if (condition < Tolerance) {
        return(list(est = fin, loglik = llglg_fin, cond = condition, conv = TRUE, iter = l))
      }
      else{
        print("The optimization process was not successful.")
        return(list(conv = FALSE))
        stop("")
      }
    }
    output <- optimum(beta0, sigma0, lambda0)
    if (format == 'simple')
       {
        if(envelope==TRUE){
           mu_est <- X %*% output$est[1:p]
           sgn <- sign(y - mu_est)
           ordresidual <- eps(output$est[1:p], output$est[p1])
           dev <- sgn * 1.4142 * sqrt((1/output$est[p2]^2) * exp(output$est[p2] * ordresidual) - (1/output$est[p2]) * ordresidual - (1/output$est[p2])^2)
           output <- list(mu = output$est[1:p], sigma = output$est[p1], lambda = output$est[p2], Knot = 0, rdev = dev, conv = output$conv, censored = FALSE)
           class(output) = "sglg"
           return(output)}
        else{
          return(list(mu = output$est[1:p],
                  sigma = output$est[p1],
                  lambda = output$est[p2]))}
    }
#---------------------------------------------------------------------------------------------------------------------------------------------------------
        llglg <- output$loglik
        mus <- output$est[1:p]
        sigma <- output$est[p1]
        lambda <- output$est[p2]
        aic <- -2 * llglg + 2 * p2
        bic <- -2 * llglg + log(n) * p2
        covar <- I_theta(sigma, lambda)
        mu_est <- X %*% mus
        ordresidual <- eps(mus, sigma)
        sgn <- sign(y - mu_est)
        ilambda <- 1/lambda
        ilambda2 <- ilambda^2
        dev <- sgn * 1.4142 * sqrt(ilambda2 * exp(lambda * ordresidual) - ilambda * ordresidual - ilambda2)
        devian <- sum(dev^2)
        part2 <- (sigma*ilambda) * (digamma(ilambda2) - log(ilambda2))
        y_est <- mu_est + part2
        scores <- U_theta(mus, sigma, lambda)
        output <- list(formula = formula,
                       size = n,
                       mu = mus,
                       sigma = sigma,
                       lambda = lambda,
                       y = y,
                       p = p,
                       X = X,
                       Knot = 0,
                       llglg = llglg,
                       scores = scores,
                       AIC = aic,
                       BIC = bic,
                       deviance = devian,
                       rdev = dev,
                       Itheta = covar,
                       convergence = output$conv,
                       condition = output$cond,
                       iterations = output$iter,
                       mu_est = mu_est,
                       y_est = y_est,
                       rord = ordresidual,
                       semi = FALSE,
                       censored = FALSE)
                       class(output) = "sglg"
        return(output)

}
