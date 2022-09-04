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
#' rows <- 300
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
#' t_lambda <- 1
#' error <- rglg(rows, 0, 1, t_lambda)
#' y1 <- error
#' y1 <- X %*%t_beta + t_sigma*error
#' data.example <- data.frame(y1,X)
#' data.example <- data.frame(y1)
#' fit <- glg(y1 ~ x1 + x2 - 1, data=data.example)
#' logLik(fit) # -449.47  # Time: 14 milliseconds
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

glg = function(formula, data, shape=-0.75, Tolerance=5e-05, Maxiter=1000, format='complete', envelope= FALSE) {
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
    frac <- 0.1
    s_fit0 <-  sort(fit0[[3]])
    breaks <- round(n*c(0.05,0.95))
    res_sam <- s_fit0[seq(breaks[1],breaks[2],length=10)]
    sigma0 <- sqrt(sum((n*frac)*(res_sam^2))/(n-p))

    ########################################################
    med <- median(y)
    mean <- mean(y)
    skew_val <- (mean - med)/med
    if(skew_val > 0.05) lambda0 <- shape
    else if(skew_val > 0.025) lambda0 <- -0.25
    else if(skew_val > - 0.055) lambda0 <- 0.25
    else lambda0 <- -shape

    ########################################################
    # Some fixed matrices
    One <- matrix(1, n, 1)
    I <- diag(1, n)

    ## Defining the components of the FSIM

    t_X <- t(X)
    t_XX <- t_X %*% X
    I_11 <- function(sigm) {
      return ((1/(sigm^2)) * t_XX)
    }

    t_XOne <- t_X %*% One
    I_12 <- function(sigm, lambd) {
        return((1/(sigm^2)) * u_lambda(lambd) * t_XOne)
    }

    I_13 <- function(sigm, lambd) {
        return((-1/(sigm * lambd)) * u_lambda(lambd) * t_XOne)
    }

    I_33 <- function(sigm, lambd) {
        return((n/(lambd^2)) * K_1(lambd))
    }

    I_theta <- function(sigm,lambd) {
        output <- cbind(I_11(sigm),I_12(sigm,lambd),I_13(sigm,lambd))
        output <- rbind(output,c(output[1:p,p1],I_22(n,sigm,lambd),I_23(n,sigm,lambd)))
        output <- rbind(output,c(output[1:p,p2],output[p1,p2],I_33(sigm,lambd)))
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

    ## Final step: Score functions

    U_beta <- function(bet, sigm, lambd) {
        return((-1/(lambd * sigm)) * t_X %*% (W(bet, sigm, lambd) %*%One))
    }

    t_One <- t(One)
    U_sigma <- function(bet, sigm, lambd) {
        return(-(1/sigm) * n - (1/(lambd * sigm)) * t_One %*% (W(bet,sigm, lambd) %*% eps(bet, sigm)))
    }

    U_lambda <- function(bet, sigm, lambd) {
        invlamb <- 1/lambd
        invlamb2 <- invlamb^2
        eta_lambd <- (invlamb) * (1 + 2 * (invlamb2) * (digamma(invlamb2) + 2 * log(abs(lambd)) - 1))
        Ds <- D(eps(bet, sigm), lambd)
        epsilons <- eps(bet, sigm)
        t_OneDs <- t_One %*% Ds
        return(n * eta_lambd - (invlamb2) * t_One %*% epsilons + (2 * invlamb*invlamb2) * t_OneDs %*% One - (invlamb2) * t_OneDs %*% epsilons)
    }

    U_theta <- function(bet, sigm, lambd) {
        return(c(U_beta(bet, sigm, lambd),U_sigma(bet, sigm, lambd),U_lambda(bet, sigm, lambd)))
    }

    ## LOG-LIKELIHOOD

    loglikglg <- function(bet, sigm, lambd) {
        epsilon <- eps(bet, sigm)
        invlamb <- 1/lambd
        return(n * log(c_l(lambd)/sigm) + (invlamb) * t_One %*% epsilon - (invlamb^2) * t_One %*% (D(epsilon, lambd) %*% One))
        #return(n * (lc_l(lambd) - log(sigm)) + (invlamb) * t_One %*% epsilon - (invlamb^2) * t_One %*% (D(epsilon, lambd) %*% One))
    }

    ## THE ESTIMATES

    newpar <- function(bet, sigm, lambd) {
        output <- matrix(0, p2, 2)
        output[, 1] <- c(bet,sigm,lambd)
        new <- output[, 1]
        output[, 2] <- output[, 1] + solve(I_theta(sigm, lambd)) %*% U_theta(bet, sigm, lambd)
        M = 2
        while (output[p1, 2] < 0 & M < Maxiter) {
            output[, 2] = 0.8 * output[, 2] + 0.2 * output[, 1]
            M = M + 1
        }
        llglg = loglikglg(output[1:p, 2], output[p1, 2], output[p2,2])
        condition = llglg - loglikglg(output[1:p, 1], output[p1, 1], output[p2, 1])
        if (condition > 0) {
            new = output[, 2]
        }
        M = 2
        while (condition < 0 & M < Maxiter) {
            output[, 2] = 0.8 * output[, 2] + 0.2 * output[, 1]
            condition = loglikglg(output[1:p, 2], output[p1, 2], output[p2, 2]) - loglikglg(output[1:p, 1], output[p1, 1], output[p2, 1])
            if (condition > 0) {
                new = output[, 2]
            }
            M = M + 1
        }
        return(new)
    }

    ## THE MAIN FUNCTION

    conv <- FALSE
    condition <- 1
    l <- 1

    optimum <- function(bet, sigm, lambd) {
        new = matrix(0, p2, Maxiter)
        new[, 1] = newpar(bet, sigm, lambd)
        ouput = new[, 1]
        new[, 2] = newpar(new[1:p, 1], new[p1, 1], new[p2, 1])
        l = 2
        llglg = loglikglg(new[1:p, l], new[p1, l], new[p2, l])
        condition = llglg - loglikglg(new[1:p, 1], new[p1, 1], new[p2, 1])
        if (condition > 0) {
            output = new[, l]
        }
        while (condition > Tolerance & l < Maxiter) {
            l = l + 1
            new[, l] = newpar(new[1:p, (l - 1)], new[p1, (l - 1)], new[p2, (l - 1)])
            llglg = loglikglg(new[1:p, l], new[p1, l], new[p2,l])
            condition = llglg - loglikglg(output[1:p], output[p1], output[p2])

            if (condition > 0) {
                output = new[, l]
            }
        }
        if (condition < Tolerance & l < Maxiter) {
            return(list(est = output, loglik = llglg, cond = condition, conv = TRUE, iter = l))
        }
        if (l >= Maxiter | scores > 0.05) {
            conv = FALSE
            return(list(conv = conv))
            stop("")
        }
    }

    output <- optimum(beta0, sigma0, lambda0)

    if (abs(output$est[p2]) < 0.2) {
        output <- optimum(beta0, sigma0, -lambda0)
    }

    if (output$conv == FALSE) {
        print("The optimization was not successful.")
        return(0)
    }

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
        ste <- sqrt(diag(solve(covar)))
        zs <- output$est/ste
        pval <- 1 - (pnorm(abs(zs)) - pnorm(-abs(zs)))
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
        delta <- 1.96 * ste
        inter <- cbind(output$est - delta, output$est + delta)
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
                       st_error = ste,
                       z_values = zs,
                       p.values = pval,
                       interval = inter,
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
