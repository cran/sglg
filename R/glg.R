#'Fitting multiple linear Generalized Log-gamma Regression Models
#'
#'\code{glg} is used to fit a multiple linear regression model suitable for analysis of data sets in which the response variable is continuous, strictly positive, and asymmetric.
#'In this setup, the location parameter of the response variable is explicitly modeled by a linear function of the parameters.
#'
#' @param formula a symbolic description of the systematic component of the model to be fitted. See details for further information.
#' @param data an optional data frame, list containing the variables in the model.
#' @param shape an optional value for the shape parameter of the error distribution of a generalized log-gamma distribution. Default value is 1.
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
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' set.seed(21)
#' rows <- 200
#' columns <- 2
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(columns, 0, 1)
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
#' y1 <- X %*%t_beta + t_sigma * error
#' data.example <- data.frame(y1,X)
#' fit <- glg(y1 ~ x1 + x2 - 1,data=data.example)
#' logLik(fit)
#' summary(fit)
#' deviance_residuals(fit)
#' #############################
#' #                           #
#' # Normal case: A limit case #
#' #                           #
#' #############################
#'
#' t_lambda <- 0.001
#' set.seed(8142031)
#' error <- rglg(rows, 0, 1, t_lambda)
#' y1 <- X %*%t_beta + t_sigma * error
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
#' @export glg

glg = function(formula, data, shape, Tolerance, Maxiter, format, envelope= FALSE) {
    if (missingArg(formula)) {
        stop("The formula argument is missing.")
    }
    if (missingArg(data)) {
        stop("The data argument is missing.")
    }

    if (missingArg(Tolerance))
        Tolerance <- 1e-04
    if (missingArg(Maxiter))
        Maxiter <- 1000
    if (missingArg(shape))
        shape <- 1
    if (missingArg(format))
        format <- 'complete'
    if (class(data) == "list")
        data <- as.data.frame(data)

    data <- model.frame(formula, data = data)

    X <- model.matrix(formula, data = data)
    y <- model.response(data)
    p <- ncol(X)
    n <- nrow(X)

    # Initial values

    fit0 <- lm(formula, data = data)
    beta0 <- coefficients(fit0)
    sigma0 <- sum(fit0$residuals^2)/fit0$df.residual
    lambda0 <- shape

    # Some fixed matrizes
    One <- matrix(1, n, 1)
    I <- diag(1, n)

    ## Defining the components of the FSIM

    t_X <- t(X)
    I_11 <- function(sigm) {
        output <- (1/(sigm^2)) * t_X %*% X
        return(output)
    }

    t_XOne <- t_X %*% One
    I_12 <- function(sigm, lambd) {
        output <- (1/(sigm^2)) * u_lambda(lambd) * t_XOne
        return(output)
    }

    I_13 <- function(sigm, lambd) {
        output <- (-1/(sigm * lambd)) * u_lambda(lambd) * t_XOne
        return(output)
    }

    I_33 <- function(sigm, lambd) {
        output <- (n/(lambd^2)) * K_1(lambd)
        return(output)
    }

    I_theta <- function(sigm,lambd) {
        output <- matrix(0, p + 2, p + 2)
        output[1:p, 1:p] = I_11(sigm)
        output[1:p, (p + 1)] = I_12(sigm,lambd)
        output[1:p, (p + 2)] = I_13(sigm,lambd)
        output[(p + 1), 1:p] = t(output[1:p, (p + 1)])
        output[(p + 1), (p + 1)] = I_22(n,sigm,lambd)
        output[(p + 1), (p + 2)] = I_23(n,sigm,lambd)
        output[(p + 2), 1:p] = t(I_13(sigm,lambd))
        output[(p + 2), (p + 1)] = t(I_23(n,sigm,lambd))
        output[(p + 2), (p + 2)] = I_33(sigm,lambd)
        return(output)
    }

    # SCORE FUNCTIONS

    ## Weight matrix for the score functions

    ### Defining mu function

    mu <- function(bet) {
        output <- X %*% bet
        return(output)
    }

    ### First step: The residuals

    eps <- function(bet, sigm) {
        epsilon <- (y - mu(bet))/sigm
        return(epsilon)
    }

    ### Second step: The matrix D

    D <- function(epsilon, lambd) {
        D_eps <- diag(as.vector(exp(lambd * epsilon)))
        return(D_eps)
    }

    ### Third step: The matrix W

    W <- function(bet, sigm, lambd) {
        output <- I - D(eps(bet, sigm), lambd)
        return(output)
    }

    ## Final step: Score functions

    U_beta <- function(bet, sigm, lambd) {
        output <- (-1/(lambd * sigm)) * t_X %*% W(bet, sigm, lambd) %*%One
        return(output)
    }

    t_One <- t(One)
    U_sigma <- function(bet, sigm, lambd) {
        output <- -(1/sigm) * n - (1/(lambd * sigm)) * t_One %*% W(bet,sigm, lambd) %*% eps(bet, sigm)
        return(output)
    }

    U_lambda <- function(bet, sigm, lambd) {
        invlamb <- 1/lambd
        eta_lambd <- (invlamb) * (1 + 2 * (invlamb^2) * (digamma(invlamb^2) + 2 * log(abs(lambd)) - 1))
        Ds <- D(eps(bet, sigm), lambd)
        epsilons <- eps(bet, sigm)
        t_OneDs <- t_One %*% Ds
        output <- n * eta_lambd - (invlamb^2) * t_One %*% epsilons + (2 * invlamb^3) * t_OneDs %*% One - (invlamb^2) * t_OneDs %*% epsilons
        return(output)
    }

    U_theta <- function(bet, sigm, lambd) {
        output <- matrix(1, p + 2, 1)
        output[1:p] <- U_beta(bet, sigm, lambd)
        output[p + 1] <- U_sigma(bet, sigm, lambd)
        output[p + 2] <- U_lambda(bet, sigm, lambd)
        return(output)
    }

    ## LOG-LIKELIHOOD

    loglikglg <- function(bet, sigm, lambd) {
        epsilon <- eps(bet, sigm)
        output <- n * log(c_l(lambd)/sigm) + (1/lambd) * t_One %*% epsilon - (1/lambd^2) * t_One %*% D(epsilon, lambd) %*% One
        return(output)
    }

    ## THE ESTIMATES

    newpar <- function(bet, sigm, lambd) {
        output <- matrix(0, p + 2, 2)
        output[1:p, 1] <- bet
        output[p + 1, 1] <- sigm
        output[p + 2, 1] <- lambd
        new <- output[, 1]
        output[, 2] <- output[, 1] + solve(I_theta(sigm, lambd)) %*% U_theta(bet,
            sigm, lambd)
        M = 2
        while (output[p + 1, 2] < 0 & M < Maxiter) {
            output[, 2] = 0.99 * output[, 2] + 0.01 * output[, 1]
            M = M + 1
        }
        llglg = loglikglg(output[1:p, 2], output[p + 1, 2], output[p + 2,
            2])
        condition = llglg - loglikglg(output[1:p, 1], output[p + 1, 1], output[p +
            2, 1])
        if (condition > 0) {
            new = output[, 2]
        }
        M = 2
        while (condition < 0 & M < Maxiter) {
            output[, 2] = 0.99 * output[, 2] + 0.01 * output[, 1]
            condition = loglikglg(output[1:p, 2], output[p + 1, 2], output[p +
                2, 2]) - loglikglg(output[1:p, 1], output[p + 1, 1], output[p +
                2, 1])
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
        new = matrix(0, p + 2, Maxiter)
        new[, 1] = newpar(bet, sigm, lambd)
        ouput = new[, 1]
        new[, 2] = newpar(new[1:p, 1], new[(p + 1), 1], new[(p + 2), 1])
        l = 2
        llglg = loglikglg(new[1:p, l], new[(p + 1), l], new[(p + 2), l])
        condition = llglg - loglikglg(new[1:p, 1], new[p + 1, 1], new[p +
            2, 1])
        if (condition > 0) {
            output = new[, l]
        }
        while (condition > Tolerance & l < Maxiter) {
            l = l + 1
            new[, l] = newpar(new[1:p, (l - 1)], new[(p + 1), (l - 1)], new[(p +
                2), (l - 1)])
            llglg = loglikglg(new[1:p, l], new[(p + 1), l], new[(p + 2),
                l])
            condition = llglg - loglikglg(output[1:p], output[(p + 1)], output[(p +
                2)])

            if (condition > 0) {
                output = new[, l]
            }
        }
        if (condition < Tolerance & l < Maxiter) {
            return(list(est = output, cond = condition, conv = TRUE, iter = l))
        }
        if (l >= Maxiter | scores > 0.05) {
            conv = FALSE
            return(list(conv = conv))
            stop("")
        }
    }

    output <- optimum(beta0, sigma0, lambda0)

    if (abs(output$est[p + 2]) < 0.09) {
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
           ordresidual <- eps(output$est[1:p], output$est[p + 1])
           dev <- sgn * sqrt(2) * sqrt((1/output$est[p + 2]^2) * exp(output$est[p + 2] * ordresidual) - (1/output$est[p + 2]) * ordresidual - (1/output$est[p + 2])^2)
           output <- list(mu = output$est[1:p], sigma = output$est[p +1], lambda = output$est[p + 2], Knot = 0, rdev = dev, conv = output$conv, censored = FALSE)
           class(output) = "sglg"
           return(output)}
        else{
            mu_est <- X %*% output$est[1:p]
            output <- list(mu = output$est[1:p], sigma = output$est[p +1], lambda = output$est[p + 2], conv = output$conv)
            class(output) = "sglg"
            return(output)}
    }

    if (output$conv == TRUE) {
        conv <- output$conv
        iter <- output$iter
        condition <- output$cond
        output <- output$est
        llglg <- loglikglg(output[1:p], output[(p + 1)], output[(p + 2)])
        aic <- -2 * llglg + 2 * (p + 2)
        bic <- -2 * llglg + log(n) * (p + 2)
        aic2 <- aic + 2 * sum(y)
        covar <- matrix(0, p + 2, p + 2)
        covar <- I_theta(output[p + 1], output[p + 2])
        ste <- sqrt(diag(solve(covar)))
        zs <- output/ste
        pval <- 1 - (pnorm(abs(zs)) - pnorm(-abs(zs)))
        mu_est <- X %*% output[1:p]
        ordresidual <- eps(output[1:p], output[p + 1])
        sgn <- sign(y - mu_est)
        dev <- sgn * sqrt(2) * sqrt((1/output[p + 2]^2) * exp(output[p + 2] * ordresidual) - (1/output[p + 2]) * ordresidual - (1/output[p + 2])^2)
        devian <- sum(dev^2)
        part2 <- ((output[p + 1])/(output[p + 2])) * (digamma((1/output[p +
            2])^2) - log((1/output[p + 2])^2))
        y_est <- mu_est + part2
        scores <- U_theta(output[1:p], output[p + 1], output[p + 2])
        inter <- matrix(0, p + 2, 2)
        inter[, 1] <- as.matrix(output - 1.96 * ste)
        inter[, 2] <- as.matrix(output + 1.96 * ste)
        good_fit <- gfit(ordresidual, output[p + 2])
        output <- list(formula = formula, size = n, mu = output[1:p], sigma = output[p +
            1], lambda = output[p + 2], y = y, p = p, X = X, Knot = 0, llglg = llglg,
            scores = scores, AIC = aic, BIC = bic, AIC2 = aic2, deviance = devian,
            rdev = dev, Itheta = covar, st_error = ste, z_values = zs, p.values = pval,
            interval = inter, goodnessoffit = good_fit, convergence = conv,
            condition = condition, iterations = iter, mu_est = mu_est, y_est = y_est, rord = ordresidual,
            semi = FALSE, censored = FALSE)
        class(output) = "sglg"
        return(output)
    }

}
