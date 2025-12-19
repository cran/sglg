#'Fitting linear generalized log-gamma regression models under the presence of right censored data.
#'
#'\code{survglg} is used to fit a multiple linear regression model in which the response variable is continuous, strictly positive, asymmetric and there are right censored observations.
#'In this setup, the location parameter of the logarithm of the response variable is modeled by a linear model of the parameters.
#'
#' @param formula a symbolic description of the systematic component of the model to be fitted. See details for further information.
#' @param data an optional data frame, list containing the variables in the model.
#' @param shape an optional value for the shape parameter of the model.
#' @param Tolerance an optional positive value, which represents the convergence criterion. Default value is 1e-05.
#' @param Maxiter an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1000.
#' @return mu a vector of parameter estimates asociated with the location parameter.
#' @return sigma estimate of the scale parameter associated with the model.
#' @return lambda estimate of the shape parameter associated with the model.
#' @return interval estimate of a 95\% confidence interval for each estimate parameters associated with the model.
#' @return Deviance the deviance associated with the model.
#' @references Carlos A. Cardozo, G. Paula and L. Vanegas. Semi-parametric accelerated failure time models with generalized log-gamma errors. In preparation.
#' @references Cardozo C.A.,  Paula G., and Vanegas L. (2022). Generalized log-gamma additive partial linear models with P-spline smoothing. Statistical Papers.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' require(survival)
#' rows  <- 200
#' columns <- 2
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#' set.seed(8142031)
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(rows, 0, 1)
#' X <- cbind(x1,x2)
#' s         <- t_sigma^2
#' a         <- 1/s
#' t_ini1    <- exp(X %*% t_beta) * rweibull(rows, scale = s, shape = a)
#' cens.time <- rweibull(rows, 0.75, 20)
#' delta1     <- ifelse(t_ini1 > cens.time, 1, 0)
#' obst1 <- t_ini1
#' obst1[delta1==1] <- cens.time[delta1==1]
#' data.example <- data.frame(obst1,delta1,X)
#' fit3 <- survglg(Surv(log(obst1),delta1) ~ x1 + x2 - 1, data=data.example, shape = 1)
#' fit3$condition
#' fit3$scores
#' summary(fit3)
#' # We can obtain the logLik (in the exp scale) from
#' fit3$llgg
#' delta2 = 1 - delta1
#' fit4 <- survreg(Surv(obst1,delta2) ~ x1 + x2 - 1, data=data.example, dist = 'weibull')
#' summary(fit4)
#' @import Formula
#' @import survival
#' @import methods
#' @export survglg
survglg = function(formula, data, shape, Maxiter=1000, Tolerance=1e-06) {
    if (missingArg(formula))
        stop("The formula argument is missing.")
    if (missingArg(data))
        stop("The data argument is missing.")
    if (!is.data.frame(data))
        stop("The data argument must be a data frame.")
    if (missingArg(shape))
        shape <- 0.5

    data <- model.frame(formula, data = data)
    X <- model.matrix(formula, data = data)
    p <- ncol(X)
    Y <- cbind(data[, 1][, 1], data[, 1][, 2])
    Delta <- as.factor(Y[, 2])
    y <- Y[, 1][Delta == 0]
    XX <- X[Delta == 0, ]
    Delta <- as.numeric(as.vector(Delta))
    per.censo <- 100 * mean(Delta)
    datus <- data.frame(y, XX)
    formula2 <- formula
    formula2 <- Formula(formula2)
    formula2 <- formula(formula2, lhs = 0)
    formula2 <- update(formula2, y ~ .)
 ############################################################################################################################################################

    fit0 <- glg(formula2, data = datus, Tolerance=1e-03, format='simple')
    beta0 <- fit0$mu
    sigma0 <- fit0$sigma
    lambda0 <- shape

    n <- nrow(X)

    # Some fixed matrizes
    p1 = p + 1
    I_n <- diag(1, n)
    One_Delta <- 1 - Delta
    # Defining mu function

    mu <- function(bet) {
        output <- X %*% bet
        return(output)
    }

    ### First step: The residuals

    eps <- function(bet, sigm) {
        epsilon <- (Y[, 1] - mu(bet))/sigm
        return(epsilon)
    }

    ### Second step: The matrix D

    S <- function(y, lambd) {
        ilambd2 <- 1/lambd^2
        s <- pgamma((ilambd2) * exp(lambd * y), ilambd2, lower.tail = FALSE)
        return(s)
    }

    ### Third step: The matrix Wd

    aa <- function(bet, sigm, lambd) {
        epsilon <- eps(bet, sigm)
        output <- (1/lambd^2) * exp(lambd * epsilon)
        return(output)
    }

    ## Final step: Score functions

    U_beta <- function(bet, sigm, lambd) {
        epsil <- eps(bet, sigm)
        aaa <- aa(bet, sigm, lambd)
        ilambd2 <- 1/lambd^2
        bb <- (aaa^(ilambd2) * exp(-aaa))/(gamma(ilambd2) * S(epsil,lambd))
        output <- t(One_Delta * (1/(lambd * sigm))*(exp(lambd * epsil) - 1) + Delta * (lambd/sigm) * bb)%*%X
        return(output)
    }

    U_sigma <- function(bet, sigm, lambd) {
        epsil <- eps(bet, sigm)
        aaa <- aa(bet, sigm, lambd)
        ilambd2 <- 1/lambd^2
        bb <- (aaa^(ilambd2) * exp(-aaa))/(gamma(ilambd2) * S(epsil,lambd))
        part1 <- (1/sigm) * (-1 - (1/lambd) * (epsil - epsil * exp(lambd * epsil)))
        part2 <- (lambd/sigm) * epsil * bb
        output <- sum(One_Delta * part1 + Delta * part2)
        return(output)
    }

    U_theta <- function(bet, sigm, lambd) {
      output <- c(U_beta(bet, sigm, lambd),U_sigma(bet, sigm, lambd))
      return(output)
    }

    # Observational Fisher Matrix

    I_beta <- function(bet, sigm, lambd) {
        epsil <- eps(bet, sigm)
        aaa <- aa(bet, sigm, lambd)
        ilambd2 <- 1/lambd^2
        bb <- (aaa^(ilambd2) * exp(-aaa))/(gamma(ilambd2) * S(epsil,lambd))
        output <- matrix(0, p, p)
        isigm2 <- 1/sigm^2
        for (l in 1:p) {
            for (j in 1:p) {
                output[l, j] <- sum(X[, l] * X[, j] * ((One_Delta) * (-isigm2) *
                  exp(lambd * epsil) + Delta * (isigm2*lambd^2) * bb * (aaa - ilambd2 - bb)))
            }
        }
        return(output)
    }

    I_sigma <- function(bet, sigm, lambd) {
        epsil <- eps(bet, sigm)
        isigm2 <- 1/sigm^2
        part1 <- isigm2 * (1 + (2/lambd) * epsil * (1 - exp(lambd * epsil)) - (epsil^2) * exp(lambd * epsil))
        aaa <- aa(bet, sigm, lambd)
        ilambd2 <- 1/lambd^2
        bb <- (aaa^(ilambd2) * exp(-aaa))/S(epsil, lambd)
        b <- aaa - (ilambd2) - bb/gamma(ilambd2)
        part2 <- (isigm2*(lambd * epsil * bb)) * (epsil * b - 2/lambd)
        output <- sum(One_Delta * part1 + Delta * part2)
        return(output)
    }

    I_sigbet <- function(bet, sigm, lambd) {
        epsil <- eps(bet, sigm)
        expep <- exp(lambd * epsil)
        H <- matrix(0, p, 1)
        for (k in 1:p) {
            H[k] <- sum(One_Delta * X[, k] * (1 - expep - lambd * epsil *expep))
        }
        isigm2 <- 1/(sigm^2)
        part1 <- (1/lambd) * isigm2 * H
        part2 <- 0 * H
        aaa <- aa(bet, sigm, lambd)
        ilambd2 <- 1/lambd^2
        bb <- (aaa^(ilambd2)) * exp(-aaa)
        deno <- gamma(ilambd2) * S(epsil, lambd)
        part21 <- -1/lambd + epsil * (aaa - (ilambd2) - bb/deno)
        for (k in 1:p) {
            H[k] <- sum(Delta * X[, k] * (((lambd*isigm2) * bb/deno) * part21))
        }
        part2 <- H
        output <- part1 + part2
        return(output)
    }

    I_tetha <- function(bet, sigm, lambd) {
        output <- matrix(0, p1, p1)
        output[1:p, 1:p] <- I_beta(bet, sigm, lambd)
        output[p1, p1] <- I_sigma(bet, sigm, lambd)
        output[1:p, p1] <- I_sigbet(bet, sigm, lambd)
        output[p1, 1:p] <- t(output[1:p, p1])
        return(output)
    }

    ## LOG-LIKELIHOOD

    loglikglg <- function(bet, sigm) {
        epsil <- eps(bet, sigm)
        output <- sum(Delta * log(S(epsil, lambda0)) + One_Delta * (log(c_l(lambda0)/sigm) + (1/lambda0) * epsil - (1/lambda0^2) * exp(lambda0 * epsil)))
        return(output)
    }

    newpar <- function(bet, sigm){
        scores <- U_theta(bet, sigm, lambda0)
        I <- I_tetha(bet, sigm, lambda0)
        dir <-  -solve(I) %*% scores
        ini <- c(bet,sigm)
        llglg_ini <- loglikglg(bet, sigm)
        new <- ini + dir
        llglg_new <- loglikglg(new[1:p], new[p1])
        condition <- llglg_new - llglg_ini
        M <- 1
        while (condition <= 0 & M <= 100) {
            new <- ini + (0.5**M)*dir
            llglg_new <- loglikglg(new[1:p], new[p1])
            condition <- llglg_new - llglg_ini
            M <- M + 1
        }
        return(new=new)
    }
    ## THE MAIN FUNCTION
    conv <- FALSE
    l <- 1
    optimum <- function(bet, sigm) {
        #########################################################################
        condition  <- 1
        l <- 1
        output <- c()
        ini <- c(bet, sigm)
        ini_llglg <- loglikglg(bet, sigm)
        while (condition > Tolerance & l <= Maxiter) {
          new <- newpar(ini[1:p], ini[p1])
          new_llglg <- loglikglg(new[1:p], new[p1])
          condition <- new_llglg - ini_llglg
          if (condition > 0) {
            ini <- new
            ini_llglg <- new_llglg
            l <- l + 1
          }
        }
        #########################################################################
        if (l <= Maxiter) {
            scores <- U_theta(ini[1:p], ini[p1],lambda0)
            return(list(est = ini, llglg = ini_llglg, cond = condition, scores = scores, conv = TRUE, iter = l))
        }
        if (l > Maxiter) {
            stop("The convergence was not successful.")
        }
    }

    gfit <- function(resid, lambd) {
      ekm <- survival::survfit(Surv(exp(resid), One_Delta) ~ 1)
      surv <- as.numeric(unlist(as.vector(summary(ekm)[6])))
      Fkm <- 1 - surv
      res <- sort((resid * One_Delta)[Delta == 0])
      Fs <- pglg(res, shape = lambd)
      r_q <- qnorm(Fs)
      diff <- abs(r_q - qnorm(Fkm))
      output <- mean(diff[-length(diff)])
      msurv <- 1 - Fs
      return(list(stat = output, msurv = msurv))
    }

    output <- optimum(beta0, sigma0)

    if (output$conv == TRUE) {
        conv <- output$conv
        iter <- output$iter
        condition <- output$cond
        scores <- output$scores
        llglg <- output$llglg
        llgg <- llglg - sum(y)
        aic <- -2 * llglg + 2 * p1
        bic <- -2 * llglg + log(n) * p1
        aic2 <- aic + 2 * sum(y)
        output <- output$est
        covar <- I_tetha(output[1:p], output[p1], lambda0)
        inter <- matrix(0, p1, 2)
        pval <- c()
        ste <-  c()
        zs <-   c()
        scovar <- solve(-covar)
        val <- diag(scovar)
        if (min(val) > 0) {
            ste <- sqrt(val)
            inter[, 1] <- as.matrix(output - 1.96 * ste)
            inter[, 2] <- as.matrix(output + 1.96 * ste)
            zs <- abs(output/ste)
            pval <- 1 - (pnorm(zs) - pnorm(-zs))
        }

        y_est <- X %*% output[1:p]
        ordresidual <- eps(output[1:p], output[p1])
        sgn <- sign(Y[, 1] - y_est)
        outputp <- lambda0
        ioutputp <- 1/outputp
        ioutputp2 <- ioutputp^2

        #######################################################################################

        dev <- sgn * sqrt(2)*(One_Delta*sqrt(ioutputp2) * exp(outputp*ordresidual) - ioutputp * ordresidual - (ioutputp2) + Delta*(-log(S(ordresidual, outputp))))

        ######################################################################################
        devian <- sum(dev^2)
        part2 <- (output[p1]*ioutputp) * (digamma(ioutputp2) - log(ioutputp2))
        y_est <- y_est + part2
        good_fit <- gfit(ordresidual, outputp)
        msurv <- good_fit$msurv
        good_fit <- good_fit$stat
        output <- list(formula = formula, size = n, per.censo = per.censo,
            p = p, mu = output[1:p], sigma = output[p1], lambda = lambda0,
            y = Y[, 1], delta = Delta, X_bar = X, y_est = y_est, rord = ordresidual,
            rdev = dev, deviance = devian,  vcov = covar, scores = scores,
             llglg = llglg, llgg = llgg, AIC = aic, BIC = bic,
            AIC2 = aic2, goodnessoffit = good_fit,modelsurv=msurv,st_error = ste, z_values = zs, p.values = pval,
            interval = inter, convergence = conv, condition = condition,
            Iterations = iter, semi = FALSE, censored = TRUE) #
        class(output) = "sglg"
        return(output)
    }

    if (output$conv == FALSE) {
        print('The optimization was not successful.')
        return(0)
    }
}
