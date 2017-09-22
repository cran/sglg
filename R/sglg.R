#'Fitting semi-parametric generalized log-gamma regression models
#'
#'\code{sglg} is used to fit a semi-parametric regression model suitable for analysis of data sets in which the response variable is continuous, strictly positive, and asymmetric.
#'In this setup, the location parameter of the response variable is explicitly modeled by semi-parametric functions, whose nonparametric components may be approximated by
#'natural cubic splines or cubic P-splines.
#'
#' @param formula a symbolic description of the systematic component of the model to be fitted. See details for further information.
#' @param npc a data frame with potential nonparametric variables of the systematic part of the model to be fitted.
#' @param basis a name of the cubic spline basis to be used in the model. Supported basis include deBoor and Gu basis
#'  which are a B-spline basis and a natural cubic spline basis, respectively.
#' @param data an optional data frame, list containing the variables in the model.
#' @param shape an optional value for the shape parameter of the error distribution of a generalized log-gamma distribution. Default value is 1.
#' @param Tolerance an optional positive value, which represents the convergence criterion. Default value is 1e-04.
#' @param Maxiter an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1e03.

#' @return mu a vector of parameter estimates asociated with the location parameter.
#' @return sigma estimate of the scale parameter associated with the model.
#' @return lambda estimate of the shape parameter associated with the model.
#' @return interval estimate of a 95\% confidence interval for each estimate parameters associated with the model.
#' @return Deviance the deviance associated with the model.

#' @references Carlos A. Cardozo, G. Paula and L. Vanegas. Semi-parametric generalized log-gamma regression models. In preparation.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' rows <- 240 # Number of observations
#' columns <- 2 # Number of parametric components

#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#' t_lambda <- 1

#' set.seed(8142031)
#' library(ssym)
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(columns, 0, 1)
#' X <- cbind(x1,x2)
#' t_knot1 <- 6
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
#' error <- robustloggamma::rloggamma(rows, 0, 1, t_lambda)
#' y2 <- X %*%t_beta + + s_N1 %*% t_g1 + t_sigma * error
#' data.example <- data.frame(y2,X,x3)
#' fit2 <- sglg(y2 ~ x1 + x2, npc=x3, data=data.example)

#' @import ssym
#' @import robustloggamma
#' @import methods
#' @export sglg

sglg = function(formula, npc, basis, data, shape, Tolerance, 
    Maxiter) {
    if (missingArg(formula)) {
        stop("The formula argument is missing.")
    }
    if (missingArg(npc)) {
        stop("This kind of model need at least one non-parametric component.")
    }
    if (missingArg(data)) {
        stop("The data argument is missing.")
    }
    
    if (missingArg(shape)) 
        shape <- 1
    if (missingArg(Tolerance)) 
        Tolerance <- 1e-04
    if (missingArg(Maxiter)) 
        Maxiter <- 1000
    
    if (missingArg(basis)) 
        basis <- rep("deBoor", dim(npc)[2])
    
    if (class(data) == "list") 
        data <- as.data.frame(data)
    
    data1 <- model.frame(formula, data = data)
    X <- model.matrix(formula, data = data1)
    y <- model.response(data1)
    
    p <- ncol(X)
    n <- nrow(X)
    
    intknt = function(x) {
        op1 <- floor(n^(1/3)) + 3
        op2 <- length(as.numeric(levels(factor(x))))
        knt <- min(op1, op2)
        if (knt < 3) {
            stop("This covariate has not at least three different values.")
        }
        return(knt)
    }
    
    k <- dim(npc)[2]
    XX <- cbind(X, npc)
    
    Knot <- matrix(0, k, 1)
    for (i in 1:k) {
        Knot[i] <- intknt(XX[, (p + i)])
    }
    Knot <- as.numeric(Knot)
    
    npoutput <- deBoor(npc, Knot)
    N <- npoutput$N
    K <- npoutput$K
    
    g0 <- function(knot) {
        g <- 0 * 1:knot
        return(g)
    }
    
    ############################################################################################################ 
    
    formula2 <- paste(".~. + ", colnames(npc))
    formula2 <- as.formula(formula2)
    formula2 <- update(formula, formula2)
    
    ############################################################################################################################################################ 
    
    # Initial values
    
    fit0 <- glg(formula2, shape = shape, data = data)
    beta0 <- fit0$mu[1:p]
    g0s <- g0(Knot)
    sigma0 <- fit0$sigma
    lambda0 <- fit0$lambda
    
    formula3 <- paste("ncs(", colnames(npc), sep = "")
    formula3 <- paste(formula3, ")", sep = "")
    formula3 <- paste(".~. + ", formula3, sep = "")
    formula3 <- as.formula(formula3)
    formula3 <- update(formula, formula3)
    
    fit00 <- ssym::ssym.l(formula3, data = data, family = "Normal")
    alpha0 <- fit00$lambdas.mu
    
    Ident <- diag(1, n)
    Ones <- matrix(1, n, 1)
    
    tNN = function(N) {
        output = t(N) %*% N
        return(output)
    }
    
    I_1 <- tNN(X)
    
    I_s <- matrix(0, Knot, Knot)
    
    for (i in 1:k) {
        if (i > 1) {
            l = 1 + sum(Knot[1:(i - 1)])
            r = sum(Knot[1:i])
        }
        if (i == 1) {
            l = 1
            r = Knot[1]
        }
        I_s[l:r, l:r] = tNN(N[, l:r])
    }
    
    ## THE FISHER INFORMATION MATRIX
    
    I_beta = function(sigm) {
        output = (1/(sigm^2)) * I_1
        return(output)
    }
    
    I_g = function(i, sigm, alph) {
        if (i > 1) {
            l = 1 + sum(Knot[1:(i - 1)])
            r = sum(Knot[1:i])
        }
        if (i == 1) {
            l = 1
            r = Knot[i]
        }
        output = (1/(sigm^2)) * I_s[l:r, l:r] + alph[i] * 
            K[l:r, l:r]
        return(output)
    }
    
    I_betg = function(i, sigm) {
        if (i > 1) {
            l = 1 + sum(Knot[1:(i - 1)])
            r = sum(Knot[1:i])
        }
        if (i == 1) {
            l = 1
            r = Knot[1]
        }
        output = (1/(sigm^2)) * t(X) %*% N[, l:r]
        return(output)
    }
    
    I_betl = function(sigm, lambd) {
        output = -(1/(lambd * sigm)) * u_lambda(lambd) * 
            t(X) %*% Ones
        return(output)
    }
    
    u_lambda = function(lambd) {
        invlamb = 1/lambd^2
        output = (1/lambd) * (digamma(1 + invlamb) - log(invlamb))
        return(output)
    }
    
    I_betsig = function(sigm, lambd) {
        output = (1/(sigm^2)) * u_lambda(lambd) * t(X) %*% 
            Ones
        return(output)
    }
    
    I_gsig = function(i, sigm, lambd) {
        if (i > 1) {
            l = 1 + sum(Knot[1:(i - 1)])
            r = sum(Knot[1:i])
        }
        if (i == 1) {
            l = 1
            r = Knot[1]
        }
        output = (1/(sigm^2)) * u_lambda(lambd) * t(N[, l:r]) %*% 
            Ones
        return(output)
    }
    
    I_gl = function(i, sigm, lambd) {
        if (i > 1) {
            l = 1 + sum(Knot[1:(i - 1)])
            r = sum(Knot[1:i])
        }
        if (i == 1) {
            l = 1
            r = Knot[1]
        }
        output = -(1/(lambd * sigm)) * u_lambda(lambd) * 
            t(N[, l:r]) %*% Ones
        return(output)
    }
    
    ## Weight matrix for the score functions
    
    ### Defining mu function
    
    mu = function(bet, g) {
        X_bar <- cbind(X, N)
        bet <- as.matrix(bet)
        g <- as.matrix(g)
        g_bar <- rbind(bet, g)
        output = X_bar %*% g_bar
        return(output)
    }
    
    ### First step: The residuals
    
    eps = function(bet, g, sigm) {
        epsilon = (y - mu(bet, g))/sigm
        return(epsilon)
    }
    
    ### Second step: The matrix D
    
    D = function(epsil, lambd) {
        w = as.vector(exp(lambd * epsil))
        D_s = diag(w, n, n)
        return(D_s)
    }
    
    ### Third step: The matrix W
    
    W = function(bet, g, sigm, lambd) {
        output = Ident - D(eps(bet, g, sigm), lambd)
        return(output)
    }
    
    ## The score functions
    
    U_beta = function(bet, g, sigm, lambd) {
        output = (-1/(lambd * sigm)) * t(X) %*% W(bet, g, 
            sigm, lambd) %*% Ones
        return(output)
    }
    
    U_g = function(bet, g, i, sigm, lambd, alph) {
        if (i > 1) {
            l = 1 + sum(Knot[1:(i - 1)])
            r = sum(Knot[1:i])
        }
        if (i == 1) {
            l = 1
            r = Knot[1]
        }
        output = (-1/(lambd * sigm)) * t(N[, l:r]) %*% W(bet, 
            g, sigm, lambd) %*% Ones - alph[i] * K[l:r, l:r] %*% 
            g[l:r]
        return(output)
    }
    
    U_sigma = function(bet, g, sigm, lambd) {
        output = -(1/sigm) * n - (1/(lambd * sigm)) * t(Ones) %*% 
            W(bet, g, sigm, lambd) %*% eps(bet, g, sigm)
        return(output)
    }
    
    U_lambda = function(bet, g, sigm, lambd) {
        invlamb = 1/lambd
        eta_lambd = (invlamb) * (1 + 2 * (invlamb^2) * (digamma(invlamb^2) + 
            2 * log(abs(lambd)) - 1))
        Ds = D(eps(bet, g, sigm), lambd)
        epsilons = eps(bet, g, sigm)
        output = n * eta_lambd - (invlamb^2) * t(Ones) %*% 
            epsilons + (2 * invlamb^3) * t(Ones) %*% Ds %*% 
            Ones - (invlamb^2) * t(Ones) %*% Ds %*% epsilons
        return(output)
    }
    
    U_sl = function(bet, g, sigm, lambd) {
        output = matrix(1, 2, 1)
        output[1] = U_sigma(bet, g, sigm, lambd)
        output[2] = U_lambda(bet, g, sigm, lambd)
        return(output)
    }
    
    U_tetha = function(bet, g, sigm, lambd, alph) {
        output = matrix(0, p + Knot + 2, 1)
        output[1:p] = U_beta(bet, g, sigm, lambd)
        for (i in 1:k) {
            if (i > 1) {
                l = 1 + sum(Knot[1:(i - 1)])
                r = sum(Knot[1:i])
            }
            if (i == 1) {
                l = 1
                r = Knot[i]
            }
            output[(p + l):(p + r)] = U_g(bet, g, i, sigm, 
                lambd, alph)
        }
        output[(p + Knot + 1):(p + Knot + 2)] = U_sl(bet, 
            g, sigm, lambd)
        return(output)
    }
    
    ## LOG-LIKELIHOOD
    
    c_l = function(lambd) {
        invlambdos = 1/lambd^2
        c = abs(lambd)/gamma(invlambdos)
        output = c * (invlambdos^invlambdos)
        return(output)
    }
    
    loglikglg = function(bet, g, sigm, lambd, alph) {
        epsilon = eps(bet, g, sigm)
        part1 = n * log(c_l(lambd)/sigm) + (1/lambd) * t(Ones) %*% 
            epsilon - (1/lambd^2) * t(Ones) %*% D(epsilon, 
            lambd) %*% Ones
        part2 = 0
        for (i in 1:k) {
            if (i > 1) {
                l = 1 + sum(Knot[1:(i - 1)])
                r = sum(Knot[1:i])
            }
            if (i == 1) {
                l = 1
                r = Knot[1]
            }
            part2 = part2 - (alph[i]/2) * g[l:r] %*% K[l:r, 
                l:r] %*% g[l:r]
        }
        output = part1 + part2
        return(output)
    }
    
    ## Estimating sigma and lambda
    
    v_lambda = function(lambd) {
        invlamb = 1/lambd^2
        output = invlamb * trigamma(1 + invlamb) + u_lambda(lambd)^2
        return(output)
    }
    
    # K_1_lambda and K_2_lambda functions
    
    K_1 = function(lambda) {
        invlamb2 = 1/lambda^2
        part1 = 4 * (1 + digamma(1 + invlamb2) - digamma(invlamb2) - 
            invlamb2 * trigamma(invlamb2))
        part2 = trigamma(1 + invlamb2) + (digamma(1 + invlamb2) - 
            log(invlamb2))^2
        output = 1 - invlamb2 * (part1 - part2)
        return(output)
    }
    
    K_2 = function(lambda) {
        invlamb = 1/lambda
        invlamb2 = invlamb^2
        output = invlamb * ((digamma(1 + invlamb2) - digamma(invlamb2)) - 
            trigamma(1 + invlamb2) - (digamma(1 + invlamb2) - 
            log(invlamb2))^2)
        return(output)
    }
    
    ## Defining the components of the FIM
    
    I_22 = function(sigm, lambd) {
        output = (n/(sigm^2)) * (1 + v_lambda(lambd))
        return(output)
    }
    
    I_23 = function(sigm, lambd) {
        output = (n/(sigm * lambd^2)) * K_2(lambd)
        return(output)
    }
    
    I_33 = function(sigm, lambd) {
        output = (n/(lambd^2)) * K_1(lambd)
        return(output)
    }
    
    I_sl = function(sigm, lambd) {
        output = matrix(0, 2, 2)
        output[1, 1] = I_22(sigm, lambd)
        output[1, 2] = I_23(sigm, lambd)
        output[2, 1] = I_23(sigm, lambd)
        output[2, 2] = I_33(sigm, lambd)
        return(output)
    }
    
    I_tetha = function(sigm, lambd, alph) {
        covar = matrix(0, p + Knot + 2, p + Knot + 2)
        covar[1:p, 1:p] = I_beta(sigm)
        for (i in 1:k) {
            if (i > 1) {
                l = 1 + sum(Knot[1:(i - 1)])
                r = sum(Knot[1:i])
            }
            if (i == 1) {
                l = 1
                r = Knot[i]
            }
            covar[1:p, (p + l):(p + r)] = I_betg(i, sigm)
            covar[(p + l):(p + r), 1:p] = t(covar[1:p, (p + 
                l):(p + r)])
            covar[(p + l):(p + r), (p + l):(p + r)] = I_g(i, 
                sigm, alph)
            covar[(p + l):(p + r), (p + Knot + 1)] = I_gsig(i, 
                sigm, lambd)
            covar[(p + Knot + 1), (p + l):(p + r)] = t(covar[(p + 
                l):(p + r), (p + Knot + 1)])
            covar[(p + l):(p + r), (p + Knot + 2)] = I_gl(i, 
                sigm, lambd)
            covar[(p + Knot + 2), (p + l):(p + r)] = t(covar[(p + 
                l):(p + r), (p + Knot + 2)])
            
        }
        covar[1:p, (p + Knot + 1)] = I_betsig(sigm, lambd)
        covar[(p + Knot + 1), 1:p] = t(covar[1:p, (p + Knot + 
            1)])
        covar[1:p, (p + Knot + 2)] = I_betl(sigm, lambd)
        covar[(p + Knot + 2), 1:p] = t(covar[1:p, (p + Knot + 
            2)])
        covar[(p + Knot + 1):(p + Knot + 2), (p + Knot + 
            1):(p + Knot + 2)] = I_sl(sigm, lambd)
        return(covar)
    }
    
    newpar = function(bet, g, sigm, lambd, alph) {
        output = matrix(0, p + Knot + 2, Maxiter)
        output[, 1] = c(bet, g, sigm, lambd)
        new = output[, 1]
        llglg = loglikglg(new[1:p], new[(p + 1):(p + Knot)], 
            new[p + Knot + 1], new[p + Knot + 2], alph)
        output[, 2] = output[, 1] + solve(I_tetha(sigm, lambd, 
            alph)) %*% U_tetha(bet, g, sigm, lambd, alph)
        condition = loglikglg(output[1:p, 2], output[(p + 
            1):(p + Knot), 2], output[p + Knot + 1, 2], output[p + 
            Knot + 2, 2], alph) - llglg
        if (condition > 0) {
            new = output[, 2]
        }
        l = 2
        while (condition > Tolerance & l < Maxiter) {
            l = l + 1
            output[, l] = output[, (l - 1)] + solve(I_tetha(output[(p + 
                Knot + 1), (l - 1)], output[(p + Knot + 2), 
                (l - 1)], alph)) %*% U_tetha(output[1:p, 
                (l - 1)], output[(p + 1):(p + Knot), (l - 
                1)], output[(p + Knot + 1), (l - 1)], output[(p + 
                Knot + 2), (l - 1)], alph)
            llglg = loglikglg(new[1:p], new[(p + 1):(p + 
                Knot)], new[p + Knot + 1], new[p + Knot + 
                2], alph)
            condition = loglikglg(output[1:p, l], output[(p + 
                1):(p + Knot), l], output[p + Knot + 1, l], 
                output[p + Knot + 2, l], alph) - llglg
            if (condition > 0) {
                new = output[, l]
                llglg = loglikglg(new[1:p], new[(p + 1):(p + 
                  Knot)], new[p + Knot + 1], new[p + Knot + 
                  2], alph)
            }
        }
        
        return(list(est = new, ll = llglg, cond = condition, 
            iter = l))
    }
    
    ## Effective degree freedom - EDF
    
    irQ = function(Q) {
        e = eigen(Q)
        V = e$vectors
        Q = solve(V %*% sqrt(diag(e$values)) %*% t(V))
        output = solve(Q)
        return(output)
    }
    
    Xbar = cbind(X, N)
    Malph = function(alph) {
        output = matrix(0, p + Knot, p + Knot)
        output[(p + 1):(p + Knot), (p + 1):(p + Knot)] = alph * 
            K
        return(output)
    }
    
    edf = function(sigm, alph) {
        XX = (1/sigm^2) * t(Xbar) %*% Xbar
        irXbar = irQ(XX)
        Malpha = Malph(alph)
        output = diag(1, p + Knot) + (sigm^2) * irXbar %*% 
            Malpha %*% irXbar
        output = sum(diag(solve(output)))
        return(output)
    }
    
    Conv = FALSE
    num.iter = 1
    masterf = function(alph) {
        news = newpar(beta0, g0s, sigma0, lambda0, alph)
        if (news$iter < Maxiter) {
            Conv = TRUE
            num.iter = news$iter
            cond = news$cond
            llglg = news$ll
        }
        
        df = edf(news$est[(p + Knot + 1)], alph)  # joint edf X and N
        
        aic = -2 * llglg + 2 * (df + 2)
        bic = -2 * llglg + log(n) * (df + 2)
        
        return(list(est = news$est, df = df, llglg = llglg, 
            AIC = aic, BIC = bic, Conv = Conv, iter = num.iter, 
            cond = cond))
    }
    
    AIC_p = function(alph) {
        tetha = masterf(alph)
        output = tetha$AIC
        output = round(output, digits = 3)
        return(output)
    }
    
    opt_alph = function(alph) {
        out = optimize(AIC_p, c(0, alph + 2), tol = 0.001)
        return(c(out$minimum, out$objective))
    }
    
    gfit = function(resid, lambd) {
        Fs = ploggamma(resid, lambda = lambd)
        equantil = qnorm(Fs)
        diff <- qqnorm(equantil, plot.it = FALSE)
        output = mean(abs(diff$x - diff$y))
        return(output)
    }
    
    total_optimum = function(start) {
        output0 = opt_alph(start)
        output1 = output0[1:k]
        output2 = output0[k + 1]
        output3 = masterf(output1)
        df <- output3$df
        dfnpc <- df - p
        output = output3$est
        llglg = output3$llglg
        aic = output3$AIC
        bic = output3$BIC
        scores = U_tetha(output[1:p], output[(p + 1):(p + 
            Knot)], output[p + Knot + 1], output[p + Knot + 
            2], output1)
        covar = I_tetha(output[p + Knot + 1], output[p + 
            Knot + 2], output1)
        inter = matrix(0, p + Knot + 2, 2)
        scovar = solve(covar)
        val = diag(scovar)
        if (min(val) > 0) {
            ste = sqrt(val)
            inter[, 1] = as.matrix(output - 1.96 * ste)
            inter[, 2] = as.matrix(output + 1.96 * ste)
            zs = abs(output/ste)
            pval = 1 - (pnorm(zs) - pnorm(-zs))
            pval2 = pval[-((p + 1):(p + Knot))]
            as = output[(p + 1):(p + Knot)]
        }
        
        y_est = X %*% output[1:p]
        for (i in 1:k) {
            if (i > 1) {
                l = 1 + sum(Knot[1:(i - 1)])
                r = sum(Knot[1:i])
            }
            if (i == 1) {
                l = 1
                r = Knot[1]
            }
            y_est = y_est + N[, l:r] %*% output[(p + l):(p + 
                r)]
        }
        ordresidual = eps(output[1:p], output[(p + 1):(p + 
            Knot)], output[p + Knot + 1])
        sgn = sign(y - y_est)
        dev = sgn * sqrt(2) * ((1/output[p + Knot + 2]^2) * 
            exp(output[p + Knot + 2] * ordresidual) - (1/output[p + 
            Knot + 2]) * ordresidual - (1/output[p + Knot + 
            2])^2)^(0.5)
        Devian = sum(dev^2)
        par(mfrow = c(1, 2))
        
        good_fit = gfit(ordresidual, output[p + Knot + 2])
        part2 = ((output[p + Knot + 1])/(output[p + Knot + 
            2])) * (digamma((1/output[p + Knot + 2])^2) - 
            log((1/output[p + Knot + 2])^2))
        y_est2 = y_est + part2
        
        return(list(formula = formula, npc = npc, size = n, 
            mu = output[1:(p + Knot)], sigma = output[p + 
                Knot + 1], lambda = output[p + Knot + 2], 
            y = y, X = X, p = p, N = N, Knot = Knot, rord = ordresidual, 
            rdev = dev, interval = inter, llglg = llglg, 
            AIC = output2, BIC = output3$BIC, scores = scores, 
            Itheta = covar, scovar = scovar, st_error = ste, 
            Z_values = zs, p.values = pval, alpha = output1, 
            d.f.model = df, d.f.npc = dfnpc, deviance = Devian, 
            goodnessoffit = good_fit, convergence = output3$Conv, 
            condition = output3$cond, iterations = output3$iter, 
            semi = TRUE, censored = FALSE, y_est2 = y_est))
    }
    output = total_optimum(alpha0)
    class(output) = "sglg"
    return(output)
}
