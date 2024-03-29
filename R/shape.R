#' shape
#'
#' Tool that supports the estimation of the shape parameter in semi-parametric or multiple linear accelerated failure time model with generalized log-gamma errors
#' under the presence of censored data. The estimation is based on the profiled likelihood function for the shape parameter of the model.

#' @param formula a symbolic description of the systematic component of the model to be fitted.
#' @param data a data frame which contains the variables in the model.
#' @param npc a data frame with potential nonparametric variables of the systematic part of the model to be fitted.
#' @param semi a logical value. TRUE means that the model has a non-parametric component. By default is FALSE.
#' @param interval an optional numerical vector of length 2. In this interval is the maximum likelihood estimate of the shape parameter of the model.
#' By default is [0.1,1.5].
#' @param step an optional positive value. This parameter represents the length of the step of the partition of the interval parameter.
#' By default is 0.1.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' rows  <- 200
#' columns <- 2
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#' t_lambda <- 1
#' set.seed(8142031)
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(columns, 0, 1)
#' X <- cbind(x1,x2)
#' s         <- t_sigma^2
#' a         <- 1/s
#' t_ini1    <- exp(X %*% t_beta) * rgamma(rows, scale = s, shape = a)
#' cens.time <- rweibull(rows, 0.3, 14)
#' delta     <- ifelse(t_ini1 > cens.time, 1, 0)
#' obst1 = t_ini1
#' for (i in 1:rows) {
#' if (delta[i] == 1) {
#'    obst1[i] = cens.time[i]
#'   }
#' }
#' example <- data.frame(obst1,delta,X)
#' lambda <- shape(Surv(log(obst1),delta) ~ x1 + x2 - 1, data=example)
#' lambda
#' # To change interval or step or both options
#' lambda <- shape(Surv(log(obst1),delta) ~ x1 + x2 - 1, data=example, interval=c(0.95,1.3), step=0.05)
#' lambda
#' @export shape

shape = function(formula, npc, data, interval, semi, step){
    if (missingArg(interval))
        interval <- c(0.1, 1.5)
    if (missingArg(step))
        step <- 0.1
    if (missingArg(semi))
        semi = FALSE
    if (semi == TRUE & missingArg(npc))
        stop("If the model is semiparametric you need to enter the non-parametric variable in the argument npc!")

    sh <- seq(interval[1], interval[2], by = step)
    output <- 0 * sh

    if (semi == FALSE) {
        for (j in 1:length(sh)) {
            out <- try(survglg(formula, data = data, shape = sh[j])$llglg,
                silent = TRUE)
            if (is.numeric(out) == FALSE) {
                out <- NA
            }
            output[j] <- out
        }
    }
    if (semi == TRUE) {
        for (j in 1:length(sh)) {
            out <- try(ssurvglg(formula, data = data, npc = npc, shape = sh[j])$llglg,
                silent = TRUE)
            if (is.matrix(out) == FALSE) {
                out <- NA
            }
            output[j] <- out
        }
    }
    index <- which.max(output)
    plot1 <- ggplot(data=as.data.frame(cbind(sh,output)), aes(sh,output)) +
    geom_point(colour="orange",alpha=1,size=1.25) +
    xlim(range(sh)) +
    geom_point(data=as.data.frame(cbind(sh[index],output[index])),aes(sh[index],output[index]), colour="blue", alpha=0.3, size=3.5) +
    labs(x = "Shape parameter", y = "Log-lik") +
    ggtitle("Profile Log-likelihood")
    suppressWarnings(grid.arrange(plot1,ncol=1))
    return(sh[index])
}

