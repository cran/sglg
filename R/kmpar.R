#'Plot simultaneously the Kaplan-Meier and parametric estimators of the survival function.
#'
#'\code{plotsurv.sglg} is used to plot simultaneously the Kaplan-Meier and parametric estimators of the survival function.
#'
#' @param fit an object of the class sglg. This object is returned from the call to survglg() or ssurvglg().
#' @references Carlos A. Cardozo, G. Paula and L. Vanegas. Semi-parametric accelerated failure time models with generalized log-gamma erros. In preparation.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' require(survival)
#' n  <- 240
#' columns <- 2
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#' t_lambda <- 1
#' set.seed(8142031)
#' x1 <- rbinom(n, 1, 0.5)
#' x2 <- runif(columns, 0, 1)
#' X <- cbind(x1,x2)
#' s         <- t_sigma^2
#' a         <- 1/s
#' t_ini1    <- exp(X %*% t_beta) * rgamma(n, scale = s, shape = a)
#' cens.time <- rweibull(n, 0.6, 14)
#' delta1     <- ifelse(t_ini1 > cens.time, 1, 0)
#' obst1 <- t_ini1
#' obst1[delta1==1] <- cens.time[delta1==1]
#' data <- data.frame(obst1,delta1,X)
#' lambda <- shape(Surv(log(obst1),delta1) ~ x1 + x2 - 1, data=data, interval=c(0.85,0.95), step=0.01)
#' fit3 <- survglg(Surv(log(obst1),delta1) ~ x1 + x2 - 1, data=data, shape=lambda)
#' plotsurv.sglg(fit3)
#' @import Formula
#' @import survival
#' @import methods
#' @export plotsurv.sglg

plotsurv.sglg = function(fit) {
    Delta <- 1 - fit$delta
    km <- survfit(Surv(exp(fit$y), Delta) ~ 1, conf.type = "log-log")
    obst <- exp(fit$y)[Delta==1]
    ###################################################################################

    a <- sort(exp(fit$y))
    b <- km$surv
    df_km  <- as.data.frame(cbind(a,b))

    c <- sort(obst)
    d <- sort(fit$modelsurv,decreasing=TRUE)
    df_fit <- as.data.frame(cbind(c,d))

    ggplot(df_km, aes(a,b)) +
    geom_line() +
    geom_line(data=df_fit,aes(c,d),colour="orange") +
      labs(x = "y", y = "Probability") +
    ggtitle("Estimated Survival Functions")
}

