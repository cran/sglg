#'Deviance Residuals for a Generalized Log-gamma Regression Model
#'
#'\code{deviance_residuals} is used to generate deviance residuals for a generalized log-gamma regression model.
#'
#' @param object an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg().
#' @param ... other arguments.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' # Example 1
#'n <- 300
#'error <- rglg(n,0,1,1)
#'y <- 0.5 + error
#'fit <- glg(y~1,data=as.data.frame(y))
#'deviance_residuals(fit)
#'# Example 2
#'n <- 300
#'error <- rglg(n,0,1,1)
#'x <- runif(n,-3,3)
#'y <- 0.5 +  2*x + error
#'fit <- glg(y~x,data=as.data.frame(y,x))
#'deviance_residuals(fit)
#'
#'@export deviance_residuals
deviance_residuals <- function(object, ...) {
  lambda <- object$lambda
  rord <- object$rord
  rdev <- object$rdev
  y_est <- object$y_est
  if (object$censored == FALSE) {
    ext <- max(abs(rdev)) + 0.5
    plot1 <- ggplot(data=as.data.frame(rdev),aes(rdev)) + ggtitle("Density Deviance Residuals") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Sample Deviances") + ylab("Density") + geom_hline(yintercept=0)
    plot2 <- ggplot(data=as.data.frame(y_est,rdev),aes(y_est,rdev)) +  ggtitle("Deviance Residuals") + geom_point(colour="blue",alpha=0.5) + xlab("y_i estimated") + ylab("Deviance values") + ylim(c(-ext,ext)) + geom_hline(yintercept=ext) + geom_hline(yintercept=-ext)
    grid.arrange(plot1, plot2, ncol=2)
  }

  if (object$censored == TRUE) {

    delta <- object$delta

    C <- rep(3, length(y_est))
    plot(y_est, rdev, main = "Deviance residuals", xlab = "Fitted values",
         ylab = "Deviance-type residuals", ylim = c(-3.2, 3.2), pch = 20)
    abline(h = C, col = 2)
    abline(h = -C, col = 2)

    ekm <- survfit(Surv(exp(rord), 1 - delta) ~ 1)
    surv <- as.numeric(unlist(as.vector(summary(ekm)[6])))
    Fkm <- 1 - surv
    # plot(ftimes, surv, xlab = 'Multiplicative error', ylab = 'Survival
    # values', main = 'Survival function', type = 'l', pch = 20)

    res <- sort((rord * (1 - delta))[delta == 0])
    Fs <- pglg(res, shape = lambda)
    r_q <- qnorm(Fs)
    diff <- abs(r_q - qnorm(Fkm))
    output <- mean(diff[-length(diff)])
    qqnorm(r_q, xlab = "Quantiles of N(0,1)", ylab = "Overall residuals",
           main = "Overall goodness-of-fit", pch = 20)
    qqline(r_q, col = 2)
  }
}
