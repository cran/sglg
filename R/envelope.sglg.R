#' envelope.sglg
#'
#'Build a Normal probability plot with simulated envelope for a deviance-type residuals
#'in semiparametric or multiple linear generalized log-gamma regression models.
#'
#' @param fit an object of the class sglg. This object is returned from the call to glg(), sglg().
#' @param Rep a positive integer. This is the number of replications on which to build the simulated envelope. Default is Rep=50.

#' @references  Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @references  Ortega, E., Paula, G. A. and Bolfarine, H. (2008) Deviance residuals in generalized log-gamma regression models with censored observations. Journal of Statistical Computation and Simulation, 78, 747-764.

#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' rows <- 120
#' columns <- 2
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 0.5
#' t_lambda <- 1
#' set.seed(8142031)
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(columns, 0, 1)
#' X <- cbind(x1,x2)
#' error <- rglg(rows, 0, 1, t_lambda)
#' y1 <- X %*%t_beta + t_sigma * error
#' data.example <- data.frame(y1,X)
#' fit <- glg(y1 ~ x1 + x2 - 1,data=data.example)
#' envelope.sglg(fit,Rep=50)
#' @import stats
#' @import graphics
#' @export envelope.sglg
#'
envelope.sglg <- function(fit, Rep) {

    if (fit$censored == FALSE) {

        if (missingArg(Rep))
            Rep <- 30

        formula <- paste("~", as.character(fit$formula)[3])
        formula <- as.formula(paste("y", formula))
        X <- fit$X
        sigma <- fit$sigma
        lambda <- fit$lambda
        n <- fit$size
        rdev <- fit$rdev
        rord <- fit$rord
        systematic_part <- rord * sigma
        e <- matrix(0, n, Rep)

        j <- 1

        if (fit$Knot >= 3) {
            npc <- fit$npc
            while (j <= Rep) {
                error <- rglg(n, shape = lambda)
                y <- systematic_part + sigma * error
                data <- data.frame(y, X)
                newfit <- try(sglg(formula, npc = npc, data = data), silent = TRUE)
                if (is.list(newfit)) {
                  if (newfit$convergence == TRUE) {
                    print(j)
                    e[, j] <- sort(newfit$rdev)
                    j <- j + 1
                  }
                }
            }
        }
        if (fit$Knot == 0) {
            while (j <= Rep) {
                error <- rglg(n, shape = lambda)
                y <- systematic_part + sigma * error
                data <- data.frame(y, X)
                newfit <- try(glg(formula, data = data), silent = TRUE)
                if (is.list(newfit)) {
                  if (newfit$convergence == TRUE) {
                    print(j)
                    e[, j] <- sort(newfit$rdev)
                    j <- j + 1
                  }
                }
            }
        }
        e1 <- numeric(n)
        e2 <- numeric(n)

        for (i in 1:n) {
            eo <- sort(e[i, ])
            e1[i] <- (eo[1] + eo[2])/2
            e2[i] <- (eo[Rep - 1] + eo[Rep])/2
        }

        med <- apply(e, 1, mean)
        faixa <- range(rdev, e1, e2)
        plot1 <- ggplot(data=as.data.frame(cbind(rdev,e1,e2,med)),aes(sample=rdev))+
        stat_qq(colour="orange",alpha=0.5) +
        stat_qq(aes(sample=e1),colour="black",alpha=0.5,geom="line") +
        stat_qq(aes(sample=e2),colour="black",alpha=0.5,geom="line") +
        stat_qq(aes(sample=med),colour="black",alpha=0.5,geom="line") +
        ggtitle("Envelope") +
        xlab("Theoretical Quantiles") +
        ylab("Desviance-type residuals")
        grid.arrange(plot1, ncol=1)
    }

    if (fit$censored == TRUE) {
        print("Sorry, for this kind of model it is not available this option.")
    }
}
