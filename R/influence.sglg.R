#' influence
#'
#' influence.sglg extracts from a object of class sglg the local influence measures and displays their graphs versus the index of the observations.
#'
#' @param model an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg().
#' @param ... other arguments.
#' @references  Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' rows <- 100
#' columns <- 2
#' t_beta  <- c(0.5, 2)
#' t_sigma <- 1
#' t_lambda <- 1
#' set.seed(8142031)
#' x1 <- rbinom(rows, 1, 0.5)
#' x2 <- runif(columns, 0, 1)
#' X <- cbind(x1,x2)
#' error <- rglg(rows, 0, 1, t_lambda)
#' y1 <- X %*%t_beta + t_sigma * error
#' data.example <- data.frame(y1,X)
#' fit1 <- glg(y1 ~ x1 + x2 - 1,data=data.example)
#' influence(fit1)
#'
#' @export
influence.sglg <- function(model, ...) {
    epsil <- model$rord
    if (model$censored == FALSE) {

        if (model$Knot == 0)
            Xbar <- model$X
        if (model$Knot > 0)
            Xbar <- model$N

        n <- dim(Xbar)[1]
        q <- dim(Xbar)[2]

        plot1 <- ggplot()
        plot2 <- ggplot()
        plot3 <- ggplot()
        plot4 <- ggplot()

        cweight <- function(model) {

            Delt_gammas <- matrix(0, q, n)
            for (i in 1:q) {
                Delt_gammas[i, ] <- Xbar[, i] * (-(1/(model$sigma * model$lambda)) *
                  (1 - exp(model$lambda * epsil)))
            }

            Delt_sw <- -(1/model$sigma) - (1/model$sigma * model$lambda) *
                epsil + (1/model$sigma * model$lambda) * epsil * exp(model$lambda *
                epsil)

            Delt_lw <- (1/model$lambda) + (2/model$lambda^3) * digamma(1/model$lambda^2) -
                (2/model$lambda^3) + 2 * log(model$lambda^2)/model$lambda^3 -
                (1/model$lambda^2) * epsil + (2/model$lambda^3) * exp(model$lambda *
                epsil) - (1/model$lambda^2) * epsil * exp(model$lambda *
                epsil)

            Delta <- Delt_gammas
            Delta <- rbind(Delta, t(Delt_sw))
            Delta <- rbind(Delta, t(Delt_lw))

            NC = t(Delta) %*% model$Itheta %*% Delta
            norm = sqrt(sum(diag(t(NC) %*% NC)))
            CNC = NC/norm
            Eigen = eigen(CNC)
            Cmax <- Eigen$vectors[, 1]

            dCNC = diag(CNC)

            vls1 <- abs(Cmax)
            df1 <- as.data.frame(vls1)
            plot1 <- ggplot(data=df1,aes(1:n,vls1)) +
            geom_point(colour="orange",alpha=0.5) +
            ggtitle("Case-weight perturbation") +
            xlab("Index") +
            ylab("Local Influence")

            vls2 <- dCNC
            df2 <- as.data.frame(vls2)
            plot2 <- ggplot(data=df2,aes(1:n,vls2)) +
            geom_point(colour="orange",alpha=0.5) +
            ggtitle("Case-weight perturbation") +
            xlab("Index") +
            ylab("Total Local Influence")

            return(list(plot1=plot1,plot2=plot2))
            }

        respert <- function(model) {

            D = function(epsilon, lambd) {
                w <- as.vector(exp(lambd * epsilon))
                D_eps <- diag(w, n, n)
                return(D_eps)
            }
            Ds <- D(epsil, model$lambda)

            Delt_gammas <- (1/model$sigma^2) * t(Xbar) %*% Ds
            Delt_sw <- (1/model$sigma * model$lambda^2) * (Ds %*% (model$lambda *
                epsil + 1) - 1)
            Delt_lw <- (1/model$sigma * (model$lambda^2)) * (-1 + Ds %*%
                (1 - model$lambda * epsil))

            Delta <- Delt_gammas
            Delta <- rbind(Delta, t(Delt_sw))
            Delta <- rbind(Delta, t(Delt_lw))

            NC = t(Delta) %*% model$Itheta %*% Delta
            norm = sqrt(sum(diag(t(NC) %*% NC)))
            CNC = NC/norm
            Eigen = eigen(CNC)
            Cmax = Eigen$vectors[, 1]
            dCNC = diag(CNC)

            vls3 <- abs(Cmax)
            df3 <- as.data.frame(vls3)
            plot3 <- ggplot(data=df3,aes(1:n,vls3)) +
              geom_point(colour="orange",alpha=0.5) +
              ggtitle("Response Perturbation") +
              xlab("Index") +
              ylab("Local Influence")

            vls4 <- dCNC
            df4 <- as.data.frame(vls4)
            plot4 <- ggplot(data=df4,aes(1:n,vls4)) +
              geom_point(colour="orange",alpha=0.5) +
              ggtitle("Response Perturbation") +
              xlab("Index") +
              ylab("Total Local Influence")

            return(list(plot3=plot3,plot4=plot4))

        }

        #par(mfrow = c(2, 2))
        plots_cw <- cweight(model)
        plots_r <- respert(model)
        grid.arrange(plots_cw$plot1, plots_cw$plot2, plots_r$plot3, plots_r$plot4, ncol=2)

    }
    if (model$censored == TRUE) {

        delta <- model$delta
        cweight <- function(model) {

            X_bar <- model$X_bar
            n <- model$size
            qs <- dim(X_bar)[2]

            aaa <- (1/model$lambda^2) * exp(model$lambda * epsil)
            S <- pgamma((1/model$lambda^2) * exp(model$lambda * epsil), 1/model$lambda^2,
                lower.tail = FALSE)
            bb <- (aaa^(1/model$lambda^2) * exp(-aaa))/(gamma(1/model$lambda^2) *
                S)


            Delt_gsw <- matrix(0, qs, n)
            for (j in 1:qs) {
                Delt_gsw[j, ] = X_bar[, j] * ((1 - delta) * (1/(model$lambda *
                  model$sigma)) * (exp(model$lambda * epsil) - 1) + delta *
                  (model$lambda/model$sigma) * bb)
            }

            part1 <- (1/model$sigma) * (-1 - (1/model$lambda) * (epsil -
                epsil * exp(model$lambda * epsil)))
            part2 <- (model$lambda/model$sigma) * epsil * bb
            Delt_sw <- (1 - delta) * part1 + delta * part2
            Delt_sw <- t(as.matrix(Delt_sw))

            Delta <- Delt_gsw
            Delta <- rbind(Delta, Delt_sw)

            NC = t(Delta) %*% model$Itheta %*% Delta
            norm = sqrt(sum(diag(t(NC) %*% NC)))
            CNC = NC/norm
            Eigen = eigen(CNC)
            Cmax <- Eigen$vectors[, 1]
            lci <- abs(Cmax)
            plot(1:n, lci, xlab = "Index", ylab = "e_max_i", main = "Case-weight perturbation",
                pch = 20)
            pinfp <- order(lci)[(n - 1):n]
            text(pinfp, lci[pinfp], label = as.character(pinfp), cex = 0.7,
                pos = 4)
            dCNC = diag(CNC)
            plot(1:n, dCNC, xlab = "Index", ylab = "B_i", main = "Case-weight perturbation",
                pch = 20)
            pinfp <- order(dCNC)[(n - 1):n]
            text(pinfp, dCNC[pinfp], label = as.character(pinfp), cex = 0.7,
                pos = 4)
        }

        respert <- function(model) {
            X <- model$X
            n <- dim(X)[1]
            p <- dim(X)[2]

            sigma <- model$sigma
            lambda <- model$lambda

            w <- as.vector(exp(lambda * epsil))
            Ds <- diag(w, n)
            aaa = (1/lambda^2) * exp(lambda * epsil)
            S = pgamma((1/lambda^2) * exp(lambda * epsil), 1/lambda^2, lower.tail = FALSE)
            bb = (aaa^(1/lambda^2) * exp(-aaa))/(gamma(1/lambda^2) * S)

            Delt_bw <- matrix(0, p, n)
            for (j in 1:p) {
                Delt_bw[j, ] = X[, j] * ((1 - delta) * (1/sigma^2) * exp(lambda *
                  epsil) + delta * (-(lambda/sigma)^2) * bb * (aaa - 1/lambda^2 -
                  bb))
            }

            Delt_sw <- (1/sigma * lambda^2) * (Ds %*% (lambda * epsil + 1) -
                1)
            expep <- exp(lambda * epsil)
            Delt_sw <- (1 - delta) * (1/(lambda * (sigma^2))) * (-1 + expep +
                lambda * epsil * expep) + delta * ((((lambda/sigma)^2) *
                bb) * (-1/lambda + epsil * (aaa - 1/lambda^2 - bb)))
            Delt_sw <- t(as.matrix(Delt_sw))
            Delta <- Delt_bw
            Delta <- rbind(Delta, Delt_sw)
            NC = t(Delta) %*% model$Itheta %*% Delta
            norm = sqrt(sum(diag(t(NC) %*% NC)))
            CNC = NC/norm
            Eigen = eigen(CNC)
            Cmax = Eigen$vectors[, 1]
            plot(1:n, abs(Cmax), xlab = "Index", ylab = "Local influence",
                main = "Response perturbation", pch = 20)
            dCNC = diag(CNC)
            plot(1:n, dCNC, xlab = "Index", ylab = "Total local influence",
                main = "Response perturbation", pch = 20)
        }

    }
}
