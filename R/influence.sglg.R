#' influence
#'
#' influence.sglg extracts from a object of class sglg the local influence measures and displays their graphs versus the index of the observations.
#'
#' @param model an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg().
#' @param ... other arguments.
#' @references  Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.

influence.sglg <- function(model, ...) {
    lambda <- model$lambda
    sigma <- model$sigma
    epsil <- model$rord
    if (model$censored == FALSE) {
        X <- model$X
        n <- dim(X)[1]
        p <- dim(X)[2]
        
        cweight <- function(model) {
            
            Delt_bw <- matrix(0, p, n)
            for (i in 1:p) {
                Delt_bw[i, ] <- X[, i] * (-(1/(sigma * lambda)) * 
                  (1 - exp(lambda * epsil)))
            }
            
            Delt_sw <- -(1/sigma) - (1/sigma * lambda) * 
                epsil + (1/sigma * lambda) * epsil * exp(lambda * 
                epsil)
            
            Delt_lw <- (1/lambda) + (2/lambda^3) * digamma(1/lambda^2) - 
                (2/lambda^3) + 2 * log(lambda^2)/lambda^3 - 
                (1/lambda^2) * epsil + (2/lambda^3) * exp(lambda * 
                epsil) - (1/lambda^2) * epsil * exp(lambda * 
                epsil)
            
            Delta <- Delt_bw
            Delta <- rbind(Delta, t(Delt_sw))
            Delta <- rbind(Delta, t(Delt_lw))
            
            NC = t(Delta) %*% model$Itheta %*% Delta
            norm = sqrt(sum(diag(t(NC) %*% NC)))
            CNC = NC/norm
            Eigen = eigen(CNC)
            Cmax <- Eigen$vectors[, 1]
            plot(1:n, abs(Cmax), xlab = "Index", ylab = "Local influence", 
                main = "Case-weight perturbation", pch = 20)
            dCNC = diag(CNC)
            plot(1:n, dCNC, xlab = "Index", ylab = "Total local influence", 
                main = "Case-weight perturbation", pch = 20)
        }
        
        respert <- function(model) {
            
            D = function(epsilon, lambd) {
                w <- as.vector(exp(lambd * epsilon))
                D_eps <- diag(w, n, n)
                return(D_eps)
            }
            Ds <- D(epsil, lambda)
            
            Delt_bw <- (1/sigma^2) * t(X) %*% Ds
            Delt_sw <- (1/sigma * lambda^2) * (Ds %*% (lambda * 
                epsil + 1) - 1)
            Delt_lw <- (1/sigma * (lambda^2)) * (-1 + Ds %*% 
                (1 - lambda * epsil))
            
            Delta <- Delt_bw
            Delta <- rbind(Delta, t(Delt_sw))
            Delta <- rbind(Delta, t(Delt_lw))
            
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
        
        par(mfrow = c(2, 2))
        cweight(model)
        respert(model)
        
    }
    if (model$censored == TRUE) {
        
        delta <- model$delta
        cweight <- function(model) {
            delta <- model$delta
            
            X_bar <- model$X_bar
            n <- model$size
            qs <- dim(X_bar)[2]
            
            aaa <- (1/lambda^2) * exp(lambda * epsil)
            S <- pgamma((1/lambda^2) * exp(lambda * epsil), 
                1/lambda^2, lower.tail = FALSE)
            bb <- (aaa^(1/lambda^2) * exp(-aaa))/(gamma(1/lambda^2) * 
                S)
            
            
            Delt_gsw <- matrix(0, qs, n)
            for (j in 1:qs) {
                Delt_gsw[j, ] = X_bar[, j] * ((1 - delta) * 
                  (1/(lambda * sigma)) * (exp(lambda * epsil) - 
                  1) + delta * (lambda/sigma) * bb)
            }
            
            part1 <- (1/sigma) * (-1 - (1/lambda) * (epsil - 
                epsil * exp(lambda * epsil)))
            part2 <- (lambda/sigma) * epsil * bb
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
            plot(1:n, lci, xlab = "Index", ylab = "e_max_i", 
                main = "", pch = 20)
            pinfp <- order(lci)[(n - 1):n]
            text(pinfp, lci[pinfp], label = as.character(pinfp), 
                cex = 0.7, pos = 4)
            dCNC = diag(CNC)
            plot(1:n, dCNC, xlab = "Index", ylab = "B_i", 
                main = "", pch = 20)
            pinfp <- order(dCNC)[(n - 1):n]
            text(pinfp, dCNC[pinfp], label = as.character(pinfp), 
                cex = 0.7, pos = 4)
        }
        
        par(mfrow = c(1, 2))
        cweight(model)
        
    }
}
