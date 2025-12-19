#'Plotting a natural cubic splines or P-splines.
#'
#'\code{plotnpc} displays a graph of a fitted nonparametric effect, either natural cubic spline or P-spline, from an object of class sglg.
#'
#' @param fit an object of the class sglg. This object is returned from the call to sglg() or ssurvglg().
#' @param conf_lev is the confidence level of the asymptotic confidence band. Default value is 0.05.
#' @references Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.
#' @references Wood, S. (2017). Additive generalized models: An R introduction. Chapman and Hall.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @import ggplot2
#' @examples
#' set.seed(1)
#' rows<- 300
#' t_beta <- c(0.5,2)
#' t_sigma <- 0.5
#' t_lambda <- 1
#' x1 <- runif(rows,-3,3)
#' x2 <- rnorm(rows,mean=2.5,sd=0.5)
#' X <- cbind(x1,x2)
#' t <- as.matrix(seq(0.01,0.99,length=rows))
#' colnames(t) <- "t"
#' f_t <- cos(4*pi*t)
#' plot(t,f_t,type='l')
#' error <- rglg(rows,0,1,t_lambda)
#' y <- X %*%t_beta + f_t + t_sigma*error
#' colnames(y) <- "y"
#' data <- data.frame(y,X,t)
#' fit1 <- sglg(y ~ x1 + x2 - 1, npc=t, data=data, basis = "Gu", alpha0=0.001)
#' # The adjusted (black) non-linear component
#' plotnpc(fit1,conf_lev=0.02)
#' @export plotnpc

plotnpc <- function(fit,conf_lev) {
    if (fit$semi == FALSE) {
        stop("Sorry, for this kind of model it is not available this option.")
    }
    if(missingArg(conf_lev))
       conf_lev <- 0.05

    y <- fit$y
    N <- fit$N
    npc <- fit$npc
    add_comp <- paste("f(",colnames(npc),sep="")
    add_comp <- paste(add_comp,")",sep="")
    add_comp <- paste("Estimated Additive Component", add_comp,sep=" ")
    p <- fit$p
    mu <- fit$mu
    f_est <- N[,-(1:p)]%*%mu[-(1:p)]
    if (fit$basis == "deBoor"){
        f_est <- f_est + mu[1]
        Knot <- fit$Knot
        var_gammas <- fit$scovar[(p+1):(p+Knot),(p+1):(p+Knot)]
        var_f_est <- N[,-(1:p)]%*%var_gammas%*%t(N[,-(1:p)])
        st_error_f_est <- sqrt(diag(var_f_est))
        f_est_low <- f_est + qnorm(0.5*conf_lev)*st_error_f_est
        f_est_up <- f_est + qnorm(1 - 0.5*conf_lev)*st_error_f_est

        df <- as.data.frame(cbind(y,N,f_est))
        plot <- ggplot(data=df,aes(npc,y))+
        geom_point(colour="blue",alpha=0.4)+
        geom_line(aes(npc,f_est),size=1.2) +
        xlab(colnames(npc))+
        ggtitle(add_comp)
        return(plot)
    }

    Knot <- fit$Knot
    var_gammas <- fit$scovar[(p+1):(p+Knot),(p+1):(p+Knot)]
    var_f_est <- N[,-(1:p)]%*%var_gammas%*%t(N[,-(1:p)])
    st_error_f_est <- sqrt(diag(var_f_est))
    f_est_low <- f_est + qnorm(0.5*conf_lev)*st_error_f_est
    f_est_up <- f_est + qnorm(1 - 0.5*conf_lev)*st_error_f_est
    df <- as.data.frame(cbind(y,N,f_est))
    plot <- ggplot(data=df,aes(npc,y))+
    geom_line(aes(npc,f_est_up),colour = "orange",size=1.2) +
    geom_line(aes(npc,f_est),size=1.2) +
    geom_line(aes(npc,f_est_low),colour = "orange",size=1.2) +
    xlab(colnames(npc))+
    ggtitle(add_comp)
    return(plot)
}


















