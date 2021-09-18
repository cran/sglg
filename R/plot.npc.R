#'Plotting a natural cubic splines or P-splines.
#'
#'\code{plotnpc} displays a graph of a fitted nonparametric effect, either natural cubic spline or P-spline, from an object of class sglg.
#'
#' @param fit an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg().
#' @param conf_lev is the confidence level of the asymptotic confidence band. Default value is 0.05.
#' @references Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.
#' @references Wood, S. (2017). Additive generalized models: An R introduction. Chapman and Hall.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @import graphics
#' @import ggplot2
#' @examples
#' set.seed(1)
#' n <- 300
#' error <- rglg(n,0,0.5,1)
#' t <- as.matrix((2*1:n - 1)/(2*n))
#' colnames(t) <- "t"
#' f_t <- cos(4*pi*t)
#' y <- 0.8 + f_t + error
#' colnames(y) <- "y"
#' data <- as.data.frame(cbind(y,1,t))
#' fit1 <- sglg(y ~ 1,npc=t,data=data,basis = "deBoor",alpha0=0.0001)
#' summary(fit1)
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
        #geom_line(aes(npc,f_est_up),colour = "orange",size=1.2) +

        geom_line(aes(npc,f_est),size=1.2) +

        #geom_line(aes(npc,f_est_low),colour = "orange",size=1.2) +

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
    #geom_point(colour="blue",alpha=0.4)+
    geom_line(aes(npc,f_est_up),colour = "orange",size=1.2) +
    geom_line(aes(npc,f_est),size=1.2) +
    geom_line(aes(npc,f_est_low),colour = "orange",size=1.2) +
    xlab(colnames(npc))+
    ggtitle(add_comp)
    return(plot)
}


















