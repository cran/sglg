#'Plotting a natural cubic splines or P-splines.
#'
#'\code{plotnpc} displays a graph of a fitted nonparametric effect, either natural cubic spline or P-spline, from an object of class sglg.
#'
#' @param fit an object of the class sglg. This object is returned from the call to sglg() or ssurvglg().
#' @param conf_lev is the confidence level of the asymptotic confidence band. Default value is 0.05.
#' @param cband is a boolean value. It indicates if the plot will contain confidence band or not. Default value is FALSE.
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
#' error <- rglg(rows,0,1,t_lambda)
#' y <- X %*%t_beta + f_t + t_sigma*error
#' colnames(y) <- "y"
#' data <- data.frame(y,X,t)
#' fit1 <- sglg(y ~ x1 + x2 - 1, npc=t, data=data, basis = "deBoor", Knot= 5, alpha0=0.1)
#' # The adjusted (black) non-linear component
#' plotnpc(fit1)
#' plotnpc(fit1,conf_lev=0.02,cband=TRUE)
#' fit2 <- sglg(y ~ x1 + x2 - 1, npc=t, data=data, basis = "Gu", Knot= 5, alpha0=0.001)
#' # The adjusted (black) non-linear component
#' plotnpc(fit2,conf_lev=0.02,cband=TRUE)
#' @export plotnpc

plotnpc <- function(fit,conf_lev,cband=FALSE) {
    if (fit$semi == FALSE) {
        stop("Sorry, for this kind of model it is not available this option.")
    }
    if(missingArg(conf_lev))
       conf_lev <- 0.05

    y <- fit$y
    X <- fit$X
    k <- fit$Knot
    N <- fit$N
    npc <- fit$npc
    add_comp <- paste("f(",colnames(npc),sep="")
    add_comp <- paste(add_comp,")",sep="")
    add_comp <- paste("Estimated Additive Component", add_comp,sep=" ")
    p <- fit$p
    mu <- fit$mu
    #conf_lev_ <- conf_lev/k

    f_est <- N[,-(1:p)]%*%mu[-(1:p)]
    df <- as.data.frame(cbind(y,N,f_est))
    if (fit$basis == "deBoor"){
        if(cband==TRUE){
          var_gammas <- fit$scovar[(p+1):(p+k),(p+1):(p+k)]
          var_f_est <- N[,-(1:p)]%*%var_gammas%*%t(N[,-(1:p)])
          st_error_f_est <- sqrt(diag(var_f_est))
          f_est_low <- f_est + qnorm(0.5*conf_lev)*st_error_f_est
          f_est_up <- f_est + qnorm(1 - 0.5*conf_lev)*st_error_f_est
          plot <- ggplot(data=df,aes(npc,y - X%*%mu[1:p]))+
          geom_point(colour="blue",alpha=0.4)+
          geom_line(colour="orange", aes(npc,f_est_up),size=0.9) +
          geom_line(aes(npc,f_est),size=0.75) +
          geom_line(colour="orange", aes(npc,f_est_low),size=0.9) +
          xlab(colnames(npc))+
          ylab('')+
          ggtitle(add_comp)
          return(plot)
        }
        plot <- ggplot(data=df,aes(npc,y - X%*%mu[1:p]))+
          geom_point(colour="blue",alpha=0.4)+
          geom_line(aes(npc,f_est),size=0.75) +
          xlab(colnames(npc))+
          ylab('')+
          ggtitle(add_comp)
        return(plot)
    }
    if(cband==TRUE){
      var_gammas <- fit$scovar[(p+1):(p+k),(p+1):(p+k)]
      var_f_est <- N[,-(1:p)]%*%var_gammas%*%t(N[,-(1:p)])
      st_error_f_est <- sqrt(diag(var_f_est))
      f_est_low <- f_est + qnorm(0.5*conf_lev)*st_error_f_est
      f_est_up <- f_est + qnorm(1 - 0.5*conf_lev)*st_error_f_est
      df <- as.data.frame(cbind(y,N,f_est))
      plot <- ggplot(data=df,aes(npc,y- X%*%mu[1:p]))+
      geom_point(colour="blue",alpha=0.4)+
      geom_line(aes(npc,f_est_up),colour = "orange",size=0.9) +
      geom_line(aes(npc,f_est),size=0.75) +
      geom_line(aes(npc,f_est_low),colour = "orange",size=0.9) +
      xlab(colnames(npc))+
      ylab('')+
      ggtitle(add_comp)
      return(plot)
    }
    plot <- ggplot(data=df,aes(npc,y- X%*%mu[1:p]))+
    geom_point(colour="blue",alpha=0.4)+
    geom_line(aes(npc,f_est),size=0.75) +
    xlab(colnames(npc))+
    ylab('')+
    ggtitle(add_comp)
    return(plot)
}


















