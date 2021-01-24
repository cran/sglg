#' Bootstrap inference for a generalized log-gamma distribution
#'
#' \code{bootglg} is used to generate bootstrap inference, such as, estimated standard errors and approximate confidence intervals for the parameters of a generalized log-gamma distribution.
#' @param y vector, a numeric random sample.
#' @param R integer, represents the number of replications. Default value is 1000.
#' @param alpha numeric, represents a confidence level for the bootstrap intervals. Default value is 0.05.
#' @param method character, indicates the type of bootstrap 'Nonparametric' or 'Parametric'. Default value is 'Parametric'.
#' @param type character, indicates the type of bootstrap confidence interval for the estimated parameters. The options are:
#' 'normal', 't_student' or 'bootstrap_t'. These intervals used the bootstrap estimated standard error of the ML estimates of the parameters.
#' Other kind of bootstrap intervals are the percentile-type intervals. We offer the option 'BCa'. It is a bias-corrected and accelerated percentile interval.
#'  The default value for the 'type' argument is 'normal'.
#' @return \code{ml_estimates} is a vector of maximum likelihood estimates associated with the location, scale, and shape parameters.
#' @return \code{boot_mean_estimates} is a vector of mean of the bootstrap estimates associated with the location, scale, and shape parameters.
#' @return \code{boot_sd_estimates} is a vector of bootstrap standard errors of the estimates associated with the location, scale, and shape parameters.
#' @return \code{type} indicates the type of confidence intervals.
#' @return \code{intervals} array of the confidence intervals of the location,scale and shape.
#' @references Cardozo C, Paula G and Vanegas L. sglg: An R package to fit semi-parametric generalized log-gamma regression models. In preparation.
#' @references Efron B and Tibshirani R (1993). An introduction to the Bootstrap. Chapman-Hall.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' set.seed(1)
#' y <- rglg(200,location=1,scale=0.5,shape=1)
#' \dontrun{
#' bootglg(y,R=300,method='Parametric',type='normal')
#' bootglg(y,R=300,method='Nonparametric',type='t_student')
#' bootglg(y,R=300,method='Parametric',type='bootstrap_t')
#' bootglg(y,R=300,method='Nonparametric',type='BCa')
#' }
#' @import moments
#' @export bootglg
bootglg = function(y, R, alpha, method, type){
  if (missingArg(y))
    y <- rglg(100,location=1,scale=0.5,shape=1)
  if (missingArg(R))
    R <- 1000
  if (missingArg(alpha))
    alpha <- 0.05
  if (missingArg(method))
    method <- 'Parametric'
  if (missingArg(type))
    method <- 'normal'

  n <- length(y)
  fit0 <- try(glg(y~1,data = as.data.frame(y),format='simple'),silent=TRUE)
  print(fit0)
  ml_estimates <- vector()

  if(is.list(fit0)){
    ml_estimates <- c(fit0$mu,fit0$sigma,fit0$lambda)
  }
  boot_coef_est <- matrix(0,nrow=R,ncol=3)
  colnames(boot_coef_est) <- c("mu","sigma","lambda")

  if(method == 'Nonparametric'){
     i <- 1
     while(i <= R){
          boot_y <- sample(y,n,replace=TRUE)
          fit <- try(suppressWarnings(glg(boot_y~1,data = as.data.frame(boot_y),format='simple')),silent=TRUE)
          if(is.list(fit)){
          boot_coef_est[i,] <- c(fit$mu,fit$sigma,fit$lambda)
          print(i)
          i <- i+1
         }
    }
  }
  if(method == 'Parametric'){
    i <- 1
    while(i <= R){
      boot_y <- rglg(n,location=fit0$mu,scale=fit0$sigma,shape=fit0$lambda)
      fit <- try(suppressWarnings(glg(boot_y~1,data = as.data.frame(boot_y),format='simple')),silent=TRUE)
      if(is.list(fit)){
        boot_coef_est[i,] <- c(fit$mu,fit$sigma,fit$lambda)
        print(i)
        i <- i+1
      }
    }
  }

  boot_coef_est <- as.data.frame(boot_coef_est)
  ext <- max(abs(boot_coef_est$mu)) + 0.5
  plot1 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'mu'])) + ggtitle("Histogram Location Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  ext <- max(abs(boot_coef_est$sigma)) + 0.5
  plot2 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'sigma'])) + ggtitle("Histogram Scale Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  ext <- max(abs(boot_coef_est$lambda)) + 0.5
  plot3 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'lambda'])) + ggtitle("Histogram Shape Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  grid.arrange(plot1, plot2, plot3,ncol=3)

  boot_mean_estimates <- apply(boot_coef_est,2,mean)
  boot_estimates_sd <- apply(boot_coef_est,2,sd)

  bts_t <- function(alpha){
           lr <- matrix(0,length(ml_estimates),2)
           for(i in 1:length(ml_estimates)){
           z_l <- stats::quantile((boot_coef_est[,i] - ml_estimates[i])/boot_estimates_sd[i],alpha/2)
           z_r <- stats::quantile((boot_coef_est[,i] - ml_estimates[i])/boot_estimates_sd[i],1-(alpha/2))
           lr[i,1] <- ml_estimates[i] + z_l*boot_estimates_sd[i]
           lr[i,2] <- ml_estimates[i] + z_r*boot_estimates_sd[i]
           }
           return(lr)
  }

  BCa <- function(alpha){
    lr <- matrix(0,length(ml_estimates),2)
    jn_coef_est <- matrix(0,nrow=n,ncol=3)
    for(j in 1:n){
        y_j <- y[-j]
        fit <- try(suppressWarnings(glg(y_j~1,data = as.data.frame(y_j),format='simple')),silent=TRUE)
        if(is.list(fit)){
          jn_coef_est[j,] <- c(fit$mu,fit$sigma,fit$lambda)
        }
    }

    a_est <- (1/6)*apply(jn_coef_est,2,moments::skewness)
    for(i in 1:length(ml_estimates)){
      z_0_est <- qnorm(mean( boot_coef_est[,i] < ml_estimates[i]))
      alpha_1 <- pnorm( z_0_est +  ( ( z_0_est + qnorm(alpha/2) )/(1 - a_est[i]*(z_0_est + qnorm(alpha/2))) ) )
      alpha_2 <- pnorm( z_0_est +  ( ( z_0_est + qnorm(1 - alpha/2) )/(1 - a_est[i]*(z_0_est + qnorm(1 - alpha/2))) ) )
      lr[i,1] <- stats::quantile(boot_coef_est[,i],alpha_1)
      lr[i,2] <- stats::quantile(boot_coef_est[,i],alpha_2)
    }
    return(lr)
  }

  type_intervals <- function(type){
    out_type_interval <- vector()
    switch(type,
           normal = matrix(c(ml_estimates + qnorm(alpha/2)*boot_estimates_sd,boot_mean_estimates + qnorm(1 - (alpha/2))*boot_estimates_sd),3,2,byrow=FALSE),
           t_student = matrix(c(ml_estimates + qt(alpha/2,n-2)*boot_estimates_sd,boot_mean_estimates + qt(1 - (alpha/2),n-2)*boot_estimates_sd),3,2,byrow=FALSE),
           bootstrap_t = bts_t(alpha),
           BCa = BCa(alpha))
  }
  type_interval <- paste(method," approximate ",type," interval")
  output <- list(method=method,ml_estimates=ml_estimates,boot_estimates_sd=boot_estimates_sd,type_interval = type_interval,intervals = type_intervals(type))
  return(output)
}
