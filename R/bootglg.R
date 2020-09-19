#' Bootstrap inference for a generalized log-gamma distribution
#'
#' \code{bootglg} is used to generate bootstrap inference, such as, estimations, standard errors, intervals for the parameters of a generalized log-gamma distribution.
#' @param y vector, a numeric random sample.
#' @param R integer, represents the number of replications. Default value is 1000.
#' @param alpha numeric, represents a confidence level for the bootstrap intervals. Default value is 0.05.
#' @param method character, indicates the type of bootstrap 'Nonparametric' or 'Parametric'. Default value is 'Parametric'.

#' @return ml_estimates is a vector of maximum likelihood estimates asociated with the location,scale and shape parameters.
#' @return boot_mean_estimates is a vector of mean of the bootstrap estimates asociated with the location,scale and shape parameters.
#' @return boot_sd_estimates is a vector of bootstrap standard errors of the estimates asociated with the location,scale and shape parameters.

#' @references Cardozo C, Paula G and Vanegas L. sglg: An R package for fitting semi-parametric generalized log-gamma regression models. In preparation.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' set.seed(1)
#' y <- rglg(100,location=1,scale=0.5,shape=1)
#' \dontrun{
#' bootglg(y,R=300,method='Parametric')
#' }
#' @export bootglg
bootglg = function(y, R, alpha,method){
  if (missingArg(y))
    y <- rglg(100,location=1,scale=0.5,shape=1)
  if (missingArg(R))
    R <- 1000
  if (missingArg(alpha))
    alpha <- 0.05
  if (missingArg(method))
    method <- 'Parametric'

  n <- length(y)
  fit0 <- try(glg(y~1,data = as.data.frame(y)),silent=TRUE)
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
          fit <- try(glg(boot_y~1,data = as.data.frame(boot_y)),silent=TRUE)
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
      fit <- try(glg(boot_y~1,data = as.data.frame(boot_y)),silent=TRUE)
      if(is.list(fit)){
        boot_coef_est[i,] <- c(fit$mu,fit$sigma,fit$lambda)
        print(i)
        i <- i+1
      }
    }
  }


  boot_coef_est <- as.data.frame(boot_coef_est)
  ext <- max(abs(boot_coef_est$mu)) + 0.5
  plot1 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'mu'])) + ggtitle("Histogram Location Paramter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  ext <- max(abs(boot_coef_est$sigma)) + 0.5
  plot2 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'sigma'])) + ggtitle("Histogram Scale Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  ext <- max(abs(boot_coef_est$lambda)) + 0.5
  plot3 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'lambda'])) + ggtitle("Histogram Shape Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(-ext,ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  grid.arrange(plot1, plot2, plot3,ncol=3)
  output <- list(method=method,ml_estimates=round(ml_estimates,4),boot_mean_estimates=round(apply(boot_coef_est,2,mean),4),boot_estimates_sd = round(apply(boot_coef_est,2,sd),4))
  return(output)
}
