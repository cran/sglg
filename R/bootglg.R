#' Bootstrap inference for a generalized log-gamma regression
#'
#' \code{bootglg} is used to generate parametric bootstrap inference, such as, estimated standard errors and approximate confidence intervals for a generalized log-gamma regression.
#' @param formula a symbolic description of the systematic component of the model to be fitted.
#' @param data data.frame, contains the variables in the formula object.
#' @param B integer, represents the number of bootstrap replications. Default value is 500.
#' @param alpha numeric, represents a confidence level for the bootstrap intervals. Default value is 0.05.
#' @param type character, indicates the type of bootstrap confidence interval for the estimated parameters. The options are:
#' 'normal', 't_student' or 'bootstrap_t'. These intervals used the bootstrap estimated standard error of the ML estimates of the parameters.
#' Other kind of bootstrap intervals are the percentile-type intervals. We offer the option 'BCa'. It is a bias-corrected and accelerated percentile interval.
#'  The default value for the 'type' argument is 'normal'.
#' @param plt_den boolean value, to request a density-type plot of the bootstrap estimates. Default value is FALSE.
#' @return \code{ml_estimates} is a vector of maximum likelihood estimates associated with the coefficients of linear structure, scale, and shape parameters.
#' @return \code{boot_mean_estimates} is a vector of mean of the bootstrap estimates associated with the coefficients of linear structure, scale, and shape parameters.
#' @return \code{boot_bias_estimates} is a vector of bootstrap estimate of bias associated with the coefficients of linear structure, scale, and shape estimators.
#' @return \code{boot_sd_estimates} is a vector of bootstrap standard errors of the estimates associated with the coefficients of linear structure, scale, and shape estimators.
#' @return \code{type} indicates the type of confidence intervals.
#' @return \code{intervals} array of the confidence intervals of the coefficients of linear structure, scale, and shape.
#' @references Cardozo C. A., Paula G. and Vanegas L. sglg: An R package to fit semi-parametric generalized log-gamma regression models. In preparation.
#' @references Efron B and Tibshirani R (1993). An introduction to the Bootstrap. Chapman & Hall, Inc.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' ##################################################################################################
#' set.seed(1)
#' n <- 300
#' x1 <- runif(n, 0, 1)
#' t_beta  <- 1.2
#' t_sigma <- 0.5
#' t_lambda <- 0.7
#' error <- rglg(n, 0, t_sigma, t_lambda)
#' y1 <- t_beta*x1 + error
#' data <- data.frame(y1, x1)
#' # The following examples are based on 50 bootstrap replications.
#' # A 90% bootstrap confidence interval with the method 'normal'.
#' bootglg(y1 ~ x1 - 1, data = data, type='normal', B = 50, alpha = 0.1)
#' # A 95% bootstrap confidence interval with the method 't_student'.
#' bootglg(y1 ~ x1 - 1, data = data, type='t_student', B = 50)
#' # A 95% bootstrap confidence interval with the method 'bootstrap_t'.
#' bootglg(y1 ~ x1 - 1, data = data, type='bootstrap_t', B = 50)
#' # A 98% bootstrap confidence interval with the method 'BCa'.
#' # bootglg(y1 ~ x1 - 1, data = data, type='BCa', B = 50, alpha = 0.02)
#' #################################################################################################
#' @import moments
#' @import progress
#' @export bootglg
bootglg = function(formula, data, B = 500, alpha = 0.05, type = 'normal', plt_den = FALSE){

  fit0 <- try(glg(formula, data = data), silent=TRUE)
  y_base <- fit0$y
  X <- as.matrix(fit0$X)
  p <- dim(X)[2]
  n <- dim(data)[1]
  loc <- fit0$mu
  scale <- fit0$sigma
  shape <- fit0$lambda
  half_alpha <- alpha/2

  if(is.list(fit0)){
    ml_estimates <- c(loc, scale, shape)
  }

  k <- length(ml_estimates)
  boot_coef_est <- matrix(0, nrow = B, ncol = k)
  par_names <- c(paste("beta_", 1:p, sep = ''), "sigma", "lambda")
  colnames(boot_coef_est) <- par_names
  progress_index <- round(quantile(seq(1:B), probs = seq(0, 1, 0.2)))[-1]
  progress_names <- names(progress_index)

  i <- 1
  print(progress_names[1])
  j <- 2
  while(i <= B){
    error <- rglg(n, 0, 1, shape)
    y <- X %*%loc + scale * error
    colnames(y) = as.character(formula)[2]
    fit <- try(suppressWarnings(glg(formula, data = as.data.frame(y,X), format='simple')), silent=TRUE)
    if(is.list(fit)){
      if(i==progress_index[j]){
          print(progress_names[j])
          j <- j+1
        }
      boot_coef_est[i,] <- c(fit$mu,fit$sigma,fit$lambda)
      i <- i+1
      }
    }

  boot_coef_est <- as.data.frame(boot_coef_est)

  #l_ext <- 0.98*min(boot_coef_est$mu)
  #u_ext <- 1.02*max(boot_coef_est$mu)
  #plot1 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'mu'])) + ggtitle("Histogram Location Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(l_ext,u_ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  #l_ext <- 0.98*min(boot_coef_est$sigma)
  #u_ext <- 1.02*max(boot_coef_est$sigma)
  #plot2 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'sigma'])) + ggtitle("Histogram Scale Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(l_ext,u_ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  #l_ext <- 0.98*min(boot_coef_est$lambda)
  #u_ext <- 1.02*max(boot_coef_est$lambda)
  #plot3 <- ggplot(data=boot_coef_est,aes(boot_coef_est[,'lambda'])) + ggtitle("Histogram Shape Parameter") + geom_density(colour="orange",fill="orange",alpha=0.25) + xlim(c(l_ext,u_ext)) + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
  #grid.arrange(plot1, plot2, plot3,ncol=3)

  #############################################################################################
  #if(plt_den == TRUE){
  #  plt <- ggplot()
  #  for(j in 1:k){
  #    l_ext <- 0.98*min(boot_coef_est[,j])
  #    u_ext <- 1.02*max(boot_coef_est[,j])
  #    x <- data.frame(b_est=boot_coef_est[,par_names[j]])
  #    plt <- ggplot(aes(b_est), data=x) + geom_density(colour="orange", fill="orange", alpha=0.25) + ggtitle(paste("Parameter ", par_names[j]))  + xlab("Bootstrap Estimates") + ylab("Density") + geom_hline(yintercept=0)
      #grid.arrange(plt, ncol= p+2)
  #  }
  #}
  #############################################################################################

  boot_mean_estimates <- apply(boot_coef_est, 2, mean)
  boot_estimates_sd <- apply(boot_coef_est, 2, sd)

  bts_t <- function(alpha){
           lr <- matrix(0,k,2)
           for(i in 1:k){
             values <- (boot_coef_est[,i] - ml_estimates[i])/boot_estimates_sd[i]
           z_l <- stats::quantile(values, half_alpha)
           z_r <- stats::quantile(values, 1 - half_alpha)
            lr[i,1] <- ml_estimates[i] + z_l*boot_estimates_sd[i]
            lr[i,2] <- ml_estimates[i] + z_r*boot_estimates_sd[i]
           }
           return(lr)
  }

  BCa <- function(alpha){
    print("Second level calculation of the BCa method")
    lr <- matrix(0, k, 2)
    jn_coef_est <- matrix(0, nrow = n, ncol = k)
    progress_indicator <- progress_bar$new(format = "  bootglg_jackknife [:bar] :percent in :elapsed",
                                           total = n,
                                           clear = FALSE,
                                           width= 60)
    for(j in 1:n){
        y_j <- as.matrix(y_base[-j])
        colnames(y_j) = as.character(formula)[2]
        X_j <- as.matrix(X[-j,])
        colnames(X_j) <- colnames(X)
        j_data <- cbind(y_j,  X_j)
        fit <- try(suppressWarnings(glg(formula, data = as.data.frame(j_data), format='simple')), silent=TRUE)
      if(is.list(fit)){
          jn_coef_est[j,] <- c(fit$mu,fit$sigma,fit$lambda)
        }
        progress_indicator$tick()
        Sys.sleep(1 / 100)
    }

    a_est <- (1/6)*apply(jn_coef_est, 2, skewness)
    for(i in 1:k){
      z_0_est <- qnorm(mean( boot_coef_est[,i] < ml_estimates[i]))
      const <- (1 - a_est[i]*(z_0_est + qnorm(half_alpha)))
      alpha_1 <- pnorm( z_0_est +  ( z_0_est + qnorm(half_alpha) )/ const)
      alpha_2 <- pnorm( z_0_est +  ( z_0_est + qnorm(1 - half_alpha) )/const)
      if( alpha_1 <  alpha_2){
        lr[i,1] <- stats::quantile(boot_coef_est[,i], alpha_1)
        lr[i,2] <- stats::quantile(boot_coef_est[,i], alpha_2)
      }
      else{
        lr[i,1] <- stats::quantile(boot_coef_est[,i], alpha_2)
        lr[i,2] <- stats::quantile(boot_coef_est[,i], alpha_1)
      }
    }
    return(lr)
  }

  type_intervals <- function(type){
    out_type_interval <- vector()
    switch(type,
           normal =    matrix(c(ml_estimates + qnorm(half_alpha)*boot_estimates_sd, ml_estimates + qnorm(1 - half_alpha)*boot_estimates_sd), k, 2, byrow=FALSE),
           t_student = matrix(c(ml_estimates + qt(half_alpha, n - 1)*boot_estimates_sd, ml_estimates + qt(1 - half_alpha, n - 1)*boot_estimates_sd), k, 2, byrow=FALSE),
           bootstrap_t = bts_t(alpha),
           BCa = BCa(alpha))
  }
  boot_estimates_bias <- boot_mean_estimates - ml_estimates
  type_interval <- paste("Approximate", paste((1-alpha)*100, '%', sep=''), type, "confidence interval")
  output <- list(method = 'Parametric',
                 ml_estimates = ml_estimates,
                 boot_replicates = B,
                 boot_estimates_bias = boot_estimates_bias,
                 boot_estimates_sd = boot_estimates_sd,
                 type_interval = type_interval,
                 intervals = type_intervals(type))
  return(output)
}
