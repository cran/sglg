#' Quantile Residuals for a Generalized Log-gamma Regression Model
#'
#' \code{quantile_sglg} is used to generate quantile residuals for a generalized log-gamma regression model.
#' @param fit is an object sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg().
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' # Example 1
#'n <- 400
#'set.seed(4)
#'error <- rglg(n,0,0.5,1)
#'y <- as.data.frame(0.5 + error)
#'names(y) <- "y"
#'fit_0 <- glg(y~1,data=y)
#'fit_0$mu
#'fit_0$sigma
#'fit_0$lambda
#'quantile_sglg(fit_0)
#'# Example 2
#'n <- 500
#'set.seed(6)
#'error <- rglg(n,0,0.5,1)
#'x1 <- runif(n,-2,2)
#'beta <- c(0.5,2)
#'y <- cbind(1,x1)%*%beta + error
#'data <- data.frame(y=y,x1=x1)
#'fit_1 <- glg(y~x1,data=data)
#'fit_1$mu
#'fit_1$sigma
#'fit_1$lambda
#'quantile_sglg(fit_1)
#'@import moments
#'@import ggplot2
#'@import gridExtra
#'@export quantile_sglg
quantile_sglg <- function(fit){
r_quantile <- quantile(as.numeric(qnorm( pglg(fit$y,location=fit$mu_est,scale=fit$sigma,shape=fit$lambda))),seq(0.02,0.98,by=0.01))
cond_1 <- r_quantile == Inf
if(sum(cond_1) > 0)
   r_quantile[r_quantile == Inf] <- 1.005*max(r_quantile[r_quantile != Inf])
cond_2 <- r_quantile == -Inf
if(sum(cond_2) > 0)
  r_quantile[r_quantile == -Inf] <- 0.995*min(r_quantile[r_quantile != -Inf])

left <- 0.995*min(r_quantile)
rigth <- 1.005*max(r_quantile)
plot1 <- ggplot(data=as.data.frame(r_quantile),aes(r_quantile)) +  ggtitle("Density Quantile Residuals") + geom_density(colour="orange",fill="orange",alpha=0.3) + xlab("Sample Quantiles") + ylab("Density") + xlim(c(left,rigth)) + geom_hline(yintercept=0)
plot2 <- ggplot(data=as.data.frame(r_quantile),aes(sample=r_quantile)) + ggtitle("Normal Q-Q Plot") + stat_qq(colour="blue",alpha=0.5) + stat_qq_line(line.p = c(0.05, 0.95)) + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
grid.arrange(plot1, plot2, ncol=2)
output <- list(mean=round(mean(r_quantile),2),sd= round(sd(r_quantile),2),skew =  round(skewness(r_quantile),2), kurt = round(kurtosis(r_quantile),2))
return(output)
}
