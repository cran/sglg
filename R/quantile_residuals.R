#' Quantile Residuals for a Generalized Log-gamma Regression Model
#'
#' \code{quantile_residuals} is used to generate quantile residuals for a generalized log-gamma regression model.
#' @param fit is an object sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg().
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' # Example 1
#'n <- 300
#'error <- rglg(n,0,1,1)
#'y <- 0.5 + error
#'fit <- glg(y~1,data=as.data.frame(y))
#'quantile_residuals(fit)
#'# Example 2
#'n <- 300
#'error <- rglg(n,0,1,1)
#'x <- runif(n,-3,3)
#'y <- 0.5 +  2*x + error
#'fit <- glg(y~x,data=as.data.frame(y,x))
#'quantile_residuals(fit)
#'@import moments
#'@import ggplot2
#'@import gridExtra
#'@export quantile_residuals
quantile_residuals <- function(fit){
r_quantile <- as.numeric(qnorm( pglg(fit$y,location=fit$mu_est,scale=fit$sigma,shape=fit$lambda)))
cond_1 <- r_quantile == Inf
if(length(cond_1) > 0)
  r_quantile[r_quantile == Inf] <- max(r_quantile[r_quantile != Inf]) + 1e-03

cond_2 <- r_quantile == -Inf
if(length(cond_2) > 0)
  r_quantile[r_quantile == -Inf] <- min(r_quantile[r_quantile != -Inf]) - 1e-03

ext <- max(abs(r_quantile)) + 0.5
plot1 <- ggplot(data=as.data.frame(r_quantile),aes(r_quantile)) +  ggtitle("Density Quantile Residuals") + geom_density(colour="orange",fill="orange",alpha=0.3) + xlab("Sample Quantiles") + ylab("Density") + xlim(c(-ext,ext)) + geom_hline(yintercept=0)
plot2 <- ggplot(data=as.data.frame(r_quantile),aes(sample=r_quantile)) + ggtitle("Normal Q-Q Plot") + stat_qq(colour="blue",alpha=0.5) + stat_qq_line() + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
grid.arrange(plot1, plot2, ncol=2)
output <- list(mean=round(mean(r_quantile),2),sd= round(sd(r_quantile),2),skew =  round(skewness(r_quantile),2), kurt = round(kurtosis(r_quantile),2))
return(output)
}
