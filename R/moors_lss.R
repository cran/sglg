#' Measures of location, scale, and shape based on quantile measures for a generalized log-gamma distribution
#'
#'\code{moors_lss} is used to obtain the median, the half interquartile range and the quantile coefficient of skewness and kurtosis for a generalized log-gamma distribution.
#' @param mu numeric, represents the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param sigma numeric, represents the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param lambda numeric, represents the shape parameter of a generalized log-gamma distribution. Default value is 1.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @references J. J. A. Moors (1988), A quantile alternative for kurtosis. The Statistician.
#'
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' moors_lss(mu = 0,sigma = 1,lambda = -1)    # Extreme value type I distribution, maximum case.
#' moors_lss(mu = 0,sigma = 1,lambda = 1)     # Extreme value type I distribution, minimum case.
#' moors_lss(mu = 0,sigma = 1,lambda = 0.05) # Standard normal distribution.
#' @export moors_lss
moors_lss = function(mu = 0, sigma = 1, lambda = 1) {
  lambda_2 <- lambda^2
  m_lambda_2 <- 1/lambda_2
  median <- qglg(0.5,mu,sigma,lambda)
  half_range <- (qglg(0.75,mu,sigma,lambda) - qglg(0.25,mu,sigma,lambda))/2
  skewness <- (qglg(0.75,mu,sigma,lambda) - 2*median + qglg(0.25,mu,sigma,lambda))/half_range
  kurtosis <- ((qglg(7/8,mu,sigma,lambda) - qglg(5/8,mu,sigma,lambda)) + (qglg(3/8,mu,sigma,lambda)-qglg(1/8,mu,sigma,lambda)))/(qglg(3/4,mu,sigma,lambda)-qglg(1/4,mu,sigma,lambda))
  return(list(median = median, half_range = half_range, moors_skewness = skewness, moors_kurtosis = kurtosis))
}
