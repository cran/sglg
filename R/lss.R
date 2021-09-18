#' Measures of location, scale and shape measures for a generalized log-gamma distribution
#'
#'\code{lss} is used to obtain the mean, median, mode, variance, coefficient of variation, skewness and kurtosis for a generalized log-gamma distribution.
#' @param mu numeric, represents the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param sigma numeric, represents the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param lambda numeric, represents the shape parameter of a generalized log-gamma distribution. Default value is 1.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @references National Institute of Standards and Technology, NIST.  Engineering Statistics Handbook.
#' https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
#'
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' lss(0,1,-1)    # Extreme value type I distribution, maximum case.
#' lss(0,1,1)     # Extreme value type I distribution, minimum case.
#' lss(0,1,0.01) # Standard normal distribution.
#' @export lss
lss = function(mu, sigma, lambda) {
    if (missingArg(mu))
        mu <- 0
    if (missingArg(sigma))
        sigma <- 1
    if (missingArg(lambda))
        lambda <- 1

    lambda_2 <- lambda^2
    m_lambda_2 <- 1/lambda_2
    mean <- mu + sigma * ((digamma(m_lambda_2) - log(m_lambda_2))/lambda)
    median <- qglg(0.5,mu,sigma,lambda)
    variance <- (sigma^2) * (trigamma(m_lambda_2)/lambda_2)
    if(abs(lambda) < 0.1)
      cv <- "It is not defined because the mean of this distribution is too close to zero!"
    else
      cv <- sqrt(variance)/mean
    skewness <- sign(lambda) * psigamma(m_lambda_2, deriv = 2)/(trigamma(m_lambda_2)^(3/2))
    kurtosis <- (psigamma(m_lambda_2, deriv = 3)/((trigamma(m_lambda_2))^2)) + 3
    return(list(mean = mean, median= median, mode = mu, variance = variance, coef_var = cv, skewness = skewness, kurtosis = kurtosis))
}
