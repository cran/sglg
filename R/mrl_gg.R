#' Mean Residual Lifetime Function for a Generalized Gamma Distribution
#'
#'\code{mrl_gg} is used to obtain the value of the mean residual lifetime function of a generalized gamma distribution at a positive value.
#' @param x numeric, represent a vector of positive values. Default value is 1.
#' @param mu numeric, represents the location parameter of a generalized gamma distribution. Default value is 0.
#' @param sigma numeric, represents the scale parameter of a generalized gamma distribution. Default value is 1.
#' @param lambda numeric, represents the shape parameter of a generalized gamma distribution. Default value is 1.
#' @return A numeric value of the mean residual lifetime of a generalized gamma distribution.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph.D. thesis. Sao Paulo University.
#' @references Jerald F. Lawless (2003). Statistical Models and Methods for Lifetime Data. Second Edition. John-Wiley & Sons
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' mrl_gg(x=0,mu=0,sigma=2,lambda=1) # Extreme value type I distribution, maximum case.
#' @export mrl_gg
mrl_gg = function(x = 1, mu=0, sigma=1, lambda=1) {
  if (sum(x < 0) > 0)
    stop('x must be a positive value or a vector of positive values!')
  surv <- function(x){
    return(1 - pglg(log(x), mu, sigma, lambda))
  }
  mean_residual <- function(x){
    output <- integrate(surv,x,Inf)$value/surv(x)
    return(output)
  }
  mean_residual_v <- Vectorize(mean_residual, vectorize.args = 'x')
  return(mean_residual_value = mean_residual_v(x))
}
