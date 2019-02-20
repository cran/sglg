#' Random Sampling of Order Statistics from a Generalized Log-gamma Distribution
#'
#'\code{order_glg} is used to obtain a random sample of order statistics from a Generalized Log-gamma Distribution.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter. Default value is 0.
#' @param sigma numeric, represents the scale parameter. Default value is 1.
#' @param lambda numeric, represents the shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statiscal Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples{
#' # A random sample of size 10 of order statistics from a Extreme Value Distribution.
#' order_glg(10,0,1,1,1,50)}
#' @importFrom stats rbeta
#' @export order_glg

order_glg <- function(size,mu,sigma,lambda,k,n){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qloggamma(initial,mu,sigma,lambda)
  pdf     <- factorial(size)*cumprod(dloggamma(sample,mu,sigma,lambda))[size]
  return(list(sample=sample,pdf=pdf))
}
