#' Density Probability Distribution of a K-th Order Statistic from a Generalized Log-gamma Distribution
#'
#'\code{dk_order_glg} is used to obtain the density probability distribution of the k-th order statistic from a generalized log-gamma distribution.
#' @param x numeric, represents a real value.
#' @param mu numeric, represents the location parameter. Default value is 0.
#' @param sigma numeric, represents the scale parameter. Default value is 1.
#' @param lambda numeric, represents the shape parameter. Default value is 1.
#' @param k numeric, represents the K-th smallest value from a sample.
#' @param n numeric, represents the size of the sample of the generalized log-gamma distribution.
#' @return A list of values of the density probability function of the k-th order statistic from a generalized log-gamma distribution.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' # The density probability distribution of 10-th order statistics at 0
#' # from a random sample of extreme value distribution with n=20.
#' dk_order_glg(0,0,1,1,k=10,n=20)
#' @importFrom stats rbeta
#' @export dk_order_glg

dk_order_glg <- function(x,mu=0,sigma=1,lambda=1,k,n){
  p <- pglg(x,mu,sigma,lambda)
  dpdf    <-  k*choose(n,k)*dglg(x,mu,sigma,lambda)*(p^(k-1))*((1 - p)^(n-k))
  return(dpdf)
}
