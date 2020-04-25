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
#' @examples
#' # A random sample of size 10 of order statistics from a Extreme Value Distribution.
#' order_glg(10,0,1,1,1,50)
#' \dontrun{ # A small comparison between two random sampling methods of order statistics
#' # Method 1
#' m <- 10
#' output <- rep(0,m)
#' order_sample <- function(m,n,k){
#' for(i in 1:m){
#' sample <- rglg(n)
#' order_sample <- sort(sample)
#' output[i] <- order_sample[k]
#' }
#' return(output)
#' }
#' N <- 10000
#' n <- 200
#' k <- 100
#' system.time(order_sample(N,n,k))
#' sample_1 <- order_sample(N,n,k)
#' hist(sample_1)
#' summary(sample_1)
#' # Method 2
#' system.time(order_glg(N,0,1,1,k,n))
#' sample_2 <- order_glg(N,0,1,1,k,n)$sample
#' hist(sample_2)
#' summary(sample_2)
#' }
#' @importFrom stats rbeta
#' @export order_glg

order_glg <- function(size,mu,sigma,lambda,k,n){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qglg(initial,mu,sigma,lambda)
  #pdf    <- factorial(size)*cumprod(dloggamma(sample,mu,sigma,lambda))[size]
  pdf    <- factorial(size)*cumprod(dglg(sample,mu,sigma,lambda))[size]
  return(list(sample=sample,pdf=pdf))
}
