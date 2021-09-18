#' Random Sampling of K-th Order Statistics from a Generalized Log-gamma Distribution
#'
#'\code{order_glg} is used to obtain a random sample of the K-th order statistics from a generalized log-gamma distribution.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter. Default value is 0.
#' @param sigma numeric, represents the scale parameter. Default value is 1.
#' @param lambda numeric, represents the shape parameter. Default value is 1.
#' @param k numeric, represents the K-th smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @return A list with a random sample of order statistics from a generalized log-gamma distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
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

order_glg <- function(size,mu,sigma,lambda,k,n,alpha=0.05){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qglg(initial,mu,sigma,lambda)
  pdf    <- factorial(size)*cumprod(dglg(sample,mu,sigma,lambda))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
