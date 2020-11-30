#' Survival, Hazard and Cumulative Hazard functions for a Generalized Gamma Distribution
#'
#'\code{survival_gg} is used to obtain the value of survival, hazard and cumulative hazard functions of a generalized gamma distribution at a positive value.
#' @param x numeric, represent a positive value. Default value is 1.
#' @param mu numeric, represents the location parameter of a generalized gamma distribution. Default value is 0.
#' @param sigma numeric, represents the scale parameter of a generalized gamma distribution. Default value is 1.
#' @param lambda numeric, represents the shape parameter of a generalized gamma distribution. Default value is 1.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @references Jerald F. Lawless (2003). Statistical Models and Methods for Lifetime Data. Second Edition. John-Wiley & Sons
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' survival_gg(0.0001,0,1,-1) # Extreme value type I distribution, maximum case.
#'  times <- seq(0.05,7,by=0.05)
#' plot(times, survival_gg(times,0,1,-1)$survival_value,type='l')
#' plot(times, survival_gg(times,0,1,-1)$hazard_value,type='l')
#' plot(times, survival_gg(times,0,1,-1)$cumulative_hazard_value,type='l')
#' @export survival_gg
survival_gg = function(x, mu, sigma, lambda) {
  if (missingArg(x))
     x <- 1
  if (sum(x < 0)>0)
     stop('x is a positive value!')
  if (missingArg(mu))
    mu <- 0
  if (missingArg(sigma))
    sigma <- 1
  if (missingArg(lambda))
    lambda <- 1

  surv <- function(x){
           return(1 - pglg(log(x), mu, sigma, lambda))
  }
  hazard <- function(x){
            return(dglg(log(x), mu, sigma, lambda)/(x*surv(x)))
     }
  cumulative_hazard <- function(x){
                       output <- integrate(hazard,0,x)
                       return(output[1][[1]])
  }
  cumulative_hazard_v <- Vectorize(cumulative_hazard,vectorize.args = 'x')
  return(list(survival_value = surv(x), hazard_value = hazard(x), cumulative_hazard_value = cumulative_hazard_v(x)))
}
