#' Density distribution function for a generalized log-gamma variable
#'
#' \code{dglg} is used to calculate the density distribution function of a generalized log-gamma variable at x.
#' @param x numeric, a real number.
#' @param location numeric, represent the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param scale numeric, represent the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param shape numeric, represent the shape parameter of a generalized log-gamma distribution. Default value is 1.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' x <- runif(60,-3,3)
#' dglg(sort(x),location=0,scale=1,shape=1)
#' plot(sort(x),dglg(sort(x),location=0.5,scale=1,shape=1),type="l",xlab="x",ylab="Density")
#' @export dglg

dglg = function(x, location, scale, shape) {
  if(missingArg(x))
    return("The real value x is missing!")
  if(missingArg(location))
    location <- 0
  if(missingArg(scale))
    scale <- 1
  if(missingArg(shape))
    shape <- 1

  c_l <- function(shape) {
    if (abs(shape) < 0.085) {
      if (shape > 0) {
        shape <- 0.085
      }
      if (shape < 0) {
        shape <- -0.085
      }
    }
    invlambdos <- 1/shape^2
    c <- abs(shape)/gamma(invlambdos)
    output <- c * (invlambdos^invlambdos)
    return(output)
  }

  base_dglg <- function(x, location, scale, shape){
    y <- (x - location)/scale
    out <- (1/shape)*y -(1/shape^2)*exp(shape*y)
    out <- exp(out)
    out <- c_l(shape)*out/scale
    return(out)
  }
  v_dglg <- Vectorize(base_dglg)
  return(v_dglg(x, location, scale, shape))
}
