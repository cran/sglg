#' Density distribution function for a generalized log-gamma variable
#'
#' \code{dglg} is used to calculate the density distribution function of a generalized log-gamma variable at x.
#' @param x numeric, a real number.
#' @param location numeric, represent the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param scale numeric, represent the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param shape numeric, represent the shape parameter of a generalized log-gamma distribution. Default value is 1.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' x <- seq(-4,4,length=100)
#' dglg(x,location=0,scale=1,shape=1)
#' plot(x,dglg(x,location=0,scale=1,shape=1),type="l",xlab="x",ylab="Density")
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

  base_dglg <- function(x, location, scale, shape){
    y <- (x - location)/scale
    out <- (1/shape)*y -(1/shape^2)*exp(shape*y)
    out <- exp(out)
    out <- c_l(shape)*out/scale
    return(out)
  }
  v_dglg <- Vectorize(base_dglg,vectorize.args = c("x","location"))
  return(v_dglg(x, location, scale, shape))
}
