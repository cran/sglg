#' Quantile function for a generalized log-gamma variable
#'
#' \code{qglg} is used to calculate the quantile function of a generalized log-gamma variable at x.
#' @param x numeric, a real number between 0 and 1.
#' @param location numeric, represent the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param scale numeric, represent the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param shape numeric, represent the shape parameter of a generalized log-gamma distribution. Default value is 1.

#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' x <- runif(3,0,1)
#' qglg(sort(x),location=0,scale=1,shape=1)
#' @export qglg

qglg = function(x, location, scale, shape) {
  if(missingArg(x))
    return("The real value x is missing!")
  if(missingArg(location))
    location <- 0
  if(missingArg(scale))
    scale <- 1
  if(missingArg(shape))
    shape <- 1

  base_qglg <- function(p, location, scale, shape){
    out <- (1/abs(shape))*log(((shape^2)/2)*qchisq(p,2/shape^2))
    return(out)
  }
  v_qglg <- Vectorize(base_qglg,vectorize.args = c("p","location"))
  return(v_qglg(x, location, scale, shape))
}
