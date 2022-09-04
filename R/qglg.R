#' Quantile function for a generalized log-gamma distribution
#'
#' \code{qglg} is used to calculate the quantile function of a generalized log-gamma variable at x.
#' @param x numeric, a vector with values between 0 and 1.
#' @param location numeric, represents the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param scale numeric, represents the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param shape numeric, represents the shape parameter of a generalized log-gamma distribution. Default value is 1.
#' @return A vector with the same size of x with the quantile values of a generalized log-gamma distribution.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' # Calculating the quartiles of a glg(0,1,-1) distribution
#' x <- c(0.25, 0.5, 0.75)
#' qglg(x, location = 0, scale = 1, shape = -1)
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
    if (shape < 0) p <- 1 - p
    out <-  (1/shape) *log((0.5 * (shape^2)) * qchisq(p, 2/shape^2))
    out <- location + scale * out
    return(out)
  }
  v_qglg <- Vectorize(base_qglg,vectorize.args = c("p","location"))
  return(v_qglg(x, location, scale, shape))
}
