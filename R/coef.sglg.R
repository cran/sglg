#' coef.sglg
#'
#' coef.sglg extracts and display the estimated coefficients associated to the location parameter from a model from an object of class 'sglg'.
#' @param object an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg() function.
#' @param ... other arguments.
#' @export
coef.sglg <- function(object, ...) {
  if (object$semi == FALSE & object$censored == FALSE) {
    return(object$mu)
  }
  if (object$semi == TRUE) {
     return(object$mu[1:(object$p)])
  }
}
