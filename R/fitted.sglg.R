#' Extract Fitted Values
#'
#' fitted.sglg extracts the fitted values from a model from an object of class 'sglg'.
#' @param object an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg() function.
#' @param ... other arguments.
#' @export
fitted.sglg <- function(object, ...) {
  return(object$y_est)
}
