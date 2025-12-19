#' vcov.sglg
#'
#' vcov.sglg extracts the estimated variance and covariance matrix associated to parameters from a model from an object of class 'sglg'.
#' @param object an object of the class sglg. This object is returned from the call to glg(), sglg(), survglg() or ssurvglg() function.
#' @param ... other arguments.
#' @export
vcov.sglg <- function(object, ...) {
    return(object$vcov)
}
