#' Random number generation for a generalized log-gamma distribution
#'
#' \code{rglg} is used to generate random numbers for a generalized log-gamma distribution.
#' @param n numeric, size of the random sample.
#' @param location numeric, represents the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param scale numeric, represents the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param shape numeric, represents the shape parameter of a generalized log-gamma distribution. Default value is 1.

#' @return A vector of size n with the generalized log-gamma random values.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' u <- rglg(100, location = 0, scale = 1, shape = -1)
#' @export rglg
rglg = function(n, location, scale, shape) {
    if (missingArg(location))
        location <- 0
    if (missingArg(scale))
        scale <- 1
    if (missingArg(shape))
        shape <- 1
    quantiles <- runif(n, 0, 1)
    pQ <- matrix(0, n, 1)

    for (i in 1:n) {
        pQ[i] <- (1/shape) * log((0.5 * shape^2) * qchisq(quantiles[i], 2/shape^2))
        pQ[i] <- location + scale * pQ[i]
    }
    return(pQ)
}
