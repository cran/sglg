#' Cumulative distribution function for a generalized log-gamma variable
#'
#' \code{pglg} is used to calculate the cumulative distribution function of a generalized log-gamma variable at x.
#' @param x numeric, a vector of real values.
#' @param location numeric, represents the location parameter of a generalized log-gamma distribution. Default value is 0.
#' @param scale numeric, represents the scale parameter of a generalized log-gamma distribution. Default value is 1.
#' @param shape numeric, represents the shape parameter of a generalized log-gamma distribution. Default value is 1.
#' @return A vector with the same size of x with the cumulative probability values of a generalized log-gamma distribution.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' x <- runif(3,-1,1)
#' pglg(sort(x),location=0, scale=1, shape=1)
#' @import pracma
#' @export pglg

pglg = function(x, location, scale, shape) {
       if(missingArg(x))
          return("The real value x is missing!")
       if(missingArg(location))
          location <- 0
       if(missingArg(scale))
          scale <- 1
       if(missingArg(shape))
         shape <- 1

       base_pglg <- function(x, location, scale, shape){
                    y <- (x - location)/scale
                    out <- as.matrix(pracma::gammainc((1/shape^2)*exp(shape*y),1/shape^2))
                    rownames(out)[3] <- NA
                    output <- out[3]
                    if(shape < 0) output <- 1 - output
                    return(output)
       }
       v_pglg <- Vectorize(base_pglg,vectorize.args = c("x","location"))
       return(v_pglg(x, location, scale, shape))
}
