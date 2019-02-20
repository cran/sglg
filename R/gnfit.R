#' gnfit
#'
#' This function provides some useful statistics to assess the quality of fit of generalized log-gamma probabilistic model, including the statistics Cramer-von Mises and Anderson-Darling. It can also calculate other goodness of fit such as Hannan-Quin Information Criterion and Kolmogorov-Smirnov test.

#' @param starts numeric vector. Initial parameters to maximize the likelihood function
#' @param data numeric vector. A sample of a generalized log-gamma distribution.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' set.seed(1)
#' sample <- rglg(16,location=0.7,scale=0.5,shape=1.2)
#' gnfit(c(0.55,0.6,1),sample)
#' @importFrom AdequacyModel goodness.fit
#' @export gnfit

gnfit <- function(starts,data){
         pdf_glg <- function(par,x){
                    a <- par[1]
                    b <- par[2]
                    c <- par[3]
                    dloggamma(x, mu = a, sigma = b, lambda= c)
         }
         cdf_glg <- function(par,x){
                    a <- par[1]
                    b <- par[2]
                    c <- par[3]
                    ploggamma(x, mu = a, sigma = b, lambda = c)
         }
         output <- suppressWarnings(goodness.fit(pdf = pdf_glg, cdf = cdf_glg, starts = c(0,1,1), data = data,
                        method = "PSO", lim_inf = c(-20,0,-5), lim_sup = c(20,10,5),mle = NULL,domain=c(-Inf,Inf)))
         x = seq(floor(min(data)),ceiling(max(data)), length.out = 500)
         hist(data, probability = TRUE,main="")
         lines(x, pdf_glg(par = output$mle, x), col = "red")
         return(list(CM=output$W,AD=output$A,KS=output$KS,HQIC=output$HQIC,BIC=output$BIC,AIC=output$AIC,MLE=output$mle))
}
