#' gnfit
#'
#' This function provides some useful statistics to assess the quality of fit of generalized log-gamma probabilistic model, including the statistics Cramer-von Mises and Anderson-Darling. It can also calculate other goodness of fit such as Hannan-Quin Information Criterion and Kolmogorov-Smirnov test.

#' @param starts numeric vector. Initial parameters to maximize the likelihood function
#' @param data numeric vector. A sample of a generalized log-gamma distribution.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>, G. Paula and L. Vanegas.
#' @examples
#' \dontrun{
#' set.seed(12)
#' sample <- rglg(100,location=0,scale=0.5,shape=0.75)
#' result <- gnfit(starts=c(0.1,0.75,1),data=sample)
#' result
#' }
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
         sample <- data
         output <- suppressWarnings(goodness.fit(pdf = pdf_glg, cdf = cdf_glg, starts = starts, data = sample,
                        method = "PSO", lim_inf = c(-20,0,-5), lim_sup = c(20,10,5),mle = NULL,domain=c(-Inf,Inf)))
         x = seq(floor(min(sample)) - 0.5, ceiling(max(sample)) + 0.5, length.out = 500)

         f_x <- pdf_glg(par = output$mle, x)
         fit <- as.data.frame(cbind(x,f_x))
         plot1 <- ggplot(data=as.data.frame(sample), aes(sample)) +
         geom_density(colour="orange",fill="orange",alpha=0.4,size=0.7) +
         xlim(range(x)) +
         geom_line(data=fit,aes(x,f_x),colour="blue",alpha=0.5,size=0.7) +
         labs(x = "The Parametric Estimated Density, PED, is the blue curve and the Non-Parametric Estimated Density, Non-PED, is the orange curve.", y = "Density", colour = "Parametric Est. Den") +
         theme(legend.position="top")+
         ggtitle("PED vs Non-PED")
         grid.arrange(plot1,ncol=1)

         return(list(CM=output$W,AD=output$A,KS=output$KS,HQIC=output$HQIC,BIC=output$BIC,AIC=output$AIC,MLE=output$mle))
}
