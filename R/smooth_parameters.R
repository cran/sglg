#' smoothp
#'
#' Tool that supports the selection of the smoothing parameters in semi-parametric generalized log-gamma models.
#' The selection is based on the AIC, BIC, or Generalized Cross Validation methods.

#' @param formula a symbolic description of the systematic component of the model to be fitted.
#' @param npc a data frame with potential nonparametric variables of the systematic part of the model to be fitted.
#' @param data a data frame which contains the variables in the model.
#' @param method There are three possible criteria to estimate the smoothing parameters: Penalized Akaike Criterion 'PAIC', Penalized Bayesian Criterion 'PBIC' and
#' Generalized Cross Validation 'GCV'. The default method is 'PAIC'.
#' @param basis a name of the cubic spline basis to be used in the model. Supported basis include deBoor and Gu basis.
#' @param interval an optional numerical vector of length 2. In this interval is the maximum likelihood estimate of the shape parameter of the model.
#' By default is [0.5,2].
#' @param step an optional positive value. This parameter represents the length of the step of the partition of the interval parameter.
#' By default is 0.5.
#' @references Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph.D. thesis. Sao Paulo University.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
#' @examples
#' set.seed(1)
#' rows<- 150
#' t_beta <- c(0.5,2)
#' t_sigma <- 0.5
#' t_lambda <- 1
#' x1 <- runif(rows,-3,3)
#' x2 <- rbinom(rows,1,0.5)
#' X <- cbind(x1,x2)
#' t <- as.matrix((2*1:rows - 1)/(2*rows))
#' colnames(t) <- "t"
#' f_t <- cos(4*pi*t)
#' error <- rglg(rows,0,1,t_lambda)
#' y <- X %*%t_beta + f_t + t_sigma*error
#' colnames(y) <- "y"
#' data <- data.frame(y,X,t)
#' fit1 <- sglg(y ~ x1 + x2 - 1,npc=t,data=data,basis = "deBoor",alpha0=1)
#' fit1$AIC
#' smoothp(y ~ x1 + x2 - 1,npc=t,data=data,basis = "deBoor")
#' fit2 <- sglg(y ~ x1 + x2 - 1,npc=t,data=data,basis = "Gu",alpha0=0.5)
#' fit2$BIC
#' smoothp(y ~ x1 + x2 - 1,npc=t,data=data,basis = "Gu")
#' #################################################
#' # An example with two non-parametric components #
#' #################################################
#' set.seed(2)
#' t_2 <- as.matrix(rnorm(rows,sd=0.5))
#' colnames(t_2) <- 't_2'
#' f_t_2 <- exp(t_2)
#' error <- rglg(rows,0,1,t_lambda)
#' y_2 <- X %*%t_beta + f_t + f_t_2 + t_sigma*error
#' colnames(y_2) <- 'y_2'
#' data2 <- data.frame(y_2,X,t,t_2)
#' npcs <- cbind(t,t_2)
#' smoothp(y ~ x1 + x2 - 1,npc=npcs,data=data, method='GCV')
#' @importFrom plot3D scatter3D
#' @importFrom TeachingSampling SupportWR
#' @export smoothp

smoothp = function(formula, npc, data, method='PAIC', basis, interval, step){
  if (missingArg(interval))
    interval <- c(0.5, 2)
  if (missingArg(step))
    step <- 0.5

  k <- dim(npc)[2]
  if(k>2)
    stop('This function only considers one or two non-parametric components!')

  if (missingArg(basis))
    basis <- rep("deBoor", k)

  N <- length(seq(interval[1],interval[2],by=step))
  alphas <-  SupportWR(N,k,ID=seq(interval[1],interval[2],by=step))
  colnames(alphas) <- paste("alpha_",1:k,sep="")

  if(k == 1){
    if(method=='PAIC'){
      values_AIC <- vector()
      J <- dim(alphas)[1]
      for(j in 1:J){
        values <- try(sglg(formula,npc,data=data,basis,alpha0=alphas[j,],format='simple'),silent=TRUE)
        if(is.list(values)){
          values_AIC[j] <- values$AIC
        }
      }
      min_AIC <- min(values_AIC,na.rm=TRUE)
      index_min_AIC <- which.min(values_AIC)
      min_alpha_AIC <- alphas[index_min_AIC,]
      output_1 <- data.frame(alphas,values_AIC)
      output_AIC <- data.frame(min_alpha_AIC,min_AIC)
      plot1 <- ggplot(data=output_1, aes(alphas,values_AIC)) +
                geom_point(colour="orange",alpha=1,size=1.25) +
                xlim(range(alphas)) +
                geom_point(data=output_AIC,aes(min_alpha_AIC,min_AIC), colour="darkblue", alpha=0.25, size=3.5) +
                labs(x = "Smooth parameters", y = "AICp") +
                ggtitle("Penalized Akaike Information Criterion")
      grid.arrange(plot1)
      return(list(output_1,optimum_AIC = output_AIC))
     }
    if(method=='PBIC'){
       values_BIC <- vector()
       J <- dim(alphas)[1]
       for(j in 1:J){
         values <- try(sglg(formula,npc,data=data,basis,alpha0=alphas[j,],format='simple'),silent=TRUE)
         if(is.list(values)){
           values_BIC[j] <- values$BIC
         }
       }
       min_BIC <- min(values_BIC,na.rm=TRUE)
       index_min_BIC <- which.min(values_BIC)
       min_alpha_BIC <- alphas[index_min_BIC,]
       output_1 <- data.frame(alphas,values_BIC)
       output_BIC <- data.frame(min_alpha_BIC,min_BIC)
       plot2 <- ggplot(data=output_1, aes(alphas,values_BIC)) +
           geom_point(colour="orange",alpha=1,size=1.25) +
           xlim(range(alphas)) +
           geom_point(data=output_BIC,aes(min_alpha_BIC,min_BIC), colour="blue", alpha=0.25, size=3.5) +
           labs(x = "Smooth parameters", y = "BICp") +
           ggtitle("Penalized Bayesian Information Criterion")
       grid.arrange(plot2)
    return(list(output_1,optimum_BIC = output_BIC))
    }
    if(method=='GCV'){
      values_GCV <- vector()
      J <- dim(alphas)[1]
      for(j in 1:J){
        values <- try(sglg(formula,npc,data=data,basis,alpha0=alphas[j,],format='simple'),silent=TRUE)
        if(is.list(values)){
          n <- length(values$y)
          values_GCV[j] <- n*sum((values$y - values$y_est)^2)/(n - (values$df + 2))^2
        }
      }
      min_GCV <- min(values_GCV,na.rm=TRUE)
      index_min_GCV <- which.min(values_GCV)
      min_alpha_GCV <- alphas[index_min_GCV,]
      output_1 <- data.frame(alphas,values_GCV)
      output_GCV <- data.frame(min_alpha_GCV,min_GCV)
      plot3 <- ggplot(data=output_1, aes(alphas,values_GCV)) +
               geom_point(colour="orange",alpha=1,size=1.25) +
               xlim(range(alphas)) +
               geom_point(data=output_GCV,aes(min_alpha_GCV,min_GCV), colour="blue", alpha=0.25, size=3.5) +
               labs(x = "Smooth parameters", y = "GCV") +
               ggtitle("Generalized Cross Validation Criterion")
      grid.arrange(plot3)
    return(list(output_1,optimum_GCV = output_GCV))
    }
  }
 else{
   if(method=='PAIC'){
     values_AIC <- vector()
     J <- dim(alphas)[1]
     for(j in 1:J){
       values <- try(sglg(formula,npc,data=data,basis,alpha0=alphas[j,],format='simple'),silent=TRUE)
       if(is.list(values)){
         values_AIC[j] <- values$AIC
       }
     }
     min_AIC <- as.matrix(min(values_AIC,na.rm=TRUE))
     colnames(min_AIC) <- "PAIC"
     index_min_AIC <- which.min(values_AIC)
     min_alpha_AIC <- alphas[index_min_AIC,]
     output_1 <- data.frame(alphas,values_AIC)
     output_AIC <- c(min_alpha_AIC,min_AIC)
     scatter3D(x = output_1$alpha_1, y = output_1$alpha_2, z = output_1$values_AIC,
               bty = 'g', pch = 20, cex = 1.5, cex.axis = 0.65, cex.lab = 0.6, ticktype = "detailed",
               xlab = "alpha_1", ylab = "alpha_2", clab = "P_AIC",
               alpha=0.6)
     return(list(output_1,optimum_AIC = output_AIC))
   }
   if(method=='PBIC'){
     values_BIC <- vector()
     J <- dim(alphas)[1]
     for(j in 1:J){
       values <- try(sglg(formula,npc,data=data,basis,alpha0=alphas[j,],format='simple'),silent=TRUE)
       if(is.list(values)){
         values_BIC[j] <- values$BIC
       }
     }
     min_BIC <- min(values_BIC,na.rm=TRUE)
     index_min_BIC <- which.min(values_BIC)
     min_alpha_BIC <- alphas[index_min_BIC,]
     output_1 <- data.frame(alphas,values_BIC)
     output_BIC <- c(min_alpha_BIC,min_BIC)
     scatter3D(x = output_1$alpha_1, y = output_1$alpha_2, z = output_1$values_BIC,
               bty = 'g', pch = 20, cex = 1.5, cex.axis = 0.65, cex.lab = 0.6, ticktype = "detailed",
               xlab = "alpha_1", ylab = "alpha_2", clab = "P_BIC",
               alpha=0.6)
     return(list(output_1,optimum_BIC = output_BIC))
   }
   if(method=='GCV'){
     values_GCV <- vector()
     J <- dim(alphas)[1]
     for(j in 1:J){
       values <- try(sglg(formula,npc,data=data,basis,alpha0=alphas[j,],format='simple'),silent=TRUE)
       if(is.list(values)){
         n <- length(values$y)
         values_GCV[j] <- n*sum((values$y - values$y_est)^2)/(n - (values$df + 2))^2
       }
     }
     min_GCV <- min(values_GCV,na.rm=TRUE)
     index_min_GCV <- which.min(values_GCV)
     min_alpha_GCV <- alphas[index_min_GCV,]
     output_1 <- data.frame(alphas,values_GCV)
     output_GCV <- c(min_alpha_GCV,min_GCV)
     scatter3D(x = output_1$alpha_1, y = output_1$alpha_2, z = output_1$values_GCV,
               bty = 'g', pch = 20, cex = 1.5, cex.axis = 0.65, cex.lab = 0.6, ticktype = "detailed",
               xlab = "alpha_1", ylab = "alpha_2", clab = "GCV",
               alpha=0.6)
     return(list(output_1,optimum_GCV = output_GCV))
   }

  }
}
