% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shape.R
\name{shape}
\alias{shape}
\title{shape}
\usage{
shape(formula, npc, data, interval, semi, step)
}
\arguments{
\item{formula}{a symbolic description of the systematic component of the model to be fitted.}

\item{npc}{a data frame with potential nonparametric variables of the systematic part of the model to be fitted.}

\item{data}{a data frame which contains the variables in the model.}

\item{interval}{an optional numerical vector of length 2. In this interval is the maximum likelihood estimate of the shape parameter of the model.
By default is [0.1,1.5].}

\item{semi}{a logical value. TRUE means that the model has a non-parametric component. By default is FALSE.}

\item{step}{an optional positive value. This parameter represents the length of the step of the partition of the interval parameter.
By default is 0.1.}
}
\description{
Tool that supports the estimation of the shape parameter in semi-parametric or multiple linear accelerated failure time model with generalized log-gamma errors
under the presence of censored data. The estimation is based on the profiled likelihood function for the shape parameter of the model.
}
\examples{
rows  <- 200
columns <- 2
t_beta  <- c(0.5, 2)
t_sigma <- 1
t_lambda <- 1
set.seed(8142031)
x1 <- rbinom(rows, 1, 0.5)
x2 <- runif(columns, 0, 1)
X <- cbind(x1,x2)
s         <- t_sigma^2
a         <- 1/s
t_ini1    <- exp(X \%*\% t_beta) * rgamma(rows, scale = s, shape = a)
cens.time <- rweibull(rows, 0.3, 14)
delta     <- ifelse(t_ini1 > cens.time, 1, 0)
obst1 = t_ini1
for (i in 1:rows) {
if (delta[i] == 1) {
   obst1[i] = cens.time[i]
  }
}
example <- data.frame(obst1,delta,X)
lambda <- shape(Surv(log(obst1),delta) ~ x1 + x2 - 1, data=example)
lambda
# To change interval or step or both options
lambda <- shape(Surv(log(obst1),delta) ~ x1 + x2 - 1, data=example, interval=c(0.95,1.3), step=0.05)
lambda
}
\references{
Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
}
