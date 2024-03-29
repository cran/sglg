% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_glg.R
\name{survival_gg}
\alias{survival_gg}
\title{Survival, Hazard, and Cumulative Hazard functions for a Generalized Gamma Distribution}
\usage{
survival_gg(x, mu, sigma, lambda)
}
\arguments{
\item{x}{numeric, represent a vector of positive values. Default value is 1.}

\item{mu}{numeric, represents the location parameter of a generalized gamma distribution. Default value is 0.}

\item{sigma}{numeric, represents the scale parameter of a generalized gamma distribution. Default value is 1.}

\item{lambda}{numeric, represents the shape parameter of a generalized gamma distribution. Default value is 1.}
}
\value{
A list of three vectors, survival, hazard, and cumulative hazard values of a generalized gamma distribution.
}
\description{
\code{survival_gg} is used to obtain the value of survival, hazard and cumulative hazard functions of a generalized gamma distribution at a positive value.
}
\examples{
survival_gg(0.0001,0,1,1) # Extreme value type I distribution, maximum case.
times <- seq(0.05,7,by=0.05)
plot(times, survival_gg(times,0,1,1)$survival_value,type='l')
plot(times, survival_gg(times,0,1,1)$hazard_value,type='l')
plot(times, survival_gg(times,0,1,1)$cumulative_hazard_value,type='l')
}
\references{
Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph.D. thesis. Sao Paulo University.

Jerald F. Lawless (2003). Statistical Models and Methods for Lifetime Data. Second Edition. John-Wiley & Sons
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
}
