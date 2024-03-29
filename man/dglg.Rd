% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dglg.R
\name{dglg}
\alias{dglg}
\title{Density distribution function for a generalized log-gamma variable}
\usage{
dglg(x, location, scale, shape)
}
\arguments{
\item{x}{numeric, a real number.}

\item{location}{numeric, represent the location parameter of a generalized log-gamma distribution. Default value is 0.}

\item{scale}{numeric, represent the scale parameter of a generalized log-gamma distribution. Default value is 1.}

\item{shape}{numeric, represent the shape parameter of a generalized log-gamma distribution. Default value is 1.}
}
\description{
\code{dglg} is used to calculate the density distribution function of a generalized log-gamma variable at x.
}
\examples{
x <- seq(-4,4,length=100)
dglg(x,location=0,scale=1,shape=1)
plot(x,dglg(x,location=0,scale=1,shape=1),type="l",xlab="x",ylab="Density")
}
\references{
Carlos Alberto Cardozo Delgado, Semi-parametric generalized log-gamma regression models. Ph. D. thesis. Sao Paulo University.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>
}
