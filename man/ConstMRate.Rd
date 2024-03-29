% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptivityMutation.R
\name{ConstMRate}
\alias{ConstMRate}
\title{Constant mutation rate.}
\usage{
ConstMRate(fit, lF)
}
\arguments{
\item{fit}{Fitness of gene.}

\item{lF}{Local configuration.}
}
\value{
Constant mutation rate.
}
\description{
Constant mutation rate.
}
\examples{
parm<-function(x){function() {return(x)}}
lF<-list()
lF$MutationRate1<-parm(0.20)
ConstMRate(100, lF)
ConstMRate(50, lF)
}
\seealso{
Other Rates: 
\code{\link{ConstCRate}()}
}
\concept{Rates}
