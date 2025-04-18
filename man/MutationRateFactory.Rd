% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptivityMutation.R
\name{MutationRateFactory}
\alias{MutationRateFactory}
\title{Configure the mutation rate function of a genetic algorithm.}
\usage{
MutationRateFactory(method = "Const")
}
\arguments{
\item{method}{A string specifying a function for the mutation rate.}
}
\value{
A mutation rate function.
}
\description{
The \code{MutationRateFactory()} implements selection
             of one of the crossover rate functions in this
             package by specifying a text string.
             The selection fails ungracefully (produces
             a runtime error), if the label does not match.
             The functions are specified locally.

             Current support:

             \enumerate{
             \item "Const" returns \code{ConstMRate()} (Default).
             \item "IV" returns \code{IAMrate()}.
                   This function gives bad genes a higher mutation rate.
             }
}
\examples{
f<-MutationRateFactory("Const")
f(10, list(MutationRate1=function() {0.2}))
}
\seealso{
Other Configuration: 
\code{\link{AcceptFactory}()},
\code{\link{ApplyFactory}()},
\code{\link{CoolingFactory}()},
\code{\link{CrossRateFactory}()},
\code{\link{TerminationFactory}()},
\code{\link{checkTerminationFactory}()},
\code{\link{xegaConfiguration}()},
\code{\link{xegaEvalPopulationFactory}()}
}
\concept{Configuration}
