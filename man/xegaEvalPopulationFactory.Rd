% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xegaEvalPopulation.R
\name{xegaEvalPopulationFactory}
\alias{xegaEvalPopulationFactory}
\title{Configures the evaluation of the population of
a genetic algorithm.}
\usage{
xegaEvalPopulationFactory(method = "EvalPopulation")
}
\arguments{
\item{method}{A string specifying the termination condition.}
}
\value{
A boolean function implementing the termination condition.
}
\description{
\code{xegaEvalPopulationFactory()} implements the selection
             of the evaluation function for the population of a 
             genetic algorithm.

Current support:

\enumerate{
  \item "EvalPopulation" returns 
       \code{xegaEvalPopulation}. (Default)
  \item "RepEvalPopulation" returns 
       \code{xegaReplEvalPopulation}. 
       For stochastic functions.
       Needs \code{lF$rep()} for the number of repetitions
       and \code{lF$apply()} for the (parallel) apply function.
   }
}
\seealso{
Other Configuration: 
\code{\link{AcceptFactory}()},
\code{\link{ApplyFactory}()},
\code{\link{CoolingFactory}()},
\code{\link{CrossRateFactory}()},
\code{\link{MutationRateFactory}()},
\code{\link{TerminationFactory}()},
\code{\link{checkTerminationFactory}()},
\code{\link{xegaConfiguration}()}
}
\concept{Configuration}
