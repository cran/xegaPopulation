% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xegaPopulation.R
\name{xegaObservePopulation}
\alias{xegaObservePopulation}
\title{Observe summary statistics of the fitness of the population.}
\usage{
xegaObservePopulation(fit, v = vector())
}
\arguments{
\item{fit}{Vector of fitness values of a population.}

\item{v}{Vector of population statistic vectors.}
}
\value{
Vector of population statistics. If position
        \code{x} modulo \code{8} equals
        \enumerate{
        \item \code{1}:  Mean fitness.
        \item \code{2}:  Min fitness.
        \item \code{3}:  Lower-hinge 
                       (approx. 1st quartile) of fitness.
        \item \code{4}:  Median fitness.
        \item \code{5}:  Upper-hinge 
                        (approx. 3rd quartile) of fitness.
        \item \code{6}:  Max fitness.
        \item \code{7}:  Variance.
        \item \code{8}: Mean absolute deviation. 
            }
}
\description{
\code{xegaObservePopulation()} reports
             summary statistics of the fitness of the population.
}
\details{
Population statistics are used for 
             \itemize{
                \item implementing individually variable operator rates and
                \item visualizing the progress of the algorithm.      
                }
}
\examples{
pop10<-xegaInitPopulation(10, lFxegaGaGene)
epop10<-xegaEvalPopulation(pop10, lFxegaGaGene)
popStats<-xegaObservePopulation(epop10$fit)
popStats<-xegaObservePopulation(epop10$fit, popStats)
matrix(popStats, ncol=8, byrow=TRUE)

}
\seealso{
Other Population Layer: 
\code{\link{xegaBestGeneInPopulation}()},
\code{\link{xegaBestInPopulation}()},
\code{\link{xegaEvalPopulation}()},
\code{\link{xegaInitPopulation}()},
\code{\link{xegaLogEvalsPopulation}()},
\code{\link{xegaNextPopulation}()},
\code{\link{xegaRepEvalPopulation}()},
\code{\link{xegaSummaryPopulation}()}
}
\concept{Population Layer}
