% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parModel.R
\name{PparLapply}
\alias{PparLapply}
\title{uses parLapply of library parallel for using workers on 
machines in a local network.}
\usage{
PparLapply(pop, EvalGene, lF)
}
\arguments{
\item{pop}{a population of genes.}

\item{EvalGene}{the function for evaluating a gene.}

\item{lF}{the local function factory which provides
all functions needed in \code{EvalGene}.}
}
\value{
Fitness vector.
}
\description{
uses parLapply of library parallel for using workers on 
machines in a local network.
}
\section{Warning}{


This section has not been properly tested.
Random number generation?
Examples?
}

\examples{
parm<-function(x) {function() {x}}
pop<-xegaInitPopulation(1000, lFxegaGaGene)
library(parallel)
clus<-makeCluster(2)
lFxegaGaGene$cluster<-parm(clus)
popnew<-PparLapply(pop, lFxegaGaGene$EvalGene, lFxegaGaGene)
stopCluster(clus)

}
\seealso{
Other Execution Model: 
\code{\link{MClapply}()},
\code{\link{MClapplyHet}()},
\code{\link{PparLapplyHet}()},
\code{\link{futureLapply}()},
\code{\link{futureLapplyHet}()}
}
\concept{Execution Model}
