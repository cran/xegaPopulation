% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cooling.R
\name{ExponentialMultiplicativeCooling}
\alias{ExponentialMultiplicativeCooling}
\title{Exponential multiplicative cooling.}
\usage{
ExponentialMultiplicativeCooling(k, lF)
}
\arguments{
\item{k}{Number of steps to discount.}

\item{lF}{Local configuration.}
}
\value{
Temperature at time k.
}
\description{
The temperature at time k is the net present value 
             of the starting temperature. The discount factor 
             is \code{lF$Alpha()}. 
             \code{lF$Alpha()} should be in \code{[0, 1]}.
}
\details{
Temperature is updated at the end of each generation
         in the main loop of the genetic algorithm.
         \code{lF$Temp0()} is the starting temperature.
         \code{lF$Alpha()} is the discount factor.
}
\examples{
parm<-function(x){function() {return(x)}}
lF<-list(Temp0=parm(100), Alpha=parm(0.99))
ExponentialMultiplicativeCooling(0, lF)
ExponentialMultiplicativeCooling(2, lF)
}
\references{
Kirkpatrick, S., Gelatt, C. D. J, and Vecchi, M. P. (1983):
     Optimization by Simulated Annealing. 
     Science, 220(4598): 671-680.
     <doi:10.1126/science.220.4598.671>
}
\seealso{
Other Cooling: 
\code{\link{ExponentialAdditiveCooling}()},
\code{\link{LogarithmicMultiplicativeCooling}()},
\code{\link{PowerAdditiveCooling}()},
\code{\link{PowerMultiplicativeCooling}()},
\code{\link{TrigonometricAdditiveCooling}()}
}
\concept{Cooling}
