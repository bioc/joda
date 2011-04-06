\name{deregulation.p.values}
\alias{deregulation.p.values}
\title{
  Calculating deregulation p-values using resampling method.
}
\description{
  Deregulation p-values based on deregulation scores.
  They are calculated as fraction of permutations that give more extreme deregulation scores than for original data.}
\usage{deregulation.p.values(data.1, beliefs.1, model.1, data.2, beliefs.2, model.2, N=100, verbose=FALSE)}

\arguments{
  \item{data.1, data.2}{Matrices of log expression ratios perturbation vs control, for the genes (rows), in the perturbations of the regulators (columns).
  See \code{\link{differential.probs}} for more details.
  }
  \item{beliefs.1, beliefs.2}{Lists of beliefs. See \code{\link{differential.probs}} for more details.
  }
  \item{model.1, model.2}{Pathway topologies. See \code{\link{differential.probs}} for more details.
  }
  \item{N}{A number of replications used to calculate p-values
  }
  \item{verbose}{When TRUE, the execution prints informative messages}
}
\details{The deregulation p-values are calculated as fraction of permutations that give more extreme deregulation scores than for original data. }
\value{ A list with two matrices. This p-values in the slot \code{deregulation.p.values} and with the original deregulation scores in the slot  \code{deregulationOrg}.
}
\references{
  http://joda.molgen.mpg.de
}
\author{
Ewa Szczurek
}

\seealso{
 \code{\link{differential.probs}}, \code{\link{regulation.scores}}, \code{\link{regulation.scores}}
}
\examples{
\dontrun{
# Step 1
library(joda)
data(damage)

deregulationObj = deregulation.p.values(data.healthy, beliefs.healthy, model.healthy, data.damage, beliefs.damage, model.damage, N=100, verbose=TRUE) 
boxplot(deregulationObj$deregulation.p.values)
}
}
