\name{regulation.scores}
\alias{regulation.scores}
\title{
  Calculating regulation scores.
}
\description{
  Regulation scores reflect the actual regulatory  influence of the regulators on the genes. For each regulator, these  are probabilities of  differential expression averaged over all perturbation experiments that affect the regulator.
}
\usage{regulation.scores(probs, model, verbose)
}

\arguments{
  \item{probs}{
    A matrix of probabilities of differential expression of the genes
    (rows) under perturbations of regulators (columns). Obtained with the \code{\link{differential.probs}} function
  }
  \item{model}{
    A pathway \code{model} is a matrix with rows and columns equal to the names of the regulators
    (i.e., the columns of probs).  Each \code{model} has an entry 1 where the
    regulator in the corresponding row influences the regulator in the
    corresponding column. A transitive closure of the input \code{model} is
    computed to get information about the experiments affecting each
    regulator. 
  }
  \item{verbose}{When TRUE, the execution prints informative messages}
}
\value{ A matrix with columns for the regulators, rows for the genes,
  and entries giving the regulation scores.
}
\references{
  http://joda.molgen.mpg.de
}
\author{
Ewa Szczurek
}
\seealso{
 \code{\link{differential.probs}}, \code{\link{deregulation.scores}} 
}
\examples{
data(damage)
	
# Get the probabilities of differential expression
# for the knockouts of ATM, RelA and p53 in healthy cells
probs.healthy= differential.probs(data.healthy, beliefs.healthy)
#Get the regulation scores for ATM, RelA and for p53
regulation.healthy= regulation.scores(probs.healthy, model.healthy, TRUE)
}
