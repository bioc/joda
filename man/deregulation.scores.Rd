\name{deregulation.scores}
\alias{deregulation.scores}
\title{
  Calculating deregulation scores.
}
\description{
  Deregulation scores  quantify the extent to which the regulatory
  eflect of each regulator changes between the two compared cell populations.}
\usage{deregulation.scores(reg.scores1, reg.scores2,verbose)}

\arguments{
  \item{reg.scores1}{
    A matrix of regulation scores of the genes  (rows) for the
    regulators (columns), compued with the \code{\link{regulation.scores}} function. Given for the first cell population.
  }
  \item{reg.scores2}{
    The same as reg.scores1 but given for the second cell population.
  }
  \item{verbose}{When TRUE, the execution prints informative messages}
}
\details{The deregulation scores are computed by subtracting reg.scores1 from  reg.scores2.}
\value{ A matrix with columns for the regulators, rows for the genes,
  and entries giving the deregulation scores.
}
\references{
  http://joda.molgen.mpg.de
}
\author{
Ewa Szczurek
}

\seealso{
 \code{\link{differential.probs}}, \code{\link{regulation.scores}} 
}
\examples{
data(damage)
	
# Step 1
# Get the probabilities of differential expression
# for the knockout of ATM in the healthy cells
probs.healthy.ATM= differential.probs(data.healthy[,"ATM",FALSE], NULL)
	
# Get the probabilities of differential expression
# for the knockout of ATM in the damaged cells
probs.damage.ATM= differential.probs(data.damage[,"ATM",FALSE], NULL)
	
# Step 2 
# Regulation scores for a dataset with only one regulator
# equal the signed probabilities
	
# Step 3 
# Get the deregulation scores
deregulation.ATM= deregulation.scores(probs.healthy.ATM, probs.damage.ATM, TRUE)

\dontrun{
# Step 1
probs.healthy= differential.probs(data.healthy, beliefs.healthy)
probs.damage= differential.probs(data.damage, beliefs.damage)

# Step 2
regulation.healthy= regulation.scores(probs.healthy, model.healthy)
regulation.damage= regulation.scores(probs.damage, model.damage)

# Step 3
deregulation= deregulation.scores(regulation.healthy, regulation.damage, TRUE)
}
}
