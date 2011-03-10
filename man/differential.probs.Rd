\name{differential.probs}
\alias{differential.probs}
\title{
  Calculating probabilities of differential expression in perturbation experiments
}
\description{
  Returns probabilities of differential expression for genes under
  perturbation of a set of regulators.
  Takes as input perturbation data and beliefs about known genes.
}
\usage{
differential.probs(data, beliefs, verbose, plot.it)
}
\arguments{
  \item{data}{
    A matrix of log expression ratios perrturbation vs control, for the
    genes (rows), in the perturbations of the regulators (columns). The data has to have row and colnames specified by the user.
  }
  \item{beliefs}{
    A list with names being a subset of the regulators (i.e., the names of \code{beliefs} have to be a subset of the columns of the \code{data}). 
    Each list entry for a given regulator is a matrix with rows corresponding to the genes that are known to respond in some way to the perturbation of
    this regulator. The rownames of the matrix must be a subset of the rows in the \code{data}.
The matrix can have either two or three columns.
    Each row  is a distribution over the differential and unchanged
    cluster (2 columns) or over down, up-regulated and unchanged cluster
    of genes  (3 columns). This distribution reflects the certainties with which a gene that corresponds to this row belongs to each of those clusters.   }
  \item{verbose}{When TRUE, for each regulator and its perturbation data, the execution prints out the parameters of the fitted model(s), indicating which components are differential and which are unchanged.} 
   \item{plot.it}{When TRUE, for each regulator and its perturbation data, the execution plots the Gaussian components of the fitted model(s), indicating which components are differential and which are unchanged.}  
}
\details{
 For each regulator, a belief-based mixture model is fitted to the observations in the \code{data}. The fitted models have the  number of model components equal to the number of columns in the corresponding beliefs.
 If no beliefs are given, unsupervised two-component mixture modeling is applied.
}
\value{ A matrix with columns for the regulators, rows for the genes,
  and entries giving the signed probabilities of differential expression.
}
\references{
  http://joda.molgen.mpg.de
}
\author{
Ewa Szczurek
}

\seealso{
 \code{\link{regulation.scores}}, \code{\link{deregulation.scores}} 
}
\examples{
	
data(damage)
	
# Get the probabilities of differential expression
# for the knockout of p53 in healthy cells
probs.healthy.p53= differential.probs(data.healthy[,"p53",FALSE],
					beliefs.healthy["p53"], TRUE,TRUE)
	
# Get the probabilities of differential expression
# for the knockout of Ste12 under pheromone treatment 	
library(bgmm)
data(Ste12)
data=as.matrix(Ste12Data)
colnames(data)="Ste12"
beliefs=list(Ste12=Ste12Beliefs)
diff.p=differential.probs(data,beliefs,TRUE,TRUE)
	
\dontrun{
probs.healthy= differential.probs(data.healthy, beliefs.healthy, TRUE,TRUE)
}
}
