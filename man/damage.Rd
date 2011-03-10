\name{damage}
\alias{data.healthy}
\alias{beliefs.healthy}
\alias{model.healthy}
\alias{data.damage}
\alias{beliefs.damage}
\alias{model.damage}
\docType{data}
\title{Knockout data for three regulators: ATM, RelA and p53 in two different cell populations: healthy (referred to as normal) cells and cells treated with neocarzinostatin (referred to as damaged cells). 
	}
\description{
  Knockout expression data for the kinase ATM and transcription factors p53 and RelA (Elkon et al., 2005), as well as knowledge about their targets and mutual regulatory relations in two different cell populations: healthy and damaged cells. 
}
\usage{data(damage)}
\format{
data.healthy (data.frame): 8463x3, beliefs.healthy (list): 2, model.healthy (matrix): 3x3,
data.damage (data.frame): 8463x3, beliefs.damage (list): 1, 

model.damage (matrix): 3x3.
}
\details{
	The \code{data.healthy} dataset contains log gene expression ratios  for the regulator knockouts versus control. 
	For the genes that are known to be targeted by RelA and genes that are targeted by p53 in normal conditions, \code{beliefs.healthy}  contains certainties (beliefs) that those targets are differentially expressed upon their regulator's knockdown.  
	The \code{model.healthy} matrix represents mutual signaling relations between the regulators in the healthy cells (here the model reflects that no regulator influences others). 
  The \code{data.damage} dataset contains log gene expression ratios for the regulator knockouts upon treatment with neocarzinostatin versus treatment with neocarzinostatin alone.
   For the genes that are known to be targeted by p53 in the damaged cells, \code{beliefs.damage} contains certainties (beliefs) that they are differentially expressed upon the knockdown of p53.
    The \code{model.damage} matrix  represents mutual signaling relations between the regulators in the damaged cells (here the model reflects ATM signaling down to RelA and p53).	
 }
\references{
	Elkon R, Rashi-Elkeles S, Lerenthal Y, Linhart C, Tenne T, Amariglio N, Rechavi G, Shamir R, Shiloh Y. Dissection of a DNA-damage-induced transcriptional network using a combination of microarrays, RNA interference and computational promoter analysis. Genome Biol. 2005;6(5):R43.
}
\author{
	Ewa Szczurek
}
\examples{
data(damage)
str(data.damage)
str(beliefs.damage)
print(model.damage)
}
\keyword{datasets}
