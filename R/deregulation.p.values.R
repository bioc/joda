deregulation.p.values <- function(data.1, beliefs.1, model.1, data.2, beliefs.2, model.2, N=100, verbose=FALSE) {
  deregulation = NULL
  
  if (verbose)   cat("Initial calculation of differential probabilities\n")
 
  probs.1 = differential.probs(data.1, beliefs.1) 
  probs.2  = differential.probs(data.2, beliefs.2) 
  regulation.1 = regulation.scores(probs.1, model.1)
  regulation.2  = regulation.scores(probs.2, model.2) 
  deregulationOrg     = deregulation.scores(regulation.1, regulation.2) 
  abs.deregulationOrg = abs(deregulationOrg)

  deregulationScore = matrix(0, nrow(deregulationOrg), ncol(deregulationOrg))
  colnames(deregulationScore) = colnames(deregulationOrg)
  rownames(deregulationScore) = rownames(deregulationOrg)

  n.probs.1 = probs.1
  n.probs.2  = probs.2

  if (verbose)   cat("Permutation steps:\n")

  tmp <- replicate(N, {
      n.probs.1 = apply(probs.1, 2, sample)
      n.probs.2 = apply(probs.2, 2, sample)
      
      regulation.1  = regulation.scores(n.probs.1, model.1)
      regulation.2  = regulation.scores(n.probs.2, model.2) 
      
      # deregulation contains a matrix of deregulation scores
      deregulation      = deregulation.scores(regulation.1, regulation.2) 
      deregulationScore <<- deregulationScore + (abs.deregulationOrg < abs(deregulation))
      if (verbose)   cat("#")
      NULL
  })
  list(deregulation.p.values = (deregulationScore/N), deregulationOrg = deregulationOrg)
}

