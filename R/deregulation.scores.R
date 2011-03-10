deregulation.scores <- function(reg.scores1,reg.scores2,verbose=FALSE){
  if(! all(dim(reg.scores1)==dim( reg.scores2))) stop("\nerror: the regulation score matrices must be the same size")
  if (!all(rownames(reg.scores1)==rownames(reg.scores2))) stop("\nerror: the regulation score matrices must have the same rownames")
  if (!all(colnames(reg.scores1)==colnames(reg.scores2))) stop("\nerror: the regulation score matrices must have the same colnames")
  if (verbose){
     cat("joda: calculating deregulation scores\n")
     cat("subtracting reg.scores1 from reg.scores2\n")
}
  reg.scores2-reg.scores1
}
