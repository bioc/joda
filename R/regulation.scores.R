regulation.scores <- function(probs, model, verbose=FALSE){
  ns=colnames(probs)
  sanity.check.model(model, ns)
  if(verbose) cat("joda: the model is correctly defined\n")

  getInfluenceStrength <- function(r, value=1){
    wh=rownames(model)[model[,r]==value]
    p=probs[,colnames(probs)%in%wh]
    
    if (length(wh)>1){
      p=rowMeans(p)
    }
    names(p)=rownames(probs)
    p
  }

  
################# create a transitive closure of the model
  model.NEL <- as(model, "graphNEL")
  model.NEL.trans <- transitive.closure(model.NEL)
  model <- as(model.NEL.trans, "matrix")
  
  if(verbose) cat("joda: getting the regulation scores...\n")
  as.matrix(sapply(ns,  getInfluenceStrength))
}

sanity.check.model <- function(model, ns){
  ## check that the model has regulators as names
  if (nrow(model)!=length(ns)) stop("\nerror: the model graph has too few rows (perturbation experiments)")
  if (ncol(model)!=length(ns)) stop("\nerror: the model graph has too few columns (regulators)")
  
  if (!all(rownames(model)%in%ns))  stop("\nerror: the model graph must have the row names equal to the names of the regulators (same as the columns of probs)")
  if (!all(colnames(model)%in%ns))  stop("\nerror: the model graph must have the column names equal to the names of the regulators (same as the columns of probs)")
}
