differential.probs <- function(data, beliefs, verbose=FALSE, plot.it=FALSE){
  sanity.check.beliefs(data,beliefs)  
  if (verbose) cat("joda: Input correctly defined \n")
  
  if(class(data)=="data.frame")
    data=as.matrix(data)
  ##Apply bgmm to infer probabilities of differential expression
  infer.probabilities(data, beliefs, verbose, plot.it)
} 

####################
## Data is a matrix with rows for all genes and columns for regulators. Each entry is the log expression ratio perturbation vs control.
## Beliefs should be a list with the names for those regulators for which we have some genes known to respond in some way to their perturbation.
## Each entry of the list is a matrix with rows for the known genes.
## Each row is a distribution over  the differential and unchanged cluster (2 columns) or over down, up-regulated and unchanged cluster of genes  (3 columns).
## The rownames of each matrix must be a subset of rows in the data.
## For each regulator, the belief-based mixture model will have the  number of model components equal to the number of columns in the corresponding beliefs.
## If no beliefs are given, unsupervised mixture modeling will be applied.
####################
sanity.check.beliefs <- function(data, beliefs){
  if (is.null(beliefs))
    beliefs=list()

  if(!class(data)%in%c("matrix", "data.frame"))
    stop(paste("\ndifferential.probs Error: data must be a matrix."))
  
  if(!class(beliefs)=="list")
    stop(paste("\ndifferential.probs Error: beliefs must be a list."))
 
  ##check the knowns
  regs=colnames(data)	
  c1 <- all(names(beliefs)%in%regs)
  if (!c1) stop(paste("\ndifferential.probs Error: beliefs must be a list of  matrices, with names as subset of the data columns (regulators)."))
  

  ##each entry in the beliefs should be a matrix  or NULL
  ck <- function(kn){
    k=beliefs[[kn]]
    OK <- TRUE
    if (!is.null(k)){
      if (!is.matrix(k)){
        cat(" \n Error:. For a given regulator, the beliefs must be given as a matrix\n")
        return(FALSE)
      }
      gens=rownames(data)
      if (! all(rownames(k)%in%gens)){
        cat("\n Error: the rownames of the belief matrices must be a subset of the rownames of the data'\n")
        return(FALSE)
      }
      if (!ncol(k)%in%c(2,3)){
        cat("\n Error: the belief matrices must have either two or three columns\n")
        return(FALSE)
      }
      OK <- all( k <= 1)
      OK <- all( k >= 0)
      OK <- all(apply(k, 1, sum)==1)
      if (!OK)    cat("\n Error: the belief matrices must have values in (0,1), each row being a distribution over the columns.\n")
    }
    OK
  }
  
  c2 <- all(sapply(regs, ck))
  if (!(c2))  stop("")
}

### Use the bgmm package to infer the probability of differential expression under the experiment, given the known examples
infer.probabilities <- function(data, beliefs,verbose=FALSE, plot.it=FALSE){
  data=data[order(rownames(data)),,drop=FALSE]
  if(plot.it){
    nr=length(colnames(data))
    if (nr<3)
      par(mfrow=c(1,nr))
    else
      par(mfrow=c(ceiling(nr/3),3))
  }
  getProbs<-function(nam){
    b <- beliefs[[nam]]
    d <- data[,nam, drop=FALSE]
    rownames(d)=rownames(data)
    p=data[,nam,drop=FALSE]
    if (is.null(b)){
      if (verbose){
        cat(paste("\nInferring probabilities of differential expression under the knockdown of ", nam,"...") ,"\n")
        cat("Applying unsupervised mixture modeling\n")
      }
      bels <- unsupervised(X=d, k=2)      
      resp=DEprobs(bels, verbose=FALSE)
      p[,nam]=resp$diff.p.X
    }else{
      if (verbose){
        cat(paste("\nInferring probabilities of differential expression under the knockdown of ", nam,"...") ,"\n")
        cat("Applying belief-based mixture modeling\n")
      }
      known=d[rownames(d)%in%rownames(b),,drop=FALSE]
      X=d[!rownames(d)%in%rownames(b), ,drop=FALSE]
      X=as.matrix(X)
      known=as.matrix(known)
      bels <- belief(X=X,knowns=known, B=b, k=2, all.possible.permutations=TRUE)
      resp=DEprobs(bels, verbose=FALSE)

      p[rownames(X),nam]=resp$diff.p.X[rownames(X)]
      p[rownames(known),nam]=resp$diff.p.knowns[rownames(known)]
    }
    if (verbose){
      cat(paste("The parameters of the model for ",nam,":",sep=""), "\n")
      di=resp$diff.c
      un=c(1,2)[-di]
      
      m=matrix(0,nrow=3,ncol=2)
      colnames(m)=c("differential","unchanged")  
      rownames(m)=c("Mixing proportions:","Means:","Variances:")
      m[1,]=c(bels$pi[di],bels$pi[un]) 
      m[2,]=c(bels$mu[di,],bels$mu[un,])
      m[3,]=c(bels$cvar[di,,],bels$cvar[un,,])
      print(m)  
    }
    if (plot.it){
      plot(bels)
      title(nam)
      if(bels$k==2)
        legend("topright",legend=c("differential","unchanged"), lty=1, col=(c(resp$diff.c, c(1,2)[-resp$diff.c])+1),inset=0.05)
    }
    p
  }
  ps=sapply(colnames(data), getProbs)
  rownames(ps)=rownames(data)
  ps  
}

