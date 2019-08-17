#' @export
#' @importFrom stats rmultinom
#' @import Rcpp
#' @importFrom dirmult simPop
#' @import doParallel
#' @import foreach
#' @useDynLib AIT
h0_p<-function(method="DMN",n,seq_depth,pai,theta,y,lambda,permutation=1000){

  if(method=="DMN"){
    sim.data=simPop(n, ncol(pai), seq_depth, pai, theta)$data
  }else if(method=="MN"){
    sim.data=t(rmultinom(n, size = seq_depth, prob =pai))
  }else{
    stop("overdispersion theta should not be samller than 0.")
  }
  x=sim.data
  y=relabel(y)
  y.u=unique(y)
  if(method=="DMN"){
    dsres.d=AIT_DMN(x,y,y.u,lambda)
    null_dirm=foreach(i=1:permutation, .combine=c) %do%{
      y=sample(y)
      x=x[order(y),]
      y=y[order(y)]
      AIT_DMN(x,y,y.u,lambda)
    }
    pvalue_dirm <- length(which(null_dirm >0)) / permutation
    return(pvalue_dirm)
  }else if(method=="MN"){
    dsres.m=AIT_MN(x,y,y.u,lambda)
    null_mult=foreach(i=1:permutation, .combine=c) %do%{
      y=sample(y)
      x=x[order(y),]
      y=y[order(y)]
      AIT_MN(x,y,y.u,lambda)
    }
    pvalue_mult <- length(which(null_mult >0)) / permutation
    return(pvalue_mult)
  }else{
    stop("'method' only has two options: 'MN' and 'DMN', please choose one from them")
  }
}

