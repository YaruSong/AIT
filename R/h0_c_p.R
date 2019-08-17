
#' @export

h0_c_p<-function(n,seq_depth,Alpha,y,z,lambda,permutation=1000){

  x<- rdirmn(size=rep(seq_depth,n), alpha=Alpha)
  y=relabel(y)
  y.u=unique(y)
  null_dirm=foreach(i=1:permutation, .combine=c) %do%{
    new.y=sample(y)
    new.x=x[order(new.y),]
    new.z=matrix(z[order(new.y),])
    new.y=new.y[order(new.y)]
    tmp_d <- try(AIT_DMN_a(new.x,new.y,new.z,lambda), silent=TRUE)
    if ('try-error' %in% class(tmp_d)) {tmp_d=NA}
    else tmp_d
    #print(c(i,tmp_d))
  }
  print(length(na.omit(null_dirm)))
  print(null_dirm)
  pvalue_dirm <- length(which(null_dirm >0)) / permutation
  return(pvalue_dirm)
  
}

