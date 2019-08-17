#' An adaptive independence test for multinomial distributed microbiome data
#'
#' An adaptive independence test based on the multinomial distribution. It tests the independence
#' between a vector of bacterial counts and a continuous variable or a categorical variable with
#' many levels.
#'
#' @param x a vector of bacterial counts follows multinomial distribution.
#' @param y a continuous or categorical outcome vector.
#' @param lambda the penalty for the number of slices.
#' @param permutation permutation time for the test (default 1000). It is recommended to permutate
#'  100 times to reduce the amount of calculation while \code{y} is a continuous variable.
#' @param slice a bool variable to determine whether return slicing scheme (default FALSE).
#' @return results of the test: \code{statistics} and \code{p_value} (if \code{slice=TRUE}, also return \code{slicing}).
#' @examples
#' set.seed(2019)
#' data(x.m)
#' y=sort(rep(1:10,10))/2
#' y.m=y
#' theta=theta_calculate(x.m)
#' lambda_m=lambda_selection("MN",x.m,y.m,1000)
#' AIT_MN_p(x.m,y.m,lambda_m,1000,slice = TRUE)
#' lambda_d=lambda_selection("DMN",x.m,y.m,1000)
#' AIT_DMN_p(x.m,y.m,lambda_d,1000,slice = TRUE)
#'
#' y.m=sort(runif(100, min=0.01, max=1))
#' lambda_m=lambda_selection("MN",x.m,y.m,100)
#' AIT_MN_p(x.m,y.m,lambda_m,1000,slice = TRUE)
#' @export
AIT_MN_p<-function(x,y,lambda,permutation=1000,slice=FALSE){
  if(dim(x)[1]!=length(y)){stop("The length of y does not match with the number of rows for x")}
  y.o=y
  x=x[order(y),]
  y=y[order(y)]
  y=relabel(y)
  y.u=unique(y)
  dsres.m=AIT_MN(x,y,y.u,lambda)
  slice_scheme=AIT_MN_slice(x,y,y.u,lambda)$slice
  colnames(slice_scheme)[1:(dim(slice_scheme)[2]-1)]=unique(y.o)
  if(dsres.m==0){pvalue_mult=1}
  else{
    i=1
    null_mult=foreach(i=1:(permutation-1), .combine=c) %do%{
      y=sample(y)
      x=x[order(y),]
      y=y[order(y)]
      AIT_MN(x,y,y.u,lambda)
    }
    null_mult[i+1]=dsres.m
    pvalue_mult <- length(which(null_mult >=dsres.m)) / permutation
  }
  if(slice){
    results=list(dsres.m,pvalue_mult,slice_scheme)
    names(results)<-c("statistics","p_value","slicing")
    return(results)
  }else{
    results=list(dsres.m,pvalue_mult)
    names(results)<-c("statistics","p_value")
    return(results)
  }
}
