#' An adaptive independence test for Dirichlet-multinomial modeled microbiome data after covariates adjustment
#'
#' An adaptive independence test based on the Dirichlet-multinomial distribution. It adjusts covariate,
#' then tests the independence between a vector of bacterial counts and a continuous variable or a
#' categorical variable with many levels.
#'
#'
#' @param x a vector of bacterial counts follows dirichlet-multinomial distribution.
#' @param y a continuous or categorical outcome vector. We use MLE to estimate coefficent, thus our method prefers
#' a pre-sliced (categorical) outcome so that we have enough samples to estimate parameters.
#' @param z the covariates.
#' @param lambda the panalty for the number of slices.
#' @param permutation permutation time for the test (default 1000). It is recommended to permutate
#'  100 times to reduce the amount of calculation while \code{y} is a continuous variable.
#' @param slice a bool variable to determine whether return slicing scheme (default FALSE).
#' @return results of the test: \code{statistics} and \code{p_value} (if \code{slice=TRUE}, also return \code{slicing}).
#' @seealso \code{\link{sim_x_c}}
#'
#' @examples
#'
#'##################################### a simulated example ###################################
#' rm(list=ls())
#' set.seed(2019)
#' p <- 6  # number of taxa
#' d <- 1  # number of covariates
#' samples <- 100
#' seq.depth=500
#' mean.proportion=rnorm(p,1,1)
#' groups.samples<-c(10,20,30,40)
#' gamma=t(matrix(sort(rep(0:3,p)),p,4))
#' coefficent=matrix(rnorm(p),1,p)
#' covariates=matrix(rnorm(samples),samples,1)
#' x.dc=sim_x_c(seq.depth,mean.proportion,gamma,groups.samples,coefficent,covariates)
#' y=sort(rep(1:10,10))/2
#' y.dc=y
#' # lambda_da=lambda_selection("DMN",x.dc,y.dc,permutation=100,covariates=covariates) #2.75
#' lambda_da=2.75
#' AIT_DMN_adj_p(x.dc,y.dc,covariates,lambda_da,permutation=100,slice=TRUE)
#'
#' @export


AIT_DMN_adj_p<-function(x,y,z,lambda,permutation=1000,slice=FALSE){
  if(dim(x)[1]!=length(y)){stop("The length of y does not match with the number of rows for x")}
  x=x[order(y),]
  z=matrix(z[order(y),])
  y=y[order(y)]
  y.o=y
  y=relabel(y)
  y.u=unique(y)
  t=AIT_DMN_a(x,y,z,lambda,slice=TRUE)
  dsres.d=round(t$dsval,5)
  slice_scheme=t$slices
  colnames(slice_scheme)[1:(dim(slice_scheme)[2]-1)]=unique(y.o)
  if(dsres.d==0){pvalue_dirm=1}
  else{
      null_dirm=foreach(i=1:(permutation-1), .combine=c) %do%{
          new.y=sample(y)
          new.x=x[order(new.y),]
          new.z=matrix(z[order(new.y),])
          new.y=new.y[order(new.y)]
          tmp_d <- try(AIT_DMN_a(new.x,new.y,new.z,lambda), silent=TRUE)
          if ('try-error' %in% class(tmp_d)) {tmp_d=NA}
          else tmp_d
          #print(c(i,tmp_d))
      }
      null_dirm[i+1]=dsres.d
      pvalue_dirm <- length(which(null_dirm >=dsres.d)) / permutation
  }
  if(slice){
    results=list(dsres.d,pvalue_dirm,slice_scheme)
    names(results)<-c("statistics","p_value","slicing")
    return(results)
  }else{
    results=list(dsres.d,pvalue_dirm)
    names(results)<-c("statistics","p_value")
    return(results)
  }
}
