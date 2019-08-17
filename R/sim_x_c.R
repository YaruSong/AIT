#'
#' Simulate microbiome composition data with covariates
#'
#' Simulate Dirichlet-multinomial distributed microbiome data subjected to covariates.
#'
#'
#' @param seq_depth sequencing depth for all samples. It is consistent for all samples.
#' @param mean_proportion a vector consisted with taxa propotion parameter \eqn{\mu} for each taxa.
#' @param gamma coefficient matrix \eqn{\gamma}, each row corresponds to each group, and each column is for taxa.
#' @param groups_samples a vector with number of samples in each group.
#' @param coefficent coefficient \eqn{c} of \code{covariates} for each taxon .
#' @param covariates covariates matrix \eqn{z} for microbiome composition.
#' @return a microbiome composition matrix.
#' @examples
#'
#' rm(list=ls())
#' set.seed(2019)
#' p <- 6  # number of taxa
#' d <- 1  # number of covariates
#' samples <- 100
#' seq.depth=500
#' mean.proportion=rnorm(p,1,1)
#' groups.samples<-c(10,20,30,40)
#' gamma=t(matrix(sort(rep(0:3,p)),p,4))
#' #gamma=t(matrix(0,p,4))
#' coefficent=matrix(rnorm(p),1,p)
#' covariates=matrix(rnorm(samples),samples,1)
#' y=sort(rep(1:10,10))
#' x.dc=sim_x_c(seq.depth,mean.proportion,gamma,groups.samples,coefficent,covariates)
#'
#'
#' @export
sim_x_c<-function(seq_depth,mean_proportion,gamma,groups_samples,coefficent,covariates)
{
  samples=sum(groups_samples)
  groups=nrow(gamma)
  p=length(mean_proportion)
  d=ncol(covariates)
  if(length(groups_samples)!=groups){
    stop("the length of gamma must be corresponded to samples_group information")
  }

  Z=matrix(0,nrow=samples,ncol=groups)
  tmp=0
  for(k in 1:groups){
    Z[(tmp+1):(tmp+groups_samples[k]),k]=1
    tmp=tmp+groups_samples[k]
  }
  Z=cbind(rep(1,samples),Z,covariates)

  B=matrix(0,nrow=d+1+groups,ncol=p)
  B[1,]=mean_proportion
  B[2:(groups+1),]=gamma
  B[(groups+2):nrow(B),]=coefficent

  Alpha <- exp(Z%*%B)
  x<- rdirmn(size=rep(seq_depth,samples), alpha=Alpha)
}

