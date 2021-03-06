#'
#' Simulate microbiome composition data
#'
#' Simulate multinomial or Dirichlet-multinomial distributed microbiome data.
#'

#' @param method two option: 'MN' or 'DMN' for multinomial or Dirichlet-multinomial model, seperately.
#' @param groups_samples a vector with number of samples in each group.
#' @param seq_depth sequencing depth for all samples. It is consistent for all samples.
#' @param mean_group a matrix consisted with taxa propotion \eqn{\phi} in each row for each group.
#' @param dispersion parameter overdispersion for Dirichlet-multinomial model. if \code{method} is 'MN',
#' then \code{dispersion} does not have effect on \code{sim_x}.
#' @return a microbiome composition matrix.
#' @examples
#' rm(list=ls())
#' set.seed(2019)
#' samples=100
#' groups.samples<-c(10,20,30,40)
#' y=sort(rep(1:10,10))/2
#' dispersion=0.2
#' seq.depth=1000
#' mean.0<-c(rep(0,6),rep(0.005,5),rep(0.01,5),rep(0.02,4),
#' rep(0.03,3),rep(0.04,2),c(0.08,0.1,0.11,0.165,0.22))
#' mean.change<-function(pi,location,alpha){
#'   mean.tmp=sort(pi)
#'   mean.tmp[0+location]= mean.tmp[0+location]+alpha*mean.tmp[20+location]
#'   mean.tmp[20+location]= (1-alpha)*mean.tmp[20+location]
#'   return(mean.tmp)
#' }
#' beta=0.6  #beta acts as the effect size
#' mean.control=mean.0
#' mean.case1=mean.change(mean.0,c(8,9,10),beta)
#' mean.case2=mean.change(mean.0,c(5,6,7),beta)
#' mean.case3=mean.change(mean.0,c(2,3,4),beta)
#' mean.group<-as.matrix(rbind(mean.control,mean.case1,mean.case2,mean.case3))
#' x.d=sim_x("DMN",groups.samples,seq.depth,mean.group,dispersion)
#' y.d=y
#' x.m=sim_x("MN",groups.samples,seq.depth,mean.group)
#' y.m=y
#'
#' @export
#'
#'
sim_x<-function(method,groups_samples,seq_depth,mean_group,dispersion=NULL)
{
  num=length(groups_samples)
  sim<-c()
  groups<-c()
  otu.cate=ncol(mean_group)
  if(method=="DMN"){
    for(i in 1:num){
      sim.tmp=simPop(groups_samples[i], otu.cate, seq_depth, mean_group[i,], dispersion)
      sim<-rbind(sim,sim.tmp$data)
    }
    return((sim))
  }
  if(method=="MN"){
    for(i in 1:num){
      sim.tmp=t(rmultinom(groups_samples[i], size = seq_depth, prob =mean_group[i,]))
      sim<-rbind(sim,sim.tmp)
    }
    return((sim))
  }
}
