#'
#' Calculate overdispersion for microbiomal data
#'
#' Calculate overdispersion \eqn{\theta} for microbiomal data. If \eqn{\theta>0}, we can assume microbiome follows
#' Dirichlet-multinomial distribution; otherwise, it follows multinomial distribution.
#'
#'
#' @param x a matrix follows multinomial or Dirichlet_multinomial distribution.
#' @return overall moment estimator for overdispersion \eqn{\theta}.
#' @seealso \code{\link{x.d}}, \code{\link{x.m}}
#' @examples
#' data(x.d)
#' data(x.m)
#' theta_calculate(x.d)
#' theta_calculate(x.m)
#' @export
theta_calculate<-function(x){
  phi<-c()
  Ni<-c()
  for(i in 1:ncol(x)){
    phi<-c(phi,sum(x[,i])/sum(x))
  }
  Ni=rowSums(x)
  N=sum(x)
  n=nrow(x)
  Nc=(N-sum(Ni^2)/N)/(n-1)
  Sj<-c()
  x.col=ncol(x)
  x.row=nrow(x)
  for(j in 1:x.col){
    tmp=0
    for(i in 1:x.row){
      tmp=tmp+Ni[i]*((x[i,j]/Ni[i]-phi[j])^2)
    }
    tmp=tmp/(n-1)
    Sj<-c(Sj,tmp)
  }
  Gj<-c()
  for(j in 1:x.col){
    tmp=0
    for(i in 1:x.row){
      tmp=tmp+Ni[i]*(x[i,j]/Ni[i])*(1-x[i,j]/Ni[i])
    }
    tmp=tmp/(N-n)
    Gj<-c(Gj,tmp)
  }
  theta=sum(Sj-Gj)/sum(Sj+(Nc-1)*Gj)
  return(theta)
}





