#' An adaptive independence test for Dirichlet-multinomial modeled microbiome data
#'
#' An adaptive independence test based on the Dirichlet-multinomial distribution. It tests
#' the independence between a vector of bacterial counts and a continuous variable or a
#' categorical variable with many levels.
#'
#' @param x a vector of bacterial counts follows dirichlet-multinomial distribution.
#' @param y a continuous or categorical outcome vector.
#' @param lambda the panalty for the number of slices.
#' @param permutation permutation time for the test (default 1000). It is recommended to permutate
#'  100 times to reduce the amount of calculation while \code{y} is a continuous variable.
#' @param slice a bool variable to determine whether return slicing scheme (default FALSE).
#' @return results of the test: \code{statistics} and \code{p_value} (if \code{slice=TRUE}, also return \code{slicing}).
#' @seealso \code{\link{x.d}},\code{\link{individual_index}},\code{\link{taxa_family}}
#' @examples
#'
#'##################################### a simulated example ###################################
#' library(GUniFrac)
#' set.seed(2019)
#' data(x.d)
#' y=sort(rep(1:10,10))/2
#' y.d=y
#' lambda_d=lambda_selection("DMN",x.d,y.d,100)
#' AIT_DMN_p(x.d,y.d,lambda_d,1000,slice = TRUE)
#'
#' set.seed(2019)
#' y.d=sort(runif(100, min=0.01, max=1))
#' #lambda_d=lambda_selection("DMN",x.d,y.d,100)  #3.125
#' lambda_d=3.125
#' AIT_DMN_p(x.d,y.d,lambda_d,100,slice = TRUE)
#'
#' ####################################### a real example #####################################
#'
#' set.seed(2019)
#' data(individual_index)
#' data(taxa_family)
#' index.tmp=subset(individual_index,individual_index$Age<=40&individual_index$Country=="Malawi")
#' age=index.tmp$Age
#' taxa_age=taxa_family[row.names(index.tmp),]
#' taxa_age_sub <- taxa_age[,(colSums(taxa_age != 0) ) >nrow(taxa_age)*0.25]
#' taxa_age_sub=taxa_age_sub[,apply(taxa_age_sub, 2, function(x){max(x)})>sum(taxa_age_sub)*0.0002]
#' dim(taxa_age_sub)
#' x=Rarefy(taxa_age_sub,depth=10000)$otu.tab.rff   #rarefy to the same sequencing depth
#' x=x[,colSums(x)>0]
#' age.group=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3)
#' y=age
#' for(i in 1:length(age)){
#'  if(age[i]>3){
#'    y[i]=ceiling(age[i])}
#'  for(j in 1:(length(age.group)-1)){
#'    if(age.group[j]<age[i] &age[i] <=age.group[j+1]){
#'      y[i]=age.group[j+1]
#'    }
#'  }
#'}
#'
#' theta_calculate(x)  # here theta>0.1, thus we use function AIT_DMN_p
#' lambda_d=lambda_selection(method="DMN",x,y,100)  #
#' test.results=AIT_DMN_p(x,y,lambda_d,1000,slice=TRUE)
#' print(test.results)
#'
#' @export
#' @importFrom GUniFrac Rarefy
#' @importFrom MGLM rdirmn
AIT_DMN_p <- function(x,y,lambda,permutation=1000,slice=FALSE){
  if(dim(x)[1]!=length(y)){
    stop("The length of y does not match with the number of rows for x")}
  x=x[order(y),]
  y=y[order(y)]
  y.o=y
  y=relabel(y)
  y.u=unique(y)
  dsres.d=AIT_DMN(x,y,y.u,lambda)
  slice_scheme=AIT_DMN_slice(x,y,y.u,lambda)$slice
  colnames(slice_scheme)[1:(dim(slice_scheme)[2]-1)]=unique(y.o)
  if(dsres.d==0){pvalue_dirm=1}
  else{
    i=1
    null_dirm=foreach(i=1:(permutation-1), .combine=c) %do%{
      y=sample(y)
      x=x[order(y),]
      y=y[order(y)]
      AIT_DMN(x,y,y.u,lambda)
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
