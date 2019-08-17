#'
#' Select a proper lambda for adaptive independence test.
#'
#' Select a proper lambda to control type I error for the independence test.
#'
#' @param method two option: 'MN' or 'DMN' for multinomial and Dirichlet-multinomial model seperately (default 'DMN').
#' @param x a matrix represents microbiomal composition.
#' @param y a continuous or categorical covariate vector.
#' @param permutation permutation time for the test (default 1000). It is recommended to permutate 100 times to reduce
#' the amount of calculation when \code{y} is continuous.
#' @param sig significance level, that's type I error.
#' @param covariates covariates for microbiome composition.
#' @return a proper penalty lambda for the the adaptive independence test.
#' @seealso \code{\link{x.m}},\code{\link{x.d}}
#' @examples
#'
#' data(x.m)
#' data(x.d)
#' y=sort(rep(1:10,10))/2
#' y.d=y.m=y
#' lambda_m=lambda_selection("MN",x.m,y.m,100)
#' lambda_d=lambda_selection("DMN",x.d,y.d,100)
#' @export
#'

lambda_selection<-function(method="DMN",x,y,permutation=100,sig=0.05,covariates=NULL){
  y=relabel(y)
  x=x[order(y),]
  n=nrow(x)

  if(length(unique(rowSums(x)))==1){
    seq_depth=rowSums(x)[1]
  }else{
    stop("please Rarefy your OTU table")
  }

  permu=permutation

  # if(length(z)>0)
  # h0_c_p(n,seq_depth,Alpha,y,z,end,permu)  # end
  #
  if(is.null(covariates)){
    y=y[order(y)]
    start=0.5
    end=(ncol(x)-1)/2
    pai=colSums(x)/sum(x)
    theta=theta_calculate(x)

    if(method=="MN"){
      theta=0
    }
    p.end=h0_p(method,n,seq_depth,pai,theta,y,end,permu)
    if(p.end>0.05){
      return(end)
    }
    p.start=h0_p(method,n,seq_depth,pai,theta,y,start,permu)
    if(p.start<=0.05){
      return(start)
    }

    p.mid=1

    while(abs(p.mid-sig)>=0.015){
      mid=(start+end)/2
      p.mid=h0_p(method,n,seq_depth,pai,theta,y,mid,permu)
      if (p.mid>=0.065){start=mid}
      else if(p.mid<=0.035){end=mid}
      else if(p.start<=0.05 & p.end<=0.05 ){return(p.start)}
      else{break}
    }
    return(mid)
  }else{
    if(method=="DMN"){
      z=matrix(covariates[order(y),])
      y=y[order(y)]
      start=0.5
      end=(ncol(x)-1)
      B_e=MGLMreg(x~z, dist="DM")@coefficients
      Alpha <- exp(cbind(rep(1,n),z)%*%B_e)

      p.end=h0_c_p(n,seq_depth,Alpha,y,z,end,permu)
      if(p.end>0.05){
        return(end)
      }

      p.start=h0_c_p(n,seq_depth,Alpha,y,z,start,permu)
      if(p.start<=0.05){
        return(start)
      }


      p.mid=1

      while(abs(p.mid-sig)>=0.015){
        mid=(start+end)/2
        p.mid=h0_c_p(n,seq_depth,Alpha,y,z,mid,permu)
        if (p.mid>=0.065){start=mid}
        else if(p.mid<=0.035){end=mid}
        else if(p.start<=0.05 & p.end<=0.05 ){return(p.start)}
        else{break}
      }
      return(mid)
    }else{
      stop("our method can't adjustment covariates for models except for multinomial distribution")
    }

  }

}
