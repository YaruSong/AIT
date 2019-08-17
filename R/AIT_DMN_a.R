#' @export
#' @importFrom MGLM MGLMreg
#' @importFrom GUniFrac Rarefy
AIT_DMN_a<-function(x,y,z,lambda,slice=FALSE){
  yu=unique(y)
  z=cbind(rep(1,nrow(x)),z)
  myfunction=regression
  leny=length(y)
  lenyu=length(yu)
  lpd = -lambda * log(leny)
  alpha=myfunction(x,z,1,nrow(x))
  ctab<-array(rep(0,(lenyu+1)*lenyu), c(lenyu+1,lenyu))
  for(i in y){
      ctab[i+1,i]=ctab[i+1,i]+1
  }
  for( i in 1:length(unique(y))){
      ctab[i+1,]=ctab[i+1,]+ctab[i,]
  }
  ctab1=rowSums(ctab)
  clumpnum=length(unique(y))+1

  score=rep(0,clumpnum)
  idx=rep(-1,clumpnum)

  # print(colMeans(alpha))

  for(i in 2:clumpnum){
      all.cutsc<-c()
      all.idx<-c()
      for(j in 1:(i-1)){
          cutsc = score[j]+lpd
          group<-c((ctab1[j]+1):ctab1[i])
          # print(group)

          alpha_tmp=myfunction(x,z,group[1],group[length(group)])
           # print(colMeans(alpha_tmp))
          for(k in group){
              for(m in 1:ncol(x)){
                  tmp_a=x[k,m]+alpha_tmp[k-group[1]+1,m]
                  if((tmp_a)>0){
                      cutsc=cutsc+lgamma(tmp_a)
                  }
                  tmp_b=alpha_tmp[k-group[1]+1,m]
                  if((tmp_b)>0){
                      cutsc=cutsc-lgamma(tmp_b)
                  }
              }
            tmp_c=sum(alpha_tmp[k-group[1]+1,])
            tmp_d=sum(alpha_tmp[k-group[1]+1,]+x[k,])
            if(tmp_c>0){
              cutsc=cutsc+lgamma(tmp_c)
            }
            if(tmp_d>0){
              cutsc=cutsc-lgamma(tmp_d)
            }

          }
          # print("cutsc")
          # print(all.cutsc)
          all.cutsc<-c(all.cutsc,cutsc)
          all.idx<-c(all.idx,j)
          # print(all.cutsc)
      }


      score[i] = max(all.cutsc)
      idx[i] = all.idx[which(all.cutsc==score[i])]
      # print("score")
      # print(score)
      # print(all.cutsc)
      # print(c(max(all.cutsc),idx[i]))
  }

  sum2= lpd

  for(i in 1:nrow(x)){
      for(j in 1:ncol(x)){
          tmp_a=x[i,j]+alpha[i,j]
          if(tmp_a>0){
              sum2=sum2+lgamma(tmp_a)
          }
          tmp_b=alpha[i,j]
          if(tmp_b>0){
              sum2=sum2-lgamma(tmp_b)
          }

      }
    tmp_c=sum(alpha[i,])
    tmp_d=sum(alpha[i,]+x[i,])
    if(tmp_c>0){
      sum2=sum2+lgamma(tmp_c)
    }
    if(tmp_d>0){
      sum2=sum2-lgamma(tmp_d)
    }
  }

  mlik = round(score[clumpnum]-sum2,5)
  mlik
  if(slice){
    slicenum = 0;
    flag = clumpnum
    # print(clumpnum)
    while(flag > 0){
      flag = idx[flag]
      #print(idx[flag])
      slicenum= slicenum + 1
    }

    slicenum=slicenum-1
    slices<-as.matrix(array(rep(0,slicenum*(lenyu+1)), c(slicenum,lenyu+1)))

    s = "s"
    rslices<-c()
    for(i in 1:slicenum){
      rslices<-c(rslices,paste(s,as.character(i),sep=""))
    }
    rownames(slices)=rslices

    cslices<-c()
    for(j in 1:lenyu){
      cslices<-c(cslices,as.character(j))
    }
    cslices<-c(cslices, "total")
    colnames(slices)=cslices

    spos=rep(0,slicenum+1)
    flag = clumpnum
    for(i in (slicenum+1):1){
      spos[i] = flag;
      flag = idx[flag];
    }

    spos[1] = 1
    for(i in 1:slicenum){
      slices[i, 1:lenyu] = ctab[spos[i+1], ] - ctab[spos[i], ];
      slices[i, lenyu+1] = sum(slices[i, 1:lenyu]);
    }
    return(list(dsval=mlik,slices=slices))
    #print(mlik)
    # print(slices)
  }
  else{
    return(mlik)
  }

}



#' @export
regression<-function(x,z,s,e){
  colnames(x)=NULL
  row.names(x)=NULL

  x_input=x[s:e,]
  z_input=z[s:e,]
  alpha_out<-matrix(0.001,nrow(x_input),ncol(x_input))
  col_zero=which(colSums(x_input!=0)<=2)
  if(length(col_zero)>=1){
    x_nonzero=as.matrix(x_input[,- as.numeric(col_zero)])
  }else{x_nonzero=x_input}


  #if(length(unique(z_input[,2]))==1){z_input=z_input[,-2]}
  B_nozero<-matrix(0,ncol(z_input),ncol(x_nonzero))

  B_e <- try(MGLMreg(x_nonzero~z_input-1, dist="DM")@coefficients, silent=TRUE)
  if ('try-error' %in% class(B_e)) {
    B_e=NA
    return(alpha_out)}

  meaning=as.numeric(substr(colnames(B_e),5,nchar(colnames(B_e))))
  B_nozero[,meaning]=B_e
  alpha_e=exp(z_input%*%B_nozero)
  alpha_e[,-meaning]=0
  if(length(col_zero)>=1){
    alpha_out[,-col_zero]=alpha_e
  }else{alpha_out=alpha_e}

  return(alpha_out)
}


