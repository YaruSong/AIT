#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List AIT_MN_slice(Rcpp::NumericMatrix x,Rcpp::NumericVector y,Rcpp::NumericVector yu, double lambda)
{
  int leny=y.size();
  int lenyu=yu.size();
  int xrow=x.nrow();
  int xcol=x.ncol();
  double lpd = -lambda * log((double)leny);
  Rcpp::NumericMatrix ctab(leny+1, lenyu);

  for(int i=0;i<leny+1;i++){
    for(int j=0;j<lenyu;j++){
      ctab(i,j)=0;
    }
  }
  for(int i=0;i<leny;i++){
    int j=y[i];
    ctab(j,j-1)=ctab(j,j-1)+1;
  }
  for(int i=1;i<lenyu+1;i++){
    for(int j=0;j<lenyu;j++){
      ctab(i,j)=ctab(i,j)+ctab(i-1,j);
    }
  }

  int lenyp=leny+1;
  Rcpp::NumericVector ctab1(lenyp);
  for(int i=0;i<leny+1;i++){
    int k=0;
    for(int j=0;j<lenyu;j++){
      k=k+ctab(i,j);
    }
    ctab1[i]=k;
  }
  int clumpnum=lenyu+1;
  Rcpp::NumericVector score(clumpnum);
  Rcpp::IntegerVector idx(clumpnum);
  for(int i=0;i<clumpnum;i++){
    score[i]=0.0;
    idx[i]=-1;
  }
  for(int i=1;i<clumpnum;i++){
    Rcpp::NumericVector allcutsc(i);
    Rcpp::IntegerVector allidx(i);
    for(int k=0;k<i;k++){
      allcutsc[k]=-100000000.0;
      allidx[k]=-1;
    }
    for(int j=0;j<=i-1;j++){
      double cutsc=score[j]+lpd;
      int l=ctab1[i]-ctab1[j]-1;
      int group[l];

      for(int m=0;m<l+1;m++){
        group[m]=ctab1[j]+m;
      }
      double nh=0.0;
      for(int k=group[0];k<=group[l];k++){
        for(int m=0;m<xcol;m++){
          nh=nh+x(k,m);
        }
      }

      Rcpp::NumericVector n_hj(xcol);
      for(int m=0;m<xcol;m++){
        n_hj[m]=0;
        for(int k=group[0];k<=group[l];k++){
          n_hj[m]=n_hj[m]+x(k,m);
        }
      }
      for(int m=0;m<xcol;m++){
        if(n_hj[m]>0){
          cutsc=cutsc+(double)n_hj[m]*log((double)n_hj[m]/nh);
        }
      }
      allcutsc[j]=cutsc;
      allidx[j]=j;
    }
    double score_max=allcutsc[0];
    int idx_max=-1;
    for(int l=0;l<i;l++){
      if(allcutsc[l]>score_max){
        score_max=allcutsc[l];
        idx_max=l;
      }
    }
    score[i]=score_max;
    idx[i]=idx_max;
  }

  Rcpp::NumericVector nj(xcol);
  for(int s=0;s<xcol;s++){
    nj[s]=0;
    for(int l=0;l<xrow;l++){
      nj[s]=nj[s]+x(l,s);
    }
  }
  double xsum=0.0;
  for (int i=0;i<xrow;i++)
  {
    for (int j=0;j<xcol;j++)
    {
      xsum=xsum+x(i,j);
    }
  }
  double sum2=lpd;;
  for(int j=0;j<xcol;j++){
    if(nj[j]>0)
      sum2=sum2+(double)nj[j]*log(double(nj[j]/xsum));
  }
  double max_likelihood = 0.0;
  max_likelihood=score[clumpnum-1]-sum2;


  int slice_num=0;
  int flag=clumpnum-1;

  while(flag>0){
    flag=idx[flag];
    slice_num++;
  }

  Rcpp::IntegerMatrix slices(slice_num, lenyu+1);
  Rcpp::CharacterVector rslice(slice_num);
  Rcpp::CharacterVector cslice(lenyu+1);
  std::string s="s";
  std::string instr;
  for(int i=0;i<slice_num;i++){
    std::stringstream ss;
    std::string str;
    ss<<i+1;
    ss>>str;
    instr=s+str;
    rslice[i]=instr;
  }
  for(int i=0;i<lenyu;i++){
    std::stringstream ss;
    std::string str;
    ss<<i+1;
    ss>>str;
    cslice[i]=str;
  }
  cslice[lenyu]="total";

  Rcpp::List dimnames = Rcpp::List::create(rslice, cslice);
  slices.attr("dimnames") = dimnames;

  int spos[slice_num];
  flag=clumpnum-1;
  for(int i=slice_num;i>=0;i--){
    spos[i]=flag;
    flag=idx[flag];

  }
  spos[0]=0;

  for(int i=0;i<slice_num;i++){
    int k=0;
    for(int j=0;j<lenyu;j++){
      slices(i,j)=ctab(spos[i+1],j)-ctab(spos[i],j);
      k=k+slices(i,j);
    }
    slices(i,lenyu)=k;
  }
  Rcpp::List slicing_res = Rcpp::List::create(Rcpp::Named("dsval")=max_likelihood, Rcpp::Named("slices")=slices);
  return slicing_res;

  }

