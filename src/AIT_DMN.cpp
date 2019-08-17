#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
double AIT_DMN(Rcpp::NumericMatrix x,Rcpp::NumericVector y,Rcpp::NumericVector yu, double lambda)
{

    int leny=y.size();
    int lenyu=yu.size();
    int xrow=x.nrow();
    int xcol=x.ncol();
    double lpd = -lambda * log((double)leny);

    double phi[xcol];
    double alpha[xcol];
    double Ni[xrow];
    double Sj[xcol];
    double Gj[xcol];
    double N=0.0;
    double Nc=0.0;
    double Ni_2_sum=0.0;
    int n=xrow;
    double theta=0.0;
    double alpha_p=0.0;
    double head=0.0;
    double tail=0.0;

    for(int j=0;j<xcol;j++){
        double tmp=0.0;
        for(int i=0;i<xrow;i++){
            tmp=tmp+x(i,j);
        }
        phi[j]=tmp;
        N=N+tmp;
    }
    for(int j=0;j<xcol;j++){
        phi[j]=phi[j]/N;
    }
    for(int i=0;i<xrow;i++){
        double tmp=0.0;
        for(int j=0;j<xcol;j++){
            tmp=tmp+x(i,j);
        }
        Ni[i]=tmp;
        Ni_2_sum=Ni_2_sum+Ni[i]*Ni[i];
    }
    Nc=(N-Ni_2_sum/N)/(n-1);

    for(int j=0;j<xcol;j++){
        double tmp_s=0.0;
        double tmp_g=0.0;
        for(int i=0;i<xrow;i++){
            tmp_s=tmp_s+Ni[i]*(x(i,j)/Ni[i]-phi[j])*(x(i,j)/Ni[i]-phi[j]);
            tmp_g=tmp_g+Ni[i]*(x(i,j)/Ni[i])*(1-x(i,j)/Ni[i]);
        }
        tmp_s=tmp_s/(n-1);
        tmp_g=tmp_g/(N-n);
        Sj[j]=tmp_s;
        Gj[j]=tmp_g;
    }

    for(int j=0;j<xcol;j++){
        head=head+Sj[j]-Gj[j];
        tail=tail+Sj[j]+(Nc-1)*Gj[j];
    }
    theta=head/tail;
    alpha_p=(1-theta)/theta;

    for(int j=0;j<xcol;j++){
        alpha[j]=phi[j]*alpha_p;
    }


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
        ctab1[i]=0;
        for(int j=0;j<lenyu;j++){
            ctab1[i]=ctab1[i]+ctab(i,j);
        }
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
            double irows=group[l]-group[0]+1;
            double pi_tmp[xcol];
            double alpha_tmp[xcol];
            double N_tmp=0.0;
            for(int p=0;p<xcol;p++){
                double tmp=0.0;
                for(int i=group[0];i<=group[l];i++){
                    tmp=tmp+x(i,p);
                }
                pi_tmp[p]=tmp;
                N_tmp=N_tmp+tmp;
            }
            for(int p=0;p<xcol;p++){
                pi_tmp[p]=pi_tmp[p]/N_tmp;
                alpha_tmp[p]=alpha_p*pi_tmp[p];
            }
            for(int k=group[0];k<=group[l];k++){
                for(int m=0;m<xcol;m++){
                    if((x(k,m)+alpha_tmp[m])>0){
                        cutsc=cutsc+lgamma(x(k,m)+alpha_tmp[m]);
                    }
                }
            }
            for(int m=0;m<xcol;m++){
                if(alpha_tmp[m]>0){
                    cutsc=cutsc-irows*lgamma(alpha_tmp[m]);
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
    double xsum=0.0;
    for (int i=0;i<xrow;i++)
    {
        for (int j=0;j<xcol;j++)
        {
            xsum=xsum+x(i,j);
        }
    }

    double sum2=lpd;
    for(int i=0;i<xrow;i++){
        for(int j=0;j<xcol;j++){
            if((x(i,j)+alpha[j])>0){
                sum2=sum2+lgamma(x(i,j)+alpha[j]);
            }
        }
    }
    for(int j=0;j<xcol;j++){
        if(alpha[j]>0){
            sum2=sum2-xrow*lgamma(alpha[j]);
        }
    }
    double max_likelihood = 0.0;
    max_likelihood=score[clumpnum-1]-sum2;

    return max_likelihood;

}


