# AIT

An adaptive independence test for microbiome community data

# Installation

Firstly, Install Rtools (https://cran.r-project.org/bin/windows/Rtools/) if your computer is under Microsoft Windows;

then, 
```
install.packages("devtools")

devtools::install_github("YaruSong/AIT")
```

# Contents

## Function

### simulation

`sim_x` simulate multinomial or Dirichlet-multinomial distributed microbiome data

`sim_x_c` simulate Dirichlet-multinomial distributed microbiome data subjected to covariates

`lambda_selection` select a proper lambda to control type I error for the independence test

### test

`AIT_MN_p` an adaptive independence test for multinomial modeled microbiome data

`AIT_DMN_p` an adaptive independence test for Dirichlet-multinomial modeled microbiome data

`AIT_DMN_adj_p` an adaptive independence test for Dirichlet-multinomial modeled microbiome data with covariates adjustment

### others

`relabel` reassigning values of a vector

`lambda_selection` select a proper lambda to control type I error for the independence test.

## data

`x.m` a simulated multinomial based microbiome data

`x.d` a simulated Dirichlet-multinomial based microbiome data

`x.d.c` a simulated Dirichlet-multinomial based microbiome data subjected to covariates

`individual_index` index information for samples

`taxa_family` Taxa counts at family level for samples

# Example

```
library(AIT)
rm(list=ls())
set.seed(2019)
samples=100
groups.samples<-c(10,20,30,40)
dispersion=0.2
seq.depth=1000
mean.0<-c(rep(0,6),rep(0.005,5),rep(0.01,5),rep(0.02,4),
rep(0.03,3),rep(0.04,2),c(0.08,0.1,0.11,0.165,0.22))
mean.change<-function(pi,location,alpha){
mean.tmp=sort(pi)
mean.tmp[0+location]= mean.tmp[0+location]+alpha*mean.tmp[20+location]
mean.tmp[20+location]= (1-alpha)*mean.tmp[20+location]
return(mean.tmp)
}
beta=0.6  #beta acts as the effect size
mean.control=mean.0
mean.case1=mean.change(mean.0,c(8,9,10),beta)
mean.case2=mean.change(mean.0,c(5,6,7),beta)
mean.case3=mean.change(mean.0,c(2,3,4),beta)
mean.group<-as.matrix(rbind(mean.control,mean.case1,mean.case2,mean.case3))
x.d=sim_x("DMN",groups.samples,seq.depth,mean.group,dispersion)

y=sort(rep(1:10,10))/2
y.d=y

lambda_d= lambda_selection(method = "DMN", x.d, y.d, 1000)
print(lambda_d)
print(AIT_DMN_p(x.d, y.d, lambda_d, 1000, slice=TRUE))

```
