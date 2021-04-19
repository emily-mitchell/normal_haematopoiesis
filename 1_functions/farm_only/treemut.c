#include <stdio.h>
#include <R.h>
#include <assert.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdint.h>
/*
res = .C("likelihood_mtr", 
         mtr=as.integer(mtr),
         depth=as.integer(depth),
         geno=as.integer(geno),
         el=as.double(el),
         p.error=as.double(p.error),
         nmuts=as.integer(nmuts),
         nsamp=as.integer(nsamp),
         nbranch=as.integer(nbranch),
         lik=as.double(nmuts*nbranch)
)
 */
//Error checking (file existence etc) is assumed to have been carried out in the calling R code
void likelihood_mtr(int *mtr,int *depth,int * geno, double * el, double * perr, int * nmuts, int * nsamp, int * nbranch,double * lik){
  int nm= *nmuts;
  int ns= *nsamp;
  int nb= *nbranch;
  double loglik=0.0;
  double logV[ns];
  double logP[ns];
  double log1minusV[ns];
  double log1minusP[ns];
  double logel[nb];
  int i,j,k;
  
 
  for(i=0;i<ns;i++){
    logV[i]=log(0.5-perr[i]);
    logP[i]=log(perr[i]);
    log1minusV[i]=log(0.5+perr[i]);
    log1minusP[i]=log(1-perr[i]);
  }
  for(i=0;i<nb;i++){
    logel[i]=log(el[i]);
  }
  for(i=0;i<nm;i++){
    if(i % 1000 ==0 ){
      printf("progress: done %d\n",i);
    }
    for(j=0;j<nb;j++){
      loglik=logel[j];
      for(k=0;k<ns;k++){
        if(geno[k*nb+j]==1){
          //INSIDE
          loglik+=mtr[k*nm+i]*logV[k]+(depth[k*nm+i]-mtr[k*nm+i])*log1minusV[k];
        }else{
          loglik+=mtr[k*nm+i]*logP[k]+(depth[k*nm+i]-mtr[k*nm+i])*log1minusP[k];
        }
      }
      lik[j*nm+i]=loglik;
    }
  }
}

double cpv(int tm,int * depth,double * probs,int n){
  
  int m;
  if(n==1){
    return pbinom( tm, depth[0],probs[0],1,0);
  }
  int D=tm>depth[0]?depth[0]:tm;
  double p=0;
  for(m=0;m<=D;m++){
    p+=dbinom(m,depth[0],probs[0],0)*cpv(tm-m,depth+1,probs+1,n-1);
  }
  return p;
}


double cpv2(int tm,int * depth,double * probs,int n){
  int m=0;
  if(n==1){
    return pbinom( tm, depth[0],probs[0],0,0)+dbinom( tm, depth[0],probs[0],0);
  }
  double p=0;
  int remaining_depth=0;
  for(int i=1;i<n;i++){
    remaining_depth+=depth[i];
  }
  int D=(tm-remaining_depth)>0?(tm-remaining_depth):0;
  for(m=D;m<=depth[0];m++){
    p+=dbinom(m,depth[0],probs[0],0)*cpv2(tm-m,depth+1,probs+1,n-1);
  }
  return p;
}



void cumulate_binomial_pval(int * totalmtr,int * depth,double * probs,int * n,int * lower_tail,double *p){
  //int m;
  int N=*n;
  int tm=*totalmtr;
//double ans;
  if(*lower_tail==1){
    *p=cpv(tm,depth,probs,N);
  }else{
    *p=cpv2(tm,depth,probs,N);
  }
  return;
}
  



