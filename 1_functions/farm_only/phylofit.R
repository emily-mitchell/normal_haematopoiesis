# Functions to run phylofit
# Author: Nick Williams

# get_phylologistic_dat
# get_traj_res
# fit_clade
# shade_between
# plot_res
# get_subsampled_tree_fixed_prop
# get_sampled_mutant_clade_from_sim
# plot_selsim_infer
# add_outgroup
# do_comparison_plots

library("rsimpop")
library("rstan")
phylofit.stan="functions{

real glogistic_lpdf(real[] t, real S, real tm, real T,real N){
int n;
real ll;
n=size(t);

ll=0.0;
for(k in 2:n){
ll=ll+log(1+exp(S*(t[k-1]-T+tm)));
ll=ll-(1/(S*N))*choose(k,2)*exp(-S*(T-tm-t[k]))*(exp(S*(t[k-1]-t[k]))-1);
ll=ll+choose(k,2)*(t[k]-t[k-1])/N+log(choose(k,2))-log(N);

}
return(ll);
}
}

data{
int N;      //num of coalescences
real t[N];  //timing of coalescences
real T;
real maxT;
real minT;
real maxLN;
real minLN;
int nmutcolony; //Could be number of mutant colonies or reads.
int nwtcolony; //Could be number of wild type colonies or reads.
}

parameters {
real<lower=0.0001,upper=2> S; //instantaneous growth rate-> exp(S)-1 per Year.
real <lower=minT,upper=maxT> tm; //midpoint
real<lower=minLN,upper=maxLN> LN; //log pop size
//real<lower=0,upper=T-t[1]> ta; //Acquisition time....
}

model {
S ~ uniform(0.001,2);
tm ~ uniform(minT,maxT);
LN ~ uniform(minLN,maxLN);
t ~ glogistic(S,tm,T,pow(10,LN));
if(nmutcolony>0)
nmutcolony ~ binomial(nmutcolony+nwtcolony,1/(1+exp(-S*(T-tm))));

}
"
PHYLOFIT_STAN_MODEL=stan_model(model_code=phylofit.stan,model_name = "phylofit")

get_phylologistic_dat=function(ultratree,node,maxt=NA){
  nc=get_all_node_children(node,ultratree)
  idx=match(nc,ultratree$edge[,2])
  nh=nodeHeights(ultratree)
  T=max(nh)
  tc=T-c(unique(sort(nh[idx,1])),T)
  if(is.na(maxt)){
    maxt=2*T
  }
  list(N=length(tc),t=tc,T=T,maxT=maxt)
}

##Gets the credibility intervals
get_traj_res=function(posterior,t,ptile=c(0.025,0.5,0.975)){
  sapply(t,function(x)
    quantile(1/(1+exp(-posterior$S*(x-posterior$tm))),ptile)
  )
}

##Fits the clade with the specified ancestral branch (node)
fit_clade=function(ultratree,node,nmutcolony,nwtcolony,nchain=3,maxt=NA,minLN=4,maxLN=7){
  
  
  dat=get_phylologistic_dat(ultratree,node,maxt = maxt)
  dat$nmutcolony=nmutcolony
  dat$nwtcolony=nwtcolony
  if(nmutcolony>=0){
    tvaf=nmutcolony/(nmutcolony+nwtcolony)
    ci=binom.test(nmutcolony,nmutcolony+nwtcolony,conf.level = 0.99)$conf.int
    #maxt=dat$T-log(tvaf)/log(1.1) #Having too wide range for maxt creates problems..
    maxt=(1/0.05)*(log(1/ci[2]-1))+dat$T
    mint=(1/2)*(log(1/ci[1]))+dat$T
    if(dat$T-dat$t[1]>mint){
      mint=mint
    }
    dat$maxT=maxt
    dat$minT=mint
    cat("maxt=",maxt,"mint=",mint,"\n")
  }else{
    dat$minT=dat$T-dat$t[1]
  }
  dat$minLN=minLN
  dat$maxLN=maxLN
  stanr=sampling(PHYLOFIT_STAN_MODEL,
                 data=dat,iter = 20000,control=list(adapt_delta=0.95),chains = nchain,cores = nchain)
  list(posterior=rstan::extract(stanr),
       res=stanr,
       ultratree=ultratree,
       dat=dat
  )
}

shade_between=function(x,y1,y2,color){
  polygon(c(x, rev(x), x[1]), c(y1, rev(y2), y1[1]),
          col = color,border = NA)
}

plot_res=function(res,maxt,conf=0.95,b_add_extra=FALSE){
  ptile=c(0.5*(1-conf),0.5,0.5*(1+conf))
  xx=seq(0.01,maxt,0.1)
  traj=get_traj_res(res$posterior,xx,ptile = ptile)
  plot(NULL,xlim=c(0,maxt),ylim=c(1e-5,1),xlab="Age (Post Conception) Years",ylab="Clonal Fraction",log="y",yaxt="n")
  yticks=c(1e-5,1e-4,1e-3,1e-2,0.1,0.5,1)
  axis(side = 2,at =yticks ,labels = sprintf("%3.2g",yticks),las=1)
  shade_between(xx,traj[1,],traj[3,],"lightgrey")
  lines(xx,traj[2,],lwd=2,col="black")
  S=quantile(exp(res$posterior$S)-1,ptile)
  text(0,1,sprintf("Estimated: S=%3.2f (%3.2f-%3.2f)",S[2],S[1],S[3]),pos = 4)
  if(b_add_extra){
    points(y=rep(1e-5,length(res$dat$t)),res$dat$T-res$dat$t,pch=15,cex=0.5)
    legend("left",c("Fitted median trajectory"),col=c("black"),lwd=2)
  }
  traj
}


get_subsampled_tree_fixed_prop=function(selsim,N=100,nmut=10){
  mutants=which(selsim$edge[,2]<=length(selsim$tip.label) & selsim$driverid>0 & selsim$state==1)
  wt=which(selsim$edge[,2]<=length(selsim$tip.label) & selsim$driverid==0 & selsim$state==1)
  outgroup=which(selsim$edge[,2]<=length(selsim$tip.label) & selsim$state==0)
  tips=selsim$edge[c(sample(mutants,nmut,replace = FALSE),sample(wt,N-nmut,replace=FALSE),outgroup),2]
  get_subsampled_tree(selsim,-1,tips)
}

get_sampled_mutant_clade_from_sim=function(selsim,nmut=-1,nsample=100){
  acf=selsim$cfg$info$population[selsim$cfg$info$driver1==1]/(sum(selsim$cfg$info$population)-1)
  if(nmut<0){
    nmut=rbinom(1,nsample,acf)
    if(nmut<4){
      nmut=4 ## Require at least 4 colonies for further analysis
    }
  }
  t1=get_subsampled_tree_fixed_prop(selsim,N = nmut,nmut=nmut)
  t1b=get_elapsed_time_tree(t1)
  t1b$edge.length=t1b$edge.length/365
  node=t1$events$node[t1$events$driverid==1]
  list(ultratree=t1b,node=node,acf=acf)
}

plot_selsim_infer=function(res,selsim){
  plot_res(res,max(selsim$timestamp/365))
  lines(selsim$timestamp/365,selsim$totaldrivercount/sum(selsim$cfg$compartment$popsize),col="red",lwd=2)
  text(0,0.5,sprintf("True S=%3.2f",exp(selsim$cfg$drivers$fitness*365*selsim$cfg$compartment$rate[2])-1),pos=4)
  points(y=rep(1e-5,length(res$dat$t)),res$dat$T-res$dat$t,pch=15,cex=0.5)
  legend("left",c("True trajectory","Fitted trajectory"),col=c("red","black"),lwd=2)
}

add_outgroup=function(tree,num.shared.var=200){
  test=tree
  test=bind.tree(test,read.tree(text="(zeros:0);"))
  test=root(test,outgroup = "zeros",resolve.root = TRUE)
  test=multi2di(test)
  root.idx=which(test$edge[,1]==(length(test$tip.label)+1))
  clone.idx=root.idx[which(test$edge[root.idx,2]!=match("zeros",test$tip.label))]
  if(length(clone.idx)!=1){
    stop("Unexpected number of non-zero direct ancestors of root")
  }
  test$edge.length[clone.idx]=num.shared.var
  test=di2multi(test)
  test
}

do_comparison_plots=function(r1,r2,lab1,lab2,xxlim=c(0,4),yylim=c(0,4)){
  r2=r2[match(r1$clade,r2$clade),]
  if(length(which(r1$clade!=r2$clade))>0){
    stop("Alignment error")
  }
  #par(mfcol=c(2,1))
  coldf=data.frame(type=c("s","d","a","da"),
                   col=RColorBrewer::brewer.pal(n = 4,"Set1"),
                   desc=c("simple","descendent","ancestral","descendent+ancestral"),stringsAsFactors = FALSE)
  if(is.null(r1$type)){
    r1$col="black"
  }else{
    r1$col=coldf$col[match(r1$type,coldf$type)]
  }
  plot(r1$S,r2$S,xlab=lab1,ylab=lab2,pch=19,col=r1$col,xlim=xxlim,ylim=yylim,main=sprintf("%s vs %s",lab2,lab1))
  segments(x0=r1$S,y0=r2$S_lb,y1=r2$S_ub,lwd=0.8,col="grey")
  segments(y0=r2$S,x0=r1$S_lb,x1=r1$S_ub,lwd=0.8,col="grey")
  text(r1$S,r2$S,label=r1$clade,cex=0.6,pos = 4);abline(a=0,b=1)
  if(!is.null(r1$type)){
    legend("topleft",coldf$desc,col=coldf$col,pch=19)
  }
  
}