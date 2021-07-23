require("rsimpop")
require("rstan")
PHYLOFIT_NITER=10000
run_benchmark_sim=function(S=0.4,N=5e4,divrate_per_year=1,nyears_driver_acquisition = 10,nyears = 33,minprop=0.05){
  run_selection_sim(0.1,
                    final_division_rate = divrate_per_year/365,
                    target_pop_size = N,
                    nyears_driver_acquisition = nyears_driver_acquisition,
                    nyears = nyears,fitness=log(1+S)/(divrate_per_year),
                    minprop = minprop)
}

get_subtree=function(selsim,ncolony){
  if(dim(selsim$cfg$info)[1]!=3){
    stop("unexpected selsim result!")
  }
  nmutcolony=-1
  kk=0
  ##browser()
  while(nmutcolony<2 && kk<100){
    st=get_subsampled_tree(selsim,ncolony)
    idx=which(st$events$driverid==1)
    kk=kk+1
    if(length(idx)!=1){
      #stop("No mutant clade")
      next()
    }
    node=st$events$node[idx]
    nmutcolony=length(get_samples_in_clade(node,st))
    cat("node=",node,"\n")
  }
  #if(nmutcolony<2){
  #  return(NULL)
  #}
  acf=selsim$cfg$info$population[3]/(selsim$cfg$info$population[2]+selsim$cfg$info$population[3])
  #st=get_subsampled_tree(selsim,ncolony)
  idx=which(st$events$driverid==1)
  if(length(idx)!=1){
    stop("No mutant clade")
  }
  node=st$events$node[idx]
  ts=st$events$ts[idx]
  idx=which(st$edge[,2]==node)
  nmutcolony=length(get_samples_in_clade(node,st))
  nwtcolony=length(st$tip.label)-1-nmutcolony ##subtract 1 for the outgroup
  T=max(st$timestamp)/365
  
  #Nhsc=selsim$cfg$compartment$popsize[2]
  #if(b.use.perfect.priors){
  #  LNmin=log(0.99*Nhsc)
  #  LNmax=log(1.01*Nhsc)
  #}
  
  st=get_elapsed_time_tree(st)
  st$edge.length=st$edge.length/365
  list(st=st,nmutcolony=nmutcolony,nwtcolony=nwtcolony,node=node,acf=acf)
}

phylofit.stan="functions{

real glogistic_lpdf(real[] t, real S, real tm, real T,real N,real offset){
int n;
real ll;
n=size(t);

ll=0.0;
for(k in 2:n){
ll=ll+log(1+exp(S*(t[k-1]-offset-T+tm)));
ll=ll-(1/(S*N))*choose(k,2)*exp(-S*(T-tm-t[k]+offset))*(exp(S*(t[k-1]-t[k]))-1);
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
real mutPerYear;
real smax;
int nmutcolony; //Could be number of mutant colonies or reads.
int nwtcolony; //Could be number of wild type colonies or reads.
#real approxMutRate;
}

parameters {
real<lower=0.0001,upper=smax> s; //instantaneous growth rate-> exp(s)-1 per Year. log(2)  0.6931472 
real <lower=minT,upper=maxT> tm; //midpoint
real<lower=minLN,upper=maxLN> LN; //log pop size
real<lower=-t[1],upper=t[N-1]> offset; //Include poisson variation in base node
//real<lower=0,upper=T-t[1]> ta; //Acquisition time....
}


model {
s ~ uniform(0.001,smax);
tm ~ uniform(minT,maxT);
LN ~ uniform(minLN,maxLN);
offset~normal(0,sqrt(mutPerYear*(T-t[1]))/mutPerYear);
t ~ glogistic(s,tm,T,pow(10,LN),offset);
if(nmutcolony>0)
nmutcolony ~ binomial(nmutcolony+nwtcolony,1/(1+exp(-s*(T-tm))));

}

generated quantities {
real S;
S=exp(s)-1;
}
"
PHYLOFIT_STAN_MODEL=stan_model(model_code=phylofit.stan,model_name = "phylofit")

get_phylologistic_dat=function(ultratree,node,maxt=NA,mutperyear=20){
  nc=get_all_node_children(node,ultratree)
  idx=match(nc,ultratree$edge[,2])
  nh=nodeHeights(ultratree)
  T=max(nh)
  tc=T-c(unique(sort(nh[idx,1])),T)
  if(is.na(maxt)){
    maxt=2*T
  }
  list(N=length(tc),t=tc,T=T,maxT=maxt,mutPerYear=mutperyear)
}

##Fits the clade with the specified ancestral branch (node)
fit_clade=function(ultratree,node,nmutcolony,nwtcolony,nchain=3,maxt=NA,minLN=4,maxLN=6,mutperyear=20,maxSYear=1,niter=PHYLOFIT_NITER,stan.control=list(adapt_delta=0.95)){
  dat=get_phylologistic_dat(ultratree,node,maxt = maxt,mutperyear = mutperyear )
  dat$nmutcolony=nmutcolony
  dat$nwtcolony=nwtcolony
  dat$smax=log(1+maxSYear)
  dat$smin=0.001
  if(nmutcolony>=0){
    tvaf=nmutcolony/(nmutcolony+nwtcolony)
    ##Get reasonable range for mid-point prior (tm).  
    #Firstly get range of likely true ACF based on colony sampling.
    ci=binom.test(nmutcolony,nmutcolony+nwtcolony,conf.level = 0.999)$conf.int
    #Assuming a maximum instantaneous growth rate (s) of 2 and a minimum of 0.05.   
    #Then solve ACF=1/(1+exp(-s(t-tm)))
    #Giving tm=(1/s)*log((1/ACF)-1)+t
    #The following considers the 4 extremal possibilities of minACF,maxACF x minS,maxS 
    #and uses the range as the bounds for the uniform prior of tm. 
    vlog1=log((1/ci[2])-1)
    vlog2=log((1/ci[1])-1)
    ismin=1/0.05
    ismax=1/2
    tmrange=c(ismin*vlog1,ismax*vlog1,ismin*vlog2,ismin*vlog2)
    mint=max(min(tmrange)+dat$T,0)
    maxt=max(max(tmrange)+dat$T,10)
    if(mint>maxt){
      stop("Inconsistency in prior range for tm!")
    }
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
                 data=dat,iter = niter,control=stan.control,chains = nchain,cores = nchain)
  list(posterior=rstan::extract(stanr),
       res=stanr,
       ultratree=ultratree,
       dat=dat
  )
}

get_phylofit_summary=function(zzz,b.full=FALSE,extra=NULL){
  if(b.full){
    res=zzz$res
  }else{
    res=NULL
  }
  list(
    S=quantile(zzz$posterior$S,prob=c(0.025,0.25,0.5,0.75,0.975)),
    Smean=mean(zzz$posterior$S),
    LN=quantile(zzz$posterior$LN,prob=c(0.025,0.25,0.5,0.75,0.975)),
    LNmean=mean(zzz$posterior$LN),
    dat=zzz$dat,
    ndivt=get_num_divergent(zzz$res),
    lowbfmichains=get_low_bfmi_chains(zzz$res),
    res=res,
    extra=extra
  )
}