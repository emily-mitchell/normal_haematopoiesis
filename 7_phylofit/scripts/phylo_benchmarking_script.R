args = commandArgs(trailingOnly=TRUE)#
if(length(args)<4){
  stop("Rscript do_phylof_benchmarking_supp.R <hnsc> <dpy> <S> <ncolony> <simid>")
}
source("phylofit_for_bench.R")
N=as.numeric(args[1])
dpy=as.numeric(args[2])
S=as.numeric(args[3])
NC=as.integer(args[4])
sid=as.integer(args[5])
OUTDIR="../sims_phylofit_final_wn/"
initSimPop(seed = as.integer(N+sid),bForce = TRUE)
nyd=5
s=log(1+S)
#Following has expected ACF=0.1
nyear=nyd+ceiling(log(N*0.1*s/(dpy+s))/s)
label=sprintf("sim%d_N%s_D%3.2f_S%3.2f_C%d",sid,N,dpy,S,NC)
maxSYear=5
PHYLOFIT_NITER=10000
N1=4
N2=6
z1=lapply(1:5,function(i){
  selsim=run_benchmark_sim(S=S,divrate_per_year = dpy,N = N,nyears_driver_acquisition = nyd,nyear=nyear,minprop = 0.02)
  st=get_subtree(selsim,NC)
  lapply(list(
    with_acf_tightN_tightS=fit_clade(st$st,st$node,st$nmutcolony,st$nwtcolony,nchain=3,maxt=NA,minLN=N1,maxLN=N2,stan.control=list(adapt_delta=0.99,max_treedepth = 15),maxSYear = 1,mutperyear = 1000),
    no_acf=fit_clade(st$st,st$node,-1,-1,nchain=3,maxt=NA,minLN=N1,maxLN=N2,stan.control=list(adapt_delta=0.99,max_treedepth = 15),maxSYear = 1,mutperyear = 10000)
  ),get_phylofit_summary,b.full=FALSE,extra=list(acf=st$acf))
  
})
saveRDS(z1,sprintf("%s/%s.RDS",OUTDIR,label))

