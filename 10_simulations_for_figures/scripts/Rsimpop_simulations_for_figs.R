## Run interactively

library(ape)
library(stringr)
library(tidyr)
library(seqinr)
library(phytools)
library(ggplot2)
library(devtools)
library(MCMCglmm)
library(phangorn)
library(spam)
library(INLA)
library(ggtree)
library(rsimpop)

#Set working directory

my_working_directory = "/lustre/scratch117/casm/team154/em16/simulations/rsimpop2/sims"

##Initialise the environment
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/summary_stats_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/foetal.filters.parallel.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/phytools_scripts.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/targeted_analysis_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/Pop_size_estimation_functions.R")
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");setwd(my_working_directory)

#Function for plotting phylodyn trajectory
plot_phylodyn=function(simtree,label=""){
  fullbnpr=BNPR(ut,lengthout=75)
  #par(mfcol=c(1,2))
  #plot_tree(ut, cex.label = 0, lwd = 0.5)
  plot_BNPR(fullbnpr,col="black",main=label, ylim = c(1E+0, 1E+9))
  abline(h=100000,col="red")
  abline(h=1000000,col="black")
}

#Specify inputs - need to edit for each individual with cut offs based on phylodyn plots

##Setup the trajectory for constant population size age 30 (Extended Fig 5b)
#-------------------------------------------------------------------------------------------------------------------
age = 30
target_pop_size= 1e5

trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
print(head(trajectory))

write.table(trajectory, "trajectory/constant_30_trajectory.txt", row.names = FALSE, quote = FALSE)

##Run the sim
sp=run_neutral_trajectory(NULL,0.5,trajectory)

#Visualise the trajectory
pdf("pdfs/constant_30_trajectory.pdf", height = 5, width = 6)
plot(sp,xlim=c(0,30),ylim= c(1E+0, 1E+8))
lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
dev.off()

## Get the tree of extant lineages
fulltree=get_tree_from_simpop(sp)
write.tree(fulltree, "trees/constant_30.tree")

st=get_subsampled_tree(fulltree,380)
write.tree(st, "trees/constant_30_st.tree")

st=get_elapsed_time_tree(st)
st$edge.length=st$edge.length/365
ut=drop.tip(st,"zeros")
ut=multi2di(ut)
class(ut)="phylo"

pdf("pdfs/constant_30_ultra.pdf", width=12,height=4, useDingbats = FALSE)
plot_tree(ut, cex.label = 0, lwd = 0.5)
dev.off()

pdf("pdfs/constant_30_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
plot_phylodyn(ut)
dev.off()

##Setup the trajectory for late-life increase in population size to age 30 to 250,000 (Extended Fig 5b)
#-------------------------------------------------------------------------------------------------------------------
age = 30
target_pop_size= 1e5
latelife = 20:age
change1=15:19

latelife_pop_change = 2.5
dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/5

trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[14]+(dif1*(1:5)),digits=0)
trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
print(head(trajectory))

write.table(trajectory, "trajectory/growth_young_2.5_trajectory.txt", row.names = FALSE, quote = FALSE)

##Run the sim
sp=run_neutral_trajectory(NULL,0.5,trajectory)

#Visualise the trajectory
pdf("pdfs/growth_young_2.5_trajectory.pdf", height = 5, width = 6)
plot(sp,xlim=c(0,30),ylim= c(1E+0, 1E+8))
lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
dev.off()

## Get the tree of extant lineages
fulltree=get_tree_from_simpop(sp)
write.tree(fulltree, "trees/growth_young_2.5.tree")

st=get_subsampled_tree(fulltree,380)
write.tree(st, "trees/growth_young_2.5_st.tree")

st=get_elapsed_time_tree(st)
st$edge.length=st$edge.length/365
ut=drop.tip(st,"zeros")
ut=multi2di(ut)
class(ut)="phylo"

pdf("pdfs/growth_young_2.5_ultra.pdf", width=12,height=4, useDingbats = FALSE)
plot_tree(ut, cex.label = 0, lwd = 0.5)
dev.off()

pdf("pdfs/growth_young_2.5_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
plot_phylodyn(ut)
dev.off()


##Setup the trajectory for late-life increase in population size to age 30 to 500,000 (Extended Fig 5b)
#-------------------------------------------------------------------------------------------------------------------
age = 30
target_pop_size= 1e5
latelife = 20:age
change1=15:19

latelife_pop_change = 5
dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/5

trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[14]+(dif1*(1:5)),digits=0)
trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
print(head(trajectory))

write.table(trajectory, "trajectory/growth_young_5_trajectory.txt", row.names = FALSE, quote = FALSE)

##Run the sim
sp=run_neutral_trajectory(NULL,0.5,trajectory)

#Visualise the trajectory
pdf("pdfs/growth_young_5_trajectory.pdf", height = 5, width = 6)
plot(sp,xlim=c(0,30),ylim= c(1E+0, 1E+8))
lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
dev.off()

## Get the tree of extant lineages
fulltree=get_tree_from_simpop(sp)
write.tree(fulltree, "trees/growth_young_5.tree")

st=get_subsampled_tree(fulltree,380)
write.tree(st, "trees/growth_young_5_st.tree")

st=get_elapsed_time_tree(st)
st$edge.length=st$edge.length/365
ut=drop.tip(st,"zeros")
ut=multi2di(ut)
class(ut)="phylo"

pdf("pdfs/growth_young_5_ultra.pdf", width=12,height=4, useDingbats = FALSE)
plot_tree(ut, cex.label = 0, lwd = 0.5)
dev.off()

pdf("pdfs/growth_young_5_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
plot_phylodyn(ut)
dev.off()


##Setup the trajectory for late-life increase in population size from 100,000 to 1,000,000 (Extended Fig 5b)
#-------------------------------------------------------------------------------------------------------------------
age =80
target_pop_size= 1e5
latelife = 61:age
change1=51:60

latelife_pop_change = 10
dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/10

trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[50]+(dif1*(1:10)),digits=0)
trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
print(head(trajectory))

write.table(trajectory, "trajectory/growth_trajectory.txt", row.names = FALSE, quote = FALSE)

##Run the sim
sp=run_neutral_trajectory(NULL,0.5,trajectory)

#Visualise the trajectory
pdf("pdfs/growth_trajectory.pdf", height = 5, width = 6)
plot(sp,xlim=c(0,80),ylim= c(1E+0, 1E+8))
lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
dev.off()

## Get the tree of extant lineages
fulltree=get_tree_from_simpop(sp)
write.tree(fulltree, "trees/growth.tree")

st=get_subsampled_tree(fulltree,380)
write.tree(st, "trees/growth_st.tree")

st=get_elapsed_time_tree(st)
st$edge.length=st$edge.length/365
ut=drop.tip(st,"zeros")
ut=multi2di(ut)
class(ut)="phylo"

pdf("pdfs/growth_ultra.pdf", width=12,height=4, useDingbats = FALSE)
plot_tree(ut, cex.label = 0, lwd = 0.5)
dev.off()

pdf("pdfs/growth_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
plot_phylodyn(ut)
dev.off()





##Setup the trajectory for constant population 10e05 size age 80 (Extended Fig 6a)
#-------------------------------------------------------------------------------------------------------------------
age = 80  
target_pop_size= 1e5

  trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
  print(head(trajectory))
  
  write.table(trajectory, "trajectory/constant_trajectory.txt", row.names = FALSE, quote = FALSE)
  
  ##Run the sim
  sp=run_neutral_trajectory(NULL,0.5,trajectory)
  
  #Visualise the trajectory
  pdf("pdfs/constant_trajectory.pdf", height = 5, width = 6)
  plot(sp,xlim=c(0,80),ylim= c(1E+0, 1E+8))
  lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
  dev.off()
  
  ## Get the tree of extant lineages
  fulltree=get_tree_from_simpop(sp)
  write.tree(fulltree, "trees/constant.tree")
  
  st=get_subsampled_tree(fulltree,380)
  write.tree(st, "trees/constant_st.tree")
  
  st=get_elapsed_time_tree(st)
  st$edge.length=st$edge.length/365
  ut=drop.tip(st,"zeros")
  ut=multi2di(ut)
  class(ut)="phylo"
  
  pdf("pdfs/constant_ultra.pdf", width=12,height=4, useDingbats = FALSE)
  plot_tree(ut, cex.label = 0, lwd = 0.5)
  dev.off()
  
  pdf("pdfs/constant_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
  plot_phylodyn(ut)
  dev.off()
  
  
  ##Setup the trajectory for bottleneck in population size between age 30-45 in 80 year old (Extended Fig 6a)
  #-------------------------------------------------------------------------------------------------------------------
  age = 80
  target_pop_size= 1e5
  midlife = 30:45
  latelife = 56:age
  change1=20:29
  change2=46:55
  
  midlife_pop_change = 0.1
  latelife_pop_change = 1
  dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/10
  dif2 = ((target_pop_size*latelife_pop_change)-(target_pop_size*midlife_pop_change))/10
  
  trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
  trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[19]+(dif1*(1:10)),digits=0)
  trajectory$target_pop_size[midlife]=round(midlife_pop_change*trajectory$target_pop_size[midlife], digits=0)
  trajectory$target_pop_size[change2]=round(trajectory$target_pop_size[45]+(dif2*(1:10)),digits =0)
  trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
  print(head(trajectory))
  
  write.table(trajectory, "trajectory/bottleneck_trajectory.txt", row.names = FALSE, quote = FALSE)
  
  ##Run the sim
  sp=run_neutral_trajectory(NULL,0.5,trajectory)
  
  #Visualise the trajectory
  pdf("pdfs/bottleneck_trajectory.pdf", height = 5, width = 6)
  plot(sp,xlim=c(0,80),ylim= c(1E+0, 1E+8))
  lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
  dev.off()
  
  ## Get the tree of extant lineages
  fulltree=get_tree_from_simpop(sp)
  write.tree(fulltree, "trees/bottleneck.tree")
  
  st=get_subsampled_tree(fulltree,380)
  write.tree(st, "trees/bottleneck_st.tree")
  
  st=get_elapsed_time_tree(st)
  st$edge.length=st$edge.length/365
  ut=drop.tip(st,"zeros")
  ut=multi2di(ut)
  class(ut)="phylo"
  
  pdf("pdfs/bottleneck_ultra.pdf", width=12,height=4, useDingbats = FALSE)
  plot_tree(ut, cex.label = 0, lwd = 0.5)
  dev.off()
  
  pdf("pdfs/bottleneck_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
  plot_phylodyn(ut)
  dev.off()
  

  
  
 
  
  
    ##Setup the trajectory for drivers at ages 20, 40, 60, 80 
    #- script run 4 times consecutively to create phylogenies in Extended Fig 10
    #-------------------------------------------------------------------------------------------------------------------
    
    
    # Add additional functions which will be incorporated into the final version of rsimpop
    plot_tree_events2=function(tree,legpos="topleft",cex.event=2,...){
      events=tree$events
      events=events[order(events$ts),]
      if(is.null(tree$events)){
        stop("tree has no events")
      }
      tree=plot_tree(tree,...)
      idx.child=match(tree$edge[,2],tree$edge[,1])
      N=dim(tree$edge)[1]
      TT=max(tree$timestamp)
      duration=ifelse(is.na(idx.child),TT-tree$tBirth,tree$tBirth[idx.child]-tree$tBirth)+1e-6
      
      
      df=data.frame(uval=unique(sort(sprintf("%s:%s",events$value,events$driverid))),stringsAsFactors = FALSE)
      
      
      cols=c(RColorBrewer::brewer.pal(9,"Set1"))#,RColorBrewer::brewer.pal(12,"Set3"))
      df$col=rep(cols,26)[1:length(df$uval)]
      df$pch=rep(c(19,17,15,25,setdiff(0:25,c(19,17,15,25))),each=length(cols))[1:length(df$uval)]
      print(df)
      
      
      events=events %>% mutate(key=sprintf("%s:%s",value,driverid)) %>% left_join(df,by=c("key"="uval"))
      events$idx=match(events$node,tree$edge[,2])
      fracs=lapply(1:N,function(x) c())
      cols =lapply(1:N,function(x) c())
      #browser()
      for(i in 1:length(events$node)){
        thisnode=events$node[i]
        idx=events$idx[i]
        info=get_edge_info(NULL,tree,thisnode)
        frac=(events$ts[i]-tree$tBirth[events$idx[i]])/duration[events$idx[i]]
        points(x=info$x,y=info$yt-frac*(info$yt-info$yb),pch=events$pch[i],col=events$col[i],cex=cex.event)
        fracs[[idx]]=c(fracs[[idx]],frac)
        cols[[idx]]=c(cols[[idx]],events$col[i])
        kids=get_all_node_children(thisnode,tree)
        for(k in match(kids,tree$edge[,2])){
          fracs[[k]]=0
          cols[[k]]=events$col[i]
        }
        
      }
      for(j in 1:N){
        info=get_edge_info(NULL,tree,tree$edge[j,2])
        ff=c(fracs[[j]][-1],1)
        segments(x0=info$x,y1=info$yt-fracs[[j]]*(info$yt-info$yb),
                 y0=info$yt-ff*(info$yt-info$yb),col=cols[[j]],lwd=1)
        
      }
      if(!is.null(legpos)){
        legend(legpos,legend=df$uval,col=df$col,pch=df$pch,cex=1)
      }
      tree
    }
    
    run_driver_process_sim_KEVIN=function(initial_division_rate = 0.1, final_division_rate = 1/365,
                                          target_pop_size = 1e+05, drivers_per_year = 0.1, nyears = 40,
                                          fitnessGen = function() {
                                            0
                                          }, bForceReseed = FALSE, offset = 0)
    {
      if (bForceReseed) {
        delay = (offset/10 + 1)
        Sys.sleep(delay)
        cat("delay=", delay, "\n")
        initSimPop(-1, bForce = TRUE)
      }
      dpcpd = (drivers_per_year/365)/target_pop_size
      k = 1
      nsd = nyears * 365
      cfg = getDefaultConfig(1e+05, rate = initial_division_rate,
                             ndriver = 1, basefit = 0)
      params = list(n_sim_days = nsd, b_stop_at_pop_size = 1, b_stop_if_empty = 0,
                    driver_rate_per_cell_per_day = dpcpd)
      growthphase = get_tree_from_simpop(sim_pop(NULL, params = params,
                                                 cfg))
      while (growthphase$status == 2 || growthphase$status == 0) {
        k = k + 1
        params = list(n_sim_days = nsd, b_stop_at_pop_size = 1,
                      b_stop_if_empty = 0, driver_rate_per_cell_per_day = dpcpd)
        if (max(growthphase$timestamp) > params$n_sim_days) {
          stop("ts>n_sim_days:unexpected bahaviour")
        }
        growthphase = addDriverEvent(growthphase, growthphase$cfg,
                                     currentCompartment = 1, fitness = fitnessGen())
        gpk = get_tree_from_simpop(sim_pop(growthphase, params = params,
                                           growthphase$cfg))
        growthphase = combine_simpops(growthphase, gpk)
      }
      tree0 = get_tree_from_simpop(growthphase)
      tree0$cfg$compartment$rate[2] = final_division_rate
      tree0$cfg$compartment$popsize[2] = target_pop_size
      years = nyears
      params[["b_stop_at_pop_size"]] = 0
      if (max(tree0$timestamp) >= nyears * 365) {
        return(tree0)
      }
      adult = sim_pop(tree0, params = params, tree0$cfg)
      adult = combine_simpops(tree0, adult)
      tree0 = get_tree_from_simpop(adult)
      if (length(which(tree0$cfg$info$fitness == 0 & tree0$cfg$info$population >
                       0)) > 2) {
        browser()
      }
      while (TRUE) {
        k = k + 1
        if (max(tree0$timestamp) >= nyears * 365) {
          return(tree0)
        }
        tree0 = addDriverEvent(tree0, tree0$cfg, 1, fitness = fitnessGen())
        tmp2 = tree0
        adult = sim_pop(tree0, params = params, tree0$cfg, b_verbose = FALSE)
        tree0 = combine_simpops(tree0, adult)
        if (length(which(tree0$cfg$info$fitness == 0 & tree0$cfg$info$population >
                         0)) > 2) {
          browser()
        }
        if (k > 1000000) {
          stop("too many iterations!")
        }
      }
    }
    
    
    
    gamma_shape = 0.73
    gamma_rate = 33
    number_drivers_per_year = 200
    fitness_threshold = 0.05
    
    ##Function to generate gamma distribution based fitness
    genGammaFitness=function(fitness_threshold,shape,rate){
      function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
    }
    fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)
    
    dps_20=run_driver_process_sim_KEVIN(0.1, final_division_rate = 0.5/(365),target_pop_size = 1e5, nyears = 20, fitness= fitnessGammaFn, drivers_per_year = number_drivers_per_year)
    
    dpst_20=get_subsampled_tree(dps_20,380)
    
    st_20=get_elapsed_time_tree(dpst_20)
    st_20$edge.length=st_20$edge.length/365
    ut_20=drop.tip(st_20,"zeros")
    ut_20=multi2di(ut_20)
    class(ut_20)="phylo"
    write.tree(ut_20, "trees/20_200_drivers8.tree")
    
    driver_counts <- dps_20$cfg$info %>% filter(population>0)
    write.table(driver_counts,"drivers/20_200_drivers8.txt")
    
    pdf("pdfs/20_200_drivers8.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree_events2(ut_20, cex.label = 0)
    dev.off()
    
    ## Continue same process to age 40
    dps_40 =continue_driver_process_sim(dps_20 ,nyears = 40 ,fitnessGen = fitnessGammaFn)
    
    dpst_40=get_subsampled_tree(dps_40,380)
    
    st_40=get_elapsed_time_tree(dpst_40)
    st_40$edge.length=st_40$edge.length/365
    ut_40=drop.tip(st_40,"zeros")
    ut_40=multi2di(ut_40)
    class(ut_40)="phylo"
    write.tree(ut_40, "trees/40_200_drivers8.tree")
    
    driver_counts <- dps_40$cfg$info %>% filter(population>0)
    write.table(driver_counts,"drivers/40_200_drivers8.txt")
    
    pdf("pdfs/40_200_drivers8.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree_events2(ut_40, cex.label = 0)
    dev.off()
    
    ## Continue same process to age 60
    dps_60 =continue_driver_process_sim(dps_40 ,nyears = 60 ,fitnessGen = fitnessGammaFn)
    
    dpst_60=get_subsampled_tree(dps_60,380)
    
    st_60=get_elapsed_time_tree(dpst_60)
    st_60$edge.length=st_60$edge.length/365
    ut_60=drop.tip(st_60,"zeros")
    ut_60=multi2di(ut_60)
    class(ut_60)="phylo"
    write.tree(ut_60, "trees/60_200_drivers8.tree")
    
    driver_counts <- dps_60$cfg$info %>% filter(population>0)
    write.table(driver_counts,"drivers/60_200_drivers8.txt")
    
    pdf("pdfs/60_200_drivers8.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree_events2(ut_60, cex.label = 0)
    dev.off()
    
    ## Continue same process to age 80
    dps_80 =continue_driver_process_sim(dps_60 ,nyears = 80 ,fitnessGen = fitnessGammaFn)
    
    dpst_80=get_subsampled_tree(dps_80,380)
    
    st_80=get_elapsed_time_tree(dpst_80)
    st_80$edge.length=st_80$edge.length/365
    ut_80=drop.tip(st_80,"zeros")
    ut_80=multi2di(ut_80)
    class(ut_80)="phylo"
    write.tree(ut_80, "trees/80_200_drivers8.tree")
    
    driver_counts <- dps_80$cfg$info %>% filter(population>0)
    write.table(driver_counts,"drivers/80_200_drivers8.txt")
    
    pdf("pdfs/80_200_drivers8.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree_events2(ut_80, cex.label = 0)
    dev.off()
    
    ## Continue same process to age 100
    dps_100 =continue_driver_process_sim(dps_80 ,nyears = 100 ,fitnessGen = fitnessGammaFn)
    
    dpst_100=get_subsampled_tree(dps_100,380)
    
    st_100=get_elapsed_time_tree(dpst_100)
    st_100$edge.length=st_100$edge.length/365
    ut_100=drop.tip(st_100,"zeros")
    ut_100=multi2di(ut_100)
    class(ut_100)="phylo"
    write.tree(ut_100, "trees/100_200_drivers8.tree")
    
    driver_counts <- dps_100$cfg$info %>% filter(population>0)
    write.table(driver_counts,"drivers/100_200_drivers8.txt")
    
    pdf("pdfs/100_200_drivers8.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree_events2(ut_100, cex.label = 0)
    dev.off()
    
    #-------------------------------------------------------------------------------------------------------------------
    
    #Function for plotting phylodyn trajectory for supplementary simulation figures (note unnecessary ablines deleted in Illustrator)
    plot_phylodyn=function(simtree,label=""){
      fullbnpr=BNPR(ut,lengthout=75)
      #par(mfcol=c(1,2))
      #plot_tree(ut, cex.label = 0, lwd = 0.5)
      plot_BNPR(fullbnpr,col="black",main=label, ylim = c(1E+0, 1E+9))
      abline(h=1000, col= "black")
      abline(h=10000, col= "red")
      abline(h=25000, col= "black")
      abline(h=50000, col= "red")
      abline(h=100000, col= "black")
      abline(h=250000,col="red")
      abline(h=500000,col="black")
    }
    
    #Specify inputs 
    
    ##Setup the trajectory for differing population sizes at age 30
    age = 30
    
    #Population size 25000 (Supplementary Figure 4)
    #-------------------------------------------------------------------------------------------------------------------
    target_pop_size= 2.5e4
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/pop_25_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/pop_25_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/pop_25.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/pop_25_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/pop_25_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/pop_25_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    # Population size 100,000 (Supplementary Figure 4)
    #-------------------------------------------------------------------------------------------------------------------
    target_pop_size= 1e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/pop_100_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/pop_100_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/pop_100.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/pop_100_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/pop_100_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/pop_100_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    # Population size 250000 (Supplementary Figure 4)
    #-------------------------------------------------------------------------------------------------------------------
    target_pop_size= 2.5e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/pop_250_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/pop_250_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/pop_250.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/pop_250_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/pop_250_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/pop_250_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    # Population size 750000 (Supplementary Figure 4)
    #-------------------------------------------------------------------------------------------------------------------
    target_pop_size= 7.5e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/pop_750_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/pop_750_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/pop_750.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/pop_750_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/pop_750_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/pop_750_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
     
    #Population size 100000 age 20 (Supplementary Figure 5)
    #-------------------------------------------------------------------------------------------------------------------
    age = 20
    target_pop_size= 1e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/constant_20_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/constant_20_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/constant_20.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/constant_20_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/constant_20_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/constant_20_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    #Population size 100000 age 40 (Supplementary Figure 5)
    #-------------------------------------------------------------------------------------------------------------------
    age = 40
    target_pop_size= 1e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/constant_40_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/constant_40_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/constant_40.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/constant_40_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/constant_40_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/constant_40_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    #Population size 100000 age 60  (Supplementary Figure 5)
    #-------------------------------------------------------------------------------------------------------------------
    age = 60
    target_pop_size= 1e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/constant_60_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/constant_60_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/constant_60.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/constant_60_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/constant_60_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/constant_60_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    #Population size 100000 age 80  (Supplementary Figure 5)
    #-------------------------------------------------------------------------------------------------------------------
    age = 80
    target_pop_size= 1e5
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/constant_80_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/constant_80_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/constant_80.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/constant_80_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/constant_80_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/constant_80_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for mid-life decrease in population size age 20 to 25000 (Supplementary Figure 6)
    #-------------------------------------------------------------------------------------------------------------------
    age = 20
    target_pop_size= 1e5
    latelife = 10:age
    change1= 8:9
    
    latelife_pop_change = 0.25
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[7]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/decline_1_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/decline_1_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/decline_1.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/decline_1_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/decline_1_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/decline_1_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for late-life decrease in population size 40 year to 25000 (Supplementary Figure 6)
    #-------------------------------------------------------------------------------------------------------------------
    age = 40
    target_pop_size= 1e5
    latelife = 20:age
    change1= 18:19
    
    latelife_pop_change = 0.25
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[17]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/decline_2_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/decline_2_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/decline_2.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/decline_2_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/decline_2_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/decline_2_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for late-life decrease in population size 60 year old to 25000 (Supplementary Figure 6)
    #-------------------------------------------------------------------------------------------------------------------
    age = 60
    target_pop_size= 1e5
    latelife = 30:age
    change1=28:29
    
    latelife_pop_change = 0.25
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[27]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/decline_3_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/decline_3_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/decline_3.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/decline_3_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/decline_3_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/decline_3_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for late-life decrease in population size age 80 to 25000 (Supplementary Figure 6)
    #-------------------------------------------------------------------------------------------------------------------
    age = 80
    target_pop_size= 1e5
    latelife = 40:age
    change1=38:39
    
    latelife_pop_change = 0.25
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[37]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/decline_4_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/decline_4_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/decline_4.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/decline_4_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/decline_4_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/decline_4_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    
    ##Setup the trajectory for late-life increase in population size age 20 to 750000 (Supplementary Figure 7)
    #-------------------------------------------------------------------------------------------------------------------
    age = 20
    target_pop_size= 1e5
    latelife = 10:age
    change1= 8:9
    
    latelife_pop_change = 7.5
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[7]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/growth_1_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/growth_1_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/growth_1.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/growth_1_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/growth_1_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/growth_1_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for late-life increase in population size 40 year to 750000 (Supplementary Figure 7)
    #-------------------------------------------------------------------------------------------------------------------
    age = 40
    target_pop_size= 1e5
    latelife = 20:age
    change1= 18:19
    
    latelife_pop_change = 7.5
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[13]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/growth_2_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/growth_2_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/growth_2.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/growth_2_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/growth_2_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/growth_2_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for late-life increase in population size 60 year old to 750000 (Supplementary Figure 7)
    #-------------------------------------------------------------------------------------------------------------------
    age = 60
    target_pop_size= 1e5
    latelife = 30:age
    change1=28:29
    
    latelife_pop_change = 7.5
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[37]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/growth_3_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/growth_3_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/growth_3.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/growth_3_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/growth_3_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/growth_3_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for late-life increase in population size age 80 to 750000 (Supplementary Figure 7)
    #-------------------------------------------------------------------------------------------------------------------
    age = 80
    target_pop_size= 1e5
    latelife = 40:age
    change1=38:39
    
    latelife_pop_change = 7.5
    dif1 = ((target_pop_size*latelife_pop_change)-target_pop_size)/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[37]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/growth_4_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/growth_4_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/growth_4.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/growth_4_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/growth_4_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/growth_4_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for bottleneck in population size age 20 to 10000 (Supplementary Figure 8)
    #-------------------------------------------------------------------------------------------------------------------
    age = 20
    target_pop_size= 1e5
    midlife = 8:12
    latelife = 15:age
    change1=6:7
    change2=13:14
    
    midlife_pop_change = 0.1
    latelife_pop_change = 1
    dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/2
    dif2 = ((target_pop_size*latelife_pop_change)-(target_pop_size*midlife_pop_change))/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[5]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[midlife]=round(midlife_pop_change*trajectory$target_pop_size[midlife], digits=0)
    trajectory$target_pop_size[change2]=round(trajectory$target_pop_size[12]+(dif2*(1:2)),digits =0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/bottleneck_1_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/bottleneck_1_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/bottleneck_1.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/bottleneck_1_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/bottleneck_1_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/bottleneck_1_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    
    ##Setup the trajectory for bottleneck in population size age 40 to 10000 (Supplementary Figure 8)
    #-------------------------------------------------------------------------------------------------------------------
    age = 40
    target_pop_size= 1e5
    midlife = 15:20
    latelife = 23:age
    change1=13:14
    change2=21:22
    
    midlife_pop_change = 0.1
    latelife_pop_change = 1
    dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/2
    dif2 = ((target_pop_size*latelife_pop_change)-(target_pop_size*midlife_pop_change))/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[12]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[midlife]=round(midlife_pop_change*trajectory$target_pop_size[midlife], digits=0)
    trajectory$target_pop_size[change2]=round(trajectory$target_pop_size[20]+(dif2*(1:2)),digits =0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/bottleneck_2_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/bottleneck_2_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/bottleneck_2.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/bottleneck_2_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/bottleneck_2_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/bottleneck_2_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    ##Setup the trajectory for bottleneck in population size age 80 to 50000 (Supplementary Figure 8)
    #-------------------------------------------------------------------------------------------------------------------
    age = 60
    target_pop_size= 1e5
    midlife = 25:30
    latelife = 33:age
    change1=23:24
    change2=31:32
    
    midlife_pop_change = 0.1
    latelife_pop_change = 1
    dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/2
    dif2 = ((target_pop_size*latelife_pop_change)-(target_pop_size*midlife_pop_change))/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[22]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[midlife]=round(midlife_pop_change*trajectory$target_pop_size[midlife], digits=0)
    trajectory$target_pop_size[change2]=round(trajectory$target_pop_size[30]+(dif2*(1:2)),digits =0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/bottleneck_3_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/bottleneck_3_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/bottleneck_3.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/bottleneck_3_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/bottleneck_3_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/bottleneck_3_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()
    
    
    ##Setup the trajectory for bottleneck in population size age 80 to 10000 (Supplementary Figure 8)
    #-------------------------------------------------------------------------------------------------------------------
    age = 80
    target_pop_size= 1e5
    midlife = 30:45
    latelife = 48:age
    change1=28:29
    change2=46:47
    
    midlife_pop_change = 0.1
    latelife_pop_change = 1
    dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/2
    dif2 = ((target_pop_size*latelife_pop_change)-(target_pop_size*midlife_pop_change))/2
    
    trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=0.5/(365))
    trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[27]+(dif1*(1:2)),digits=0)
    trajectory$target_pop_size[midlife]=round(midlife_pop_change*trajectory$target_pop_size[midlife], digits=0)
    trajectory$target_pop_size[change2]=round(trajectory$target_pop_size[45]+(dif2*(1:2)),digits =0)
    trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
    print(head(trajectory))
    
    write.table(trajectory, "trajectory/bottleneck_4_trajectory.txt", row.names = FALSE, quote = FALSE)
    
    ##Run the sim
    sp=run_neutral_trajectory(NULL,0.5,trajectory)
    
    #Visualise the trajectory
    pdf("pdfs/bottleneck_4_trajectory.pdf", height = 5, width = 6)
    plot(sp,xlim=c(0,age),ylim= c(1E+0, 1E+8))
    lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
    dev.off()
    
    ## Get the tree of extant lineages
    fulltree=get_tree_from_simpop(sp)
    write.tree(fulltree, "trees/bottleneck_4.tree")
    
    st=get_subsampled_tree(fulltree,380)
    write.tree(st, "trees/bottleneck_4_st.tree")
    
    st=get_elapsed_time_tree(st)
    st$edge.length=st$edge.length/365
    ut=drop.tip(st,"zeros")
    ut=multi2di(ut)
    class(ut)="phylo"
    
    pdf("pdfs/bottleneck_4_ultra.pdf", width=12,height=4, useDingbats = FALSE)
    plot_tree(ut, cex.label = 0, lwd = 0.5)
    dev.off()
    
    pdf("pdfs/bottleneck_4_phylodyn.pdf", width = 5, height = 4, useDingbats = FALSE)
    plot_phylodyn(ut)
    dev.off()