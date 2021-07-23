####################################################################################################

# Modelling functions

# Author: Kevin Dawson
####################################################################################################

# dps = tree for whole simulated population
# st = subsampled tree

####################################################################################################
#
# Selection coefficient of top 5 drivers
# This can be interpreted as:
#
# Top 5 individual drivers  
#  simtree phylo/simpop. 
# This can be a sub-sampled tree or the full population
get_top_drivers=function(simtree,top.n=5){
  drivers=with(simtree$cfg,drivers[match(simtree$events$driverid[which(simtree$events$driverid>0)],drivers$driver),])
  head(drivers[order(drivers$fitness,decreasing = TRUE),],top.n)
}
####################################################################################################
#
# Selection coefficient of top 5 drivers
# This can be interpreted as:
#
# Top 5 clades. Here a clade can have more than one driver, currently the selection is calculated additively we could make it 
# min.pop=100
# top.n=5
get_top_clades=function(simtree,top.n=5,min.pop=1){
  inf=simtree$cfg$info %>% filter(population>=min.pop & fitness>0)
  ##Every clade its has "id" set to the most recently acquired driver.
  head(inf[order(inf$fitness,decreasing = TRUE),1:4],top.n)
}
####################################################################################################
#
# Number drivers on shared branches
#
# Number drivers on private branches
#
# Total number of drivers in tree
#
get_driver_counts=function(dps){
  drivernodes=dps$events$node[which(dps$events$driverid>0)]
  shared_drivers=length(setdiff(drivernodes,1:length(dps$tip.label)))
  private_drivers=length(intersect(drivernodes,1:length(dps$tip.label)))
  total_drivers=shared_drivers+private_drivers
  list(drivercount_on_shared=shared_drivers,drivercount_on_private=private_drivers,drivercount_total=total_drivers)
}
####################################################################################################
 
# Temporary fix for rsimpop limit on number of drivers in population
###########################################################################################################

# continue_driver_process_sim=function(insim,nyears,fitnessGen){
continue_driver_process_sim_KEVIN=function(insim,nyears,fitnessGen){
  tree0=insim
  params=insim$params
  params$maxt=NULL
  params$n_sim_days=nyears*365
  if(length(which(tree0$cfg$info$fitness==0 & tree0$cfg$info$population>0))>2){
    browser()
  }
  k=1
  while(TRUE){
    k=k+1
    ## Deal with end of last sim and next driver being too close together. Very rarely kicks in.
    ##params[["n_sim_days"]]=min(nyears*365,max(driverdt[k],max(tree0$timestamp)+0.1)) ##years of simulation
    if(max(tree0$timestamp)>=nyears*365){
      return(tree0)
    }
    tree0=addDriverEvent(tree0,tree0$cfg,1,fitness=fitnessGen())
    tmp2=tree0
    adult=sim_pop(tree0,params=params,tree0$cfg,b_verbose = FALSE)
    tree0=combine_simpops(tree0,adult)
    #tree0=get_tree_from_simpop(adult)
    #print(tree0$cfg$info %>% dplyr::filter(population>0))
    if(length(which(tree0$cfg$info$fitness==0 & tree0$cfg$info$population>0))>2){
      browser()
    }
    if(k > 1000000){
      stop("too many iterations!")
    }
  }
}

###########################################################################################################

# run_driver_process_sim=function (initial_division_rate = 0.1, final_division_rate = 1/365,
#    target_pop_size = 1e+05, drivers_per_year = 0.1, nyears = 40,
#    fitnessGen = function() {
#        0
#    }, bForceReseed = FALSE, offset = 0)
# {
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

###########################################################################################################


