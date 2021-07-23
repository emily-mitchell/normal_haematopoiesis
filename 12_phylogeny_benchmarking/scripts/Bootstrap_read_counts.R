#!/software/R-3.6.1/bin/Rscript

#Author: Mike Spencer Chapman

library(stringr)
library(seqinr)
library(parallel)
library(ids)

#Settings for this individual
exp_ID="KX003"
filtering="vaf_filtered"
version="5_01"
filtering_settings="standard_rho01"
nsim=10

#Set the file paths from this
my_working_directory=paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/tree_bootstraps")
mats_and_params_file=paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/mats_and_params_",exp_ID,"_",version)
tree_data_file=paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/tree_",exp_ID,"_",version,"_",filtering_settings,".tree")
filtered_muts_data_file=paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/annotated_mut_set_",exp_ID,"_",version,"_",filtering_settings)

#Source functions needed for the script
R_function_files = list.files("/lustre/scratch119/realdata/mdt1/team154/ms56/fetal_HSC/Phylogeny_of_foetal_haematopoiesis/R_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[c(1,3)],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

system("mkdir -p mpboot_files")
system("mkdir -p final_trees")

for(i in 1:nsim) {
  print(paste("Bootstrap",i))
  
  ##SET RUN_ID AND FILEPATHS
  Run_ID = paste0("Bootstrap_",exp_ID,"_",ids::random_id(1,8))
  dna_string_file = paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/tree_bootstraps/mpboot_files/DNA_string_file_bootstrap_",exp_ID,"_",Run_ID,".fa")
  mpboot_tree_file = paste0(dna_string_file,".treefile")
  tree_file_path = paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering,"/tree_bootstraps/final_trees/Tree_",exp_ID,"_",Run_ID,".tree")
  
  #Load the reduced mats and params file from previous filtering
  load(mats_and_params_file)
  rownames(COMB_mats$NR)=rownames(COMB_mats$NV)=COMB_mats$mat$mut_ref
  
  #Reduce more stringently to remove germline mutations and recurrent artefacts (but less stringent than the final filtering)
  COMB_mats=list_subset(COMB_mats,select_vector = (filter_params$germline_pval<0.001 & filter_params$bb_rhoval>0.1)) #Pre-apply the beta-binomial filter, as boot-strapping systematically increases the count dispersion
  filter_params<-filter_params[COMB_mats$mat$mut_ref,]
  
  calculate_vaf=function(NV,NR){
    NR[NR==0]<-1
    vaf=NV/NR
    dimnames(vaf)=dimnames(NV)
    return(vaf)
  }
  
  #THE BOOTSTRAPPING STEEP - The current vaf becomes the binomial probability, and the depth is the number of draws. Do this for all mutations, sample by sample.
  vaf=calculate_vaf(NV = COMB_mats$NV,NR=COMB_mats$NR)
  NV_sim=sapply(1:ncol(vaf),function(i) {
    NV_vec=rbinom(n=nrow(vaf),size = as.numeric(COMB_mats$NR[,i]),prob = as.numeric(vaf[,i]))
  })
  dimnames(NV_sim)<-dimnames(COMB_mats$NV)
  COMB_mats$NV<-NV_sim #Now assign these new bootstrapped counts to the main matrix & continue in usual manner.
  
  ##SET FILTER PARAMETERS - just need to have "min_variant_reads" to be included in analysis for first section (to remove null mutations)
  min_variant_reads_auto = 3
  min_variant_reads_xy = 2
  
  #Remove mutations with no samples meeting variant read threshold
  if(COMB_mats$gender=="male") {
    null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto|(COMB_mats$NV >= min_variant_reads_xy & COMB_mats$mat$Chrom %in% c("X","Y"))) == 0
  } else if(COMB_mats$gender=="female") {
    null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto) == 0
  }
  
  COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
  print(paste(sum(null_remove),"mutations removed as no sample meets minimum threshold of number of variant reads"))
  
  #Single list of functions to apply over data, so that these can then be done in parallel (these functions are all sourced from the foetal.filters.parallel.R script)
  COMB_mats$NV=as.matrix(COMB_mats$NV)
  COMB_mats$NR=as.matrix(COMB_mats$NR)
  fun.list = c(pval_matrix,
               germline.binomial.filter,
               #beta.binom.filter,
               get_mean_depth,
               low_vaf_in_pos_samples_dp3,
               get_max_depth_in_pos)
  
  #Run functions in parallel
  fun.out = mclapply(fun.list, function(f) {f(COMB_mats)},mc.preschedule = F)
  
  #Assign output of the first function to the main set of matrices - this is the "PVal matrix" for the final filter
  #Then delete it from the function output
  COMB_mats$PVal <- fun.out[[1]]
  fun.out[[1]] <- NULL
  
  #The remaining output of the function is the parameters list for all except the "max_pval_in_pos" filters
  filter_params <- as.data.frame(Reduce(cbind,fun.out)) 
  
  colnames(filter_params) <- c("germline_pval",
                               #"bb_rhoval",
                               "mean_depth",
                               "pval_within_pos_dp3",
                               "max_depth_in_pos_samples")
  rownames(filter_params) <- COMB_mats$mat$mut_ref
  
  #Now run the pval in pos filter. This has to be run separately as it relies on the "PVal" element of the COMB_mats matrix set.
  filter_params$max_pval_in_pos_sample <- get_max_pval_in_pos(COMB_mats)
  
  #FILTER THE MUTSET USING SAME PARAMETERS AS FOR THE DATA
  XY_low_depth_cutoff = 4; XY_high_depth_cutoff = 20; AUTO_low_depth_cutoff = 8; AUTO_high_depth_cutoff = 40 #for real analysis
  
  #This is the main function - applies the set cut-offs to the mutation set, filtering the mutations, assigning genotypes for each to all the samples, and building the dummy dna strings for tree building.
  filtered_muts = get_filtered_mut_set(input_set_ID = Run_ID,  #the Run_ID of the unfiltered mutation set used as input - though won't account for any removal of samples from set
                                       COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                       gender = COMB_mats$gender, #patient's gender
                                       
                                       #These parameters decide whether a mutation is retained in the "true somatic mutation" set
                                       retain_muts = NULL,  #any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                       germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                       rho = NA,  #rho cutoff for the beta-binomial filter, a measure of how "over-dispersed" the counts are compared to a binomial distribution
                                       mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order 1. lower threshold for autosomes, 2. upper threshold for autosomes, 3. lower threshold for XY, 4. upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                       pval_dp2=NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads
                                       pval_dp3=0.01,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                       min_depth = c(6,4), #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                       min_pval_for_true_somatic = 0.1,   #Default: 0.1. the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                       min_vaf = NA, #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY)
                                       
                                       #These parameters decide the genotype for each sample for each "true somatic mutation".  These may be less stringent than the initial parameters.
                                       min_variant_reads_SHARED = 2,  #the minimum number of reads for samples to be assigned a positive genotype
                                       min_pval_for_true_somatic_SHARED = 0.05,  #the p-value for coming from "true somatic mutation" read distribution to be assigned a positive genotype
                                       min_vaf_SHARED = NA) #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype
  
  
  write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file)
  system(paste0("/lustre/scratch117/casm/team154/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/mpboot -s ", dna_string_file," -bb 1000"))
  
  #Import the tree into R using ape
  tree <- read.tree(mpboot_tree_file)
  tree <- drop.tip(tree,"Ancestral")
  tree <- multi2di(tree)
  tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work
  
  #ASSIGN MUTATIONS TO THE TREE USING THE TREE_MUT PACKAGE
  df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape
  
  #Get matrices in order, and run the main assignment functions
  mtr = filtered_muts$COMB_mats.tree.build$NV; mtr = as.matrix(mtr)
  depth = filtered_muts$COMB_mats.tree.build$NR; depth = as.matrix(depth)
  p.error = rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR))
  res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  treefit_pval_cutoff = 1e-3
  
  poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
  
  #Keep the known PVV which will always be poor fit
  poor_fit[which(filtered_muts$COMB_mats.tree.build$mat$mut_ref=="22-37053571-G-A")] <- FALSE
  print(sum(poor_fit))
  
  filter_poor_fit = "no"
  if(filter_poor_fit == "yes") {
    filtered_muts$COMB_mats.tree.build <- list_subset(filtered_muts$COMB_mats.tree.build,select_vector = !poor_fit)
    mtr = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
    depth = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
    p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
    res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
  }
  
  #Assign edge lengths from the res object
  tree$edge.length <- res$df$df$edge_length
  
  #Add node and pval information to the filtered_muts object
  filtered_muts$COMB_mats.tree.build$mat$node <- tree$edge[res$summary$edge_ml,2]
  filtered_muts$COMB_mats.tree.build$mat$pval <- res$summary$pval
  
  #Set up a "mut list"
  mut_list=filtered_muts$COMB_mats.tree.build$mat[,c("mut_ref","node")]
  mut_list$shared=ifelse(mut_list$node%in%1:length(tree$tip.label),"No","Yes")
  
  #Save the res file and the tree file
  write.tree(tree, file = tree_file_path)
  #write.table(mut_list,sep = "\t",quote = F,row.names = F,file=paste0("mut_list_files/mut_list_",Run_ID))
}
