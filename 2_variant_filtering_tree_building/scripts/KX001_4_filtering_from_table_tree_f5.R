#This script is designed to be interactive & quick to allow rapid exploring of filtering parameters
#run interactively on farm
#bsub -Is -G team78-grp -q yesterday -R 'select[mem>=30000] rusage[mem=30000]' -M30000 R

#module load vagrent
#bsub to farm5
#ensure working directory correct at 2 points in script

#Specify options to parse from the command line
args = commandArgs(TRUE)

##SET RUN_ID AND FILEPATHS
Run_ID = as.character(args[1])
output_directory = toString(args[2])

library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(data.table)
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/foetal.filters.parallel.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/phytools_scripts.R")

#Set file paths for saved files
ID = "KX001" #Edit
Iteration = "KX001_4" #Edit
Run_ID = "KX001_4_01" #Edit
filtering_ID = "standard_rho01" #Edit
setwd(paste0("/lustre/scratch117/casm/team154/em16/",ID,"/",Iteration,"/trees/")) #Ensure this directory exists and copy mats_and_params file there
mats_and_param_file = paste0("mats_and_params_", Run_ID)
filtered_muts_file = paste0("filtered_muts_",Run_ID,"_",filtering_ID)
dna_string_file = paste0("Filter_", Run_ID,"_", filtering_ID,".fa")
mpboot_tree_file = paste0(dna_string_file,".treefile")
res_file = paste0("res_file_",Run_ID,"_",filtering_ID)
tree_file_path = paste0("tree_", Run_ID,"_",filtering_ID, ".tree")
file_annot = paste0("annotated_mut_set_", Run_ID,"_",filtering_ID)


#If have previous run analysis - load up files
previously_run = FALSE #Edit
if(previously_run) {load(file_annot); load(tree_file_path)}

#Load up the mats and params file
load(mats_and_param_file)

#Can specify coverage cut-off here, to exclude low coverage samples from analysis - samples too difficult to even place on tree
min_sample_mean_cov = 0
other_samples_to_remove = NULL


#Remove the low coverage samples and their private mutations
if(min_sample_mean_cov > 0) {
  output = remove_low_coverage_samples(COMB_mats = COMB_mats,
                                       filter_params = filter_params,
                                       min_sample_mean_cov = min_sample_mean_cov,
                                       other_samples_to_remove = other_samples_to_remove,
                                       min_variant_reads_auto = 3, #these parameters are to remove mutations from the matrix that are no longer positive in any samples
                                       min_variant_reads_xy = 2)
  COMB_mats= output$COMB_mats
  filter_params = output$filter_params
}

#REVIEW MEAN DEPTH HISTOGRAMS TO DECIDE MANUAL CUT-OFFS FOR EXCLUDING OUTLIERS
#hist(filter_params$mean_depth, breaks = 100, xlim = c(0,60))
#hist(filter_params$mean_depth[COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 200, xlim = c(0,30))
#hist(filter_params$mean_depth[!COMB_mats$mat$Chrom %in% c("X","Y")], breaks = 200, xlim = c(0,30))
XY_low_depth_cutoff = 4; XY_high_depth_cutoff = 20; AUTO_low_depth_cutoff = 8; AUTO_high_depth_cutoff = 40

#Get the filtered mutation set - returns object with the filtered mat/ NV/NR matrices, as well as a full genotype matrix, and matrix of shared mutations only
filtered_muts = get_filtered_mut_set(input_set_ID = Run_ID,  #the Run_ID of the unfiltered mutation set used as input - just gets stored with output as a record
                                     COMB_mats = COMB_mats,  #the main full mutation matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                     filter_params = filter_params,  #the filter_params matrix as outputed by the "HSC_filtering_treebuild_table.R" script
                                     gender = COMB_mats$gender, #patients gender. Must be "male" or "female". This is determined in the "HSC_filtering_treebuild_table.R" script
                                     retain_muts = NA,  #the mut_refs (i.e. Chr-Pos-Ref-Alt) of any mutations that should be manually retained, despite not meeting filtering criteria, NULL by default
                                     
                                     #PARAMETERS FOR FILTERING THE MUTATION SET. If parameter is set to null, filter will not be applied.
                                     germline_pval = -10,  #the log10 p-value cutoff for mutations coming from an expected germline distribution
                                     rho = 0.1,  #rho cutoff for the beta-binomial filter
                                     mean_depth = c(AUTO_low_depth_cutoff,AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff),   #Numeric vector of length 4 defining mean depth at mutation site cut-offs. This is in the order (1) lower threshold for autosomes, (2) upper threshold for autosomes, (3) lower threshold for XY, (4) upper threshold for XY. This removes mis-mapping/ low reliability loci.
                                     min_depth = NA, #Numeric vector of length 2 defining minimum depths that at least one positive sample must have for mutation to be retained (AUTO and XY)
                                     min_vaf = c(0.2,0.4), #Numeric vector of length 2 defining minimum vaf in at least one sample for mutation to be retained (AUTO and XY). Mutations with value > cutoff are retained.
                                     
                                     #These filters are only appropriate for pure clonal samples i.e. single-cell colonies
                                     pval_dp2 = NA,  #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 2 reads
                                     pval_dp3 = 0.001,   #the p-value cut-off if using the "pval within pos" filter, with positive samples defined as having >= 3 reads (allows for more index hopping)
                                     min_pval_for_true_somatic = NA,   #FOR  PURE CLONAL SAMPLES ONLY (not LCM): the minimum p-value that at least one sample must have for the variant:normal read distribution coming from that expected for a true somatic
                                     
                                     #PARAMETERS FOR DECIDING GENOTYPE FOR EACH SAMPLE FOR EACH RETAINED MUTATION - should be equal to, or more relaxed than the above. At least one parameter must be used.
                                     min_variant_reads_SHARED = 2,  #the minimum number of reads for subsequent samples to be assigned a positive genotype
                                     min_pval_for_true_somatic_SHARED = NA, #FOR PURE CLONAL SAMPLES ONLY (not LCM): the minimum p-value for coming from "true somatic mutation" read distribution for subsequent samples to be assigned a positive genotype
                                     min_vaf_SHARED = c(0.15,0.3)) #Numeric vector of length 2, defining minimum vaf to be assigned a positive genotype

#Decide an ID for this filtered set, depending on approach taken, and save
save(filtered_muts, file = filtered_muts_file)
#load(filtered_muts_file)

#Write a fasta file of the dummy DNA strings
write.fasta(filtered_muts$dna_strings, names=names(filtered_muts$dna_strings), dna_string_file)

#BUILD TREE with MPBoot
cat("Starting tree building with MPBoot/n")
system(paste0("/lustre/scratch117/casm/team154/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/mpboot -s ", dna_string_file," -bb 1000"))

#Import the tree into R using ape
tree <- read.tree(mpboot_tree_file)
tree <- drop.tip(tree,"Ancestral")
tree <- multi2di(tree)
tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work

#Assign mutations back to the tree
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R"); setwd(paste0("/lustre/scratch117/casm/team154/em16/",ID,"/",Iteration,"/trees/")) #R scripts have to be sourced within containing directory or outputs an error.
df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape

#Get matrices in order, and run the main assignment functions
mtr = as.matrix(filtered_muts$COMB_mats.tree.build$NV)
depth = as.matrix(filtered_muts$COMB_mats.tree.build$NR)
p.error = c(rep(0.01, ncol(filtered_muts$COMB_mats.tree.build$NR)))
res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
treefit_pval_cutoff = 1e-3
poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
sum(poor_fit)

#Remove poor_fit muts if they look dodgy
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
#Add node information to the filtered_muts object
filtered_muts$COMB_mats.tree.build$mat$node <- tree$edge[res$summary$edge_ml,2]

early_muts = filtered_muts$COMB_mats.tree.build$mat$mut_ref[res$summary$edge_ml %in% which(nodeHeights(tree = tree)[,2] < 20)]
length(early_muts)

#Save the res file and the tree file

save(res, file = res_file)
#load(res_file)

write.tree(tree, file = tree_file_path)

#Write vcf files for VariantCaller analysis - separate out "All mutations", "Shared mutations", "Private mutations"
create_vcf_files = function(mat, select_vector = NULL) {
  if(is.null(select_vector)) {vcf_file = mat[,2:5]} else {vcf_file = mat[select_vector,2:5]}
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  vcf_file = vcf_file[,c(1,2,8,3,4,7,6,5)]
  return(vcf_file)
}

vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat)
early_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat, select_vector = res$summary$edge_ml %in% which(nodeHeights(tree)[,2] < 20))
shared_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat, select_vector = filtered_muts$COMB_mats.tree.build$mat$mut_ref %in% rownames(filtered_muts$Genotype_shared_bin))
private_vcf_file = create_vcf_files(filtered_muts$COMB_mats.tree.build$mat, select_vector = !filtered_muts$COMB_mats.tree.build$mat$mut_ref %in% rownames(filtered_muts$Genotype_shared_bin))

#Write files
write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("mutations_all_", Run_ID, ".vcf"), row.names = FALSE)
write.table(early_vcf_file, sep = "\t", quote = FALSE, file = paste0("mutations_early_", Run_ID, ".vcf"), row.names = FALSE)
write.table(shared_vcf_file, sep = "\t", quote = FALSE, file = paste0("mutations_shared_", Run_ID, ".vcf"), row.names = FALSE)
write.table(private_vcf_file, sep = "\t", quote = FALSE, file = paste0("mutations_private_", Run_ID, ".vcf"), row.names = FALSE)

#Run VAGRENT on mutation set to annotate the filtered mutations
vcf_header_path = "/lustre/scratch117/casm/team154/em16/R_scripts/VCF_header_for_VaGrent.txt"
vcf_path = paste0("mutations_all_", Run_ID, ".vcf")
vagrent_input_path = paste0("mutations_all_", Run_ID, "_header.vcf")
vagrent_output_path = paste0(vagrent_input_path,".annot")

#1. paste vcf file to a dummy header file
system(paste0("cat ",vcf_header_path," ",vcf_path," > ", vagrent_input_path))
#2. commands to run vagrent
system(paste0("AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp Human -as NCBI37 -c /lustre/scratch117/casm/team154/em16/R_scripts/vagrent.cache.gz"))
#3. import vagrent output
vagrent_output = fread(vagrent_output_path,skip = "#CHROM")
annot_info = as.data.frame(str_split(vagrent_output$INFO, pattern = ";",simplify = TRUE), stringsAsFactors = FALSE)
colnames(annot_info) <- c("VT","VD","VC","VW")

annot_info$VC <- gsub(x=annot_info$VC, pattern = "VC=", replacement = "")
annot_info$VT <- gsub(x=annot_info$VT, pattern = "VT=", replacement = "")
annot_info$VW <- gsub(x=annot_info$VW, pattern = "VW=", replacement = "")
annot_info$VD <- gsub(x=annot_info$VD, pattern = "VD=", replacement = "")


#Attempt to functionalize neatly
split_vagrent_output = function(df,split_col,col_IDs = c("Gene","Transcript","RNA","CDS","Protein","Type","SO_codes")) {
  col = df[[split_col]]
  output = matrix(nrow = nrow(df), ncol = length(col_IDs))
  for(i in 1:length(col_IDs)) {
    output[,i] = str_split(col, pattern = "\\|", simplify = TRUE)[,i]
  }
  colnames(output) = col_IDs
  return(as.data.frame(output))
}

filtered_muts$COMB_mats.tree.build$mat <- cbind(filtered_muts$COMB_mats.tree.build$mat,split_vagrent_output(df = annot_info,split_col = "VD"))


#Save the annotated filtered_muts files (post tree filtering)
save(filtered_muts, file = file_annot)