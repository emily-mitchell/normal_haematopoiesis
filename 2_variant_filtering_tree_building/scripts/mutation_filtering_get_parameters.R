#Set up matrices with information about each variant 
#This script is submitted to the Sanger farm
#Allows subsequent filtering of variants using the "filtering_from_table_tree" scripts

#Author: Mike Spencer Chapman

#!/software/R-3.6.1/bin/Rscript
library(stringr)
library(seqinr)
library(optparse)
library(parallel)

source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/filters_parallel_functions.R") #This should point to where ever this script is saved

#Specify options to parse from the command line
option_list = list(
  make_option(c("-r", "--runid"), action="store", default='HSPC_filter', type='character', help="Run ID for this filtering run"),
  make_option(c("-s", "--snvfile"), action="store", type='character', help="File path for the merged SNV output matrix from cgpVAF"),
  make_option(c("-i", "--indelfile"), action="store", default=NULL, type='character', help="File path for the merged indel output matrix from cgpVAF"),
  make_option(c("-c", "--covcut"), action="store", default=0, type='numeric', help="Remove samples with mean coverage below specified cut-off"),
  make_option(c("-o", "--output"), action="store", default=NULL, type='character', help="output directory for files"),
  make_option(c("-m", "--mixed_remove"), action="store_true", default=FALSE, type='logical', help="option to automatically check for and remove mixed colonies based on VAF plots. This will also print the VAF plots for review."),
  make_option(c("-x", "--exclude_samples"), action="store", default=NULL, type='character', help="option to manually exclude certain samples from the analysis")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

##SET RUN_ID AND FILEPATHS
Run_ID = opt$r
filepath_SNV = opt$s
filepath_INDEL = opt$i
#Mean coverage cut-off to include samples in analysis (mean coverage calculated from loci in bed-file only)
min_sample_mean_cov = opt$c
output_directory = opt$o
exclude_samples = opt$x

#Set the output directory to the current directly if it is null
if(is.null(output_directory)) {output_directory<-getwd()}

#Make the output directory if it doesn't exist
system(paste0("mkdir -p ",output_directory))

##SET UP MATRICES FOR FILTERING
#Import cgpVAF Caveman output matrix
cat("Importing SNV matrix\n")
mat_SNV <- read.delim(filepath_SNV, stringsAsFactors = FALSE)
mat_SNV <- mat_SNV[, !grepl("PDv37is", colnames(mat_SNV))] #Remove reference genome columns

#Remove samples that have a mean coverage below specified cut-off (if specified)
if(min_sample_mean_cov > 0|!is.null(exclude_samples)) {
  mean_sample_cov = colMeans(mat_SNV[,grepl("DEP", colnames(mat_SNV))])
  low_cov_names = str_split(names(mean_sample_cov)[mean_sample_cov < min_sample_mean_cov], pattern = "_", simplify = TRUE)[,1]
  remove_cols <- which(str_split(colnames(mat_SNV), pattern = "_", simplify = TRUE)[,1] %in% c(low_cov_names,exclude_samples))
  mat_SNV <- mat_SNV[,-remove_cols]
}

mat_SNV$mut_ref <- paste(mat_SNV$Chrom, mat_SNV$Pos, mat_SNV$Ref, mat_SNV$Alt, sep = "-") #Create mut_ref column
mat_SNV <- mat_SNV[,c(ncol(mat_SNV), 1:(ncol(mat_SNV) -1))] #Swap mut_ref column to the beginning of matrix
mat_SNV$Mut_type = "SNV"
NV_SNV <- mat_SNV[,grepl("MTR", colnames(mat_SNV))]
NR_SNV <- mat_SNV[,grepl("DEP", colnames(mat_SNV))]
mat_SNV <- mat_SNV[,!grepl("MTR", colnames(mat_SNV)) & !grepl("DEP", colnames(mat_SNV))]

SNV_mats <- list(mat_SNV, NV_SNV, NR_SNV)
names(SNV_mats) = c("mat", "NV", "NR")

#Import cgpVAF Pindel output matrix if filepath is specified (optional)
if(!is.null(filepath_INDEL)) {
  cat("Importing INDEL matrix\n")
  mat_INDEL <- read.delim(filepath_INDEL, stringsAsFactors = FALSE)
  mat_INDEL <- mat_INDEL[, !grepl("PDv37is", colnames(mat_INDEL))] #Remove reference genome columns
  #filter the low coverage samples if cut-off is specified
  if(min_sample_mean_cov > 0|!is.null(exclude_samples)) {
    remove_cols <- which(str_split(colnames(mat_INDEL), pattern = "_", simplify = TRUE)[,1] %in% c(low_cov_names,exclude_samples))
    mat_INDEL <- mat_INDEL[,-remove_cols]
  }
  mat_INDEL$mut_ref <- paste(mat_INDEL$Chrom, mat_INDEL$Pos, mat_INDEL$Ref, mat_INDEL$Alt, sep = "-") #Create mut_ref column
  mat_INDEL <- mat_INDEL[,c(ncol(mat_INDEL), 1:(ncol(mat_INDEL)-1))] #Swap mut_ref column to the beginning of matrix
  mat_INDEL <- mat_INDEL[,!colnames(mat_INDEL)== "Type"]
  mat_INDEL$Mut_type = "INDEL"
  NV_INDEL <- mat_INDEL[,grepl("MTR", colnames(mat_INDEL))]
  NR_INDEL <- mat_INDEL[,grepl("DEP", colnames(mat_INDEL))]
  mat_INDEL <- mat_INDEL[,!grepl("MTR", colnames(mat_INDEL)) & !grepl("DEP", colnames(mat_INDEL))]
  
  INDEL_mats <- list(mat_INDEL, NV_INDEL, NR_INDEL)
  names(INDEL_mats) = c("mat", "NV", "NR")
  
  #Make combined matrices
  mat_COMB = rbind(mat_SNV, mat_INDEL)
  NV_COMB = rbind(NV_SNV, NV_INDEL)
  NR_COMB = rbind(NR_SNV, NR_INDEL)
} else {
  mat_COMB = mat_SNV
  NV_COMB = NV_SNV
  NR_COMB = NR_SNV
}

COMB_mats = list(mat_COMB, NV_COMB, NR_COMB)
names(COMB_mats) = c("mat", "NV", "NR")

#ASSIGN GENDER by looking for Y chromosome in "Chrom" column. Occasional reads may be mismapped to the Y chromosome, so set threshold of 1%
if(sum(COMB_mats$mat$Chrom == "Y")/sum(COMB_mats$mat$Chrom == "X") > 0.01) {
  COMB_mats$gender <- "male"
} else {
  COMB_mats$gender <- "female"
}

##SET FILTER PARAMETERS - just need to have "min_variant_reads" to be included in analysis for first section (to remove null mutations)
if(COMB_mats$gender == "male") {
  min_variant_reads_auto = 3
  min_variant_reads_xy = 2
} else if(COMB_mats$gender == "female") {
  min_variant_reads_auto = min_variant_reads_xy = 3
}

#Extract number of samples in the set (i.e. number of colonies included in analysis) from number of columns in NV matrix
COMB_mats$nsamp = ncol(COMB_mats$NV)

#Remove mutations with no samples meeting variant read threshold (was likely in bedfile due to excluded sample)
null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto|(COMB_mats$NV >= min_variant_reads_xy & COMB_mats$mat$Chrom %in% c("X","Y"))) == 0
COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
print(paste(sum(null_remove),"mutations removed as no sample meets minimum threshold of number of variant reads"))

#Now apply the filters.  If you are removing mixed colonies, run the germline & binomial first, use these to recognize and remove mixed colonies before running the other filters. Otherwise, can apply all the filters in one go.
if(opt$m) {
	#Initial filters used to recognise and remove mixed colonies
  fun.list = c(germline.binomial.filter,
               beta.binom.filter)
   #Subsequent filters
  fun.list.2 = c(pval_matrix,
                 get_mean_depth,
                 low_vaf_in_pos_samples_dp2,
                 low_vaf_in_pos_samples_dp3,
                 get_max_depth_in_pos)
  
  #Run fun.list functions in parallel
  fun.out = mclapply(fun.list, function(f) {f(COMB_mats)},mc.preschedule = F)
  
  #Use these parameters to recognise and remove mixed colony samples
  basic_params <- as.data.frame(Reduce(cbind,fun.out))
  colnames(basic_params) <- c("germline_pval","bb_rhoval")
  
  #Now can use the results of these to determine any mixed colonies
  colnames(COMB_mats$NV) <- gsub(pattern = "_MTR", replacement = "",x = colnames(COMB_mats$NV))
  sample_peak_vaf = sapply(colnames(COMB_mats$NV), check_peak_vaf, COMB_mats=COMB_mats, filter_params = basic_params)
  mixed_samples = names(sample_peak_vaf[sample_peak_vaf <0.43])
  pass_samples = names(sample_peak_vaf[sample_peak_vaf >=0.43])
  print(paste("Removing sample", mixed_samples,"as peak VAF is", round(sample_peak_vaf[sample_peak_vaf <0.43], digits = 2), "suggesting a non-clonal sample"))
  
  #Save the VAF plots for these discarded samples
  pdf(file = paste0(output_directory, "/mixed_colony_vaf_plots_",Run_ID,".pdf"))
  sapply(mixed_samples, vaf_density_plot, COMB_mats, filter_params = basic_params)
  dev.off()
  pdf(file = paste0(output_directory, "/pass_colony_vaf_plots_",Run_ID,".pdf"))
  sapply(pass_samples, vaf_density_plot, COMB_mats, filter_params = basic_params)
  dev.off() 
  
  #Now remove these sample columns from the matrix before further filter parameters are calculated
  COMB_mats$NR <- COMB_mats$NR[,-which(colnames(COMB_mats$NV) %in% mixed_samples)]
  COMB_mats$NV <- COMB_mats$NV[,-which(colnames(COMB_mats$NV) %in% mixed_samples)]
  null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto|(COMB_mats$NV >= min_variant_reads_xy & COMB_mats$mat$Chrom %in% c("X","Y"))) == 0
  COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
  basic_params = basic_params[!null_remove,]
  print(paste(sum(null_remove),"mutations removed as no sample meets minimum threshold of number of variant reads now that mixed colonies have been removed"))
  
  fun.out.2 = mclapply(fun.list.2, function(f) {f(COMB_mats)})
  
  #Assign output of the first function to the main set of matrices - this is the "PVal matrix" for the final filter
  #Then delete it from the function output
  COMB_mats$PVal <- fun.out.2[[1]]
  fun.out.2[[1]] <- NULL
  
  #The remaining output of the function is the parameters list for all except the "max_pval_in_pos" filters
  basic_params.2 <- as.data.frame(Reduce(cbind,fun.out.2))
  filter_params <-cbind(basic_params, basic_params.2)
  
} else {
  #Single list of functions to apply over data, so that these can then be done in parallel (these functions are all sourced from the foetal.filters.parallel.R script)
  fun.list = c(pval_matrix,
               germline.binomial.filter,
               beta.binom.filter,
               get_mean_depth,
               low_vaf_in_pos_samples_dp2,
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
}

colnames(filter_params) <- c("germline_pval","bb_rhoval","mean_depth","pval_within_pos_dp2","pval_within_pos_dp3","max_depth_in_pos_samples")
rownames(filter_params) <- COMB_mats$mat$mut_ref

#Now run the pval in pos filter. This has to be run separately as it relies on the "PVal" element of the COMB_mats matrix set.
filter_params$max_pval_in_pos_sample <- get_max_pval_in_pos(COMB_mats)

#Now run the final filter. Simple "maximum vaf" within positive sample measure.
filter_params$max_mut_vaf <- get_max_vaf(COMB_mats)

#Now save the matrices to the specified output folder
cat("Saving matrices for downstream analysis\n")
save(COMB_mats, filter_params, file = paste0(output_directory,"/mats_and_params_", Run_ID))

#Create reduced set (to make more manageable file sizes), where samples that are highly likely to be SNPs are pre-filtered: this is the very stringent criteria of meeting both binomial & beta-binomial filter criteria to very stringent cut-offs.
reduced = filter_params$germline_pval < 0.01|filter_params$bb_rhoval >0.01
filter_params = filter_params[reduced,]
COMB_mats = list_subset(COMB_mats, reduced)
save(COMB_mats, filter_params, file = paste0(output_directory, "/mats_and_params_", Run_ID,"_reduced"))
