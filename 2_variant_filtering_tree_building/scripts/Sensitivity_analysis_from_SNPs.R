#!/software/R-3.6.1/bin/Rscript

# Author: Mike Spencer Chapman

library(optparse)
library(data.table)

##This script uses the sensitivity of calling germline variants (SNVs or INDELs) to estimate overall sensitivity for calling variants of these types
##It has 2 main steps:
#1. Define a list of germline SNVs and INDELs in the format Chrom-Pos-Ref-Alt.  Can do this from the output of the "Mutation_filtering_get_parameters.R" script.  Germline variants
# should have both a relatively high "germline filter" p-value (defined here as > 10^-6), and a low "beta-binomial filter" rho value (defined here as <0.05).
#2. Look through each of the "pass" mutations for each individual sample (i.e. the caveman "pass" SNVs that have also been through the low-input pipeline filtering OR
# the "pass" pindel calls) to see what proportion of germline variants have been called in each sample.
#This is outputed as a sensitivity dataframe where each sample is a row & there are separate columns for SNV and INDEL sensitivity

option_list = list(
  make_option(c("-m", "--mats"), action="store", type='character', help="Path for the mats_and_params file"),
  make_option(c("-o", "--output"), action="store", type='character', help="Directory for output_file"),
  make_option(c("-n", "--name"), action="store", type='character', help="Name for output_file"),
  make_option(c("-i", "--indels_folder"), action="store", default="pindel_raw", type='character', help="Directory for indel files - must have the extension '.pindel.annot.vcf_pass_flags'"),
  make_option(c("-s", "--snv_folder"), action="store", default="MS_filters/output_files", type='character', help="Directory SNV files - must have the extension '_complete_final_retained_3.vcf'")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
print(opt)

#To run this script, must have access to the final sets of filtered mutation calls from running the sample against an unmatched reference genome.
print("Defining objects")
mats_and_params_file = opt$m
output_dir=opt$o
output_file_name=opt$n
SNV_calls_folder=opt$s
INDEL_calls_folder =opt$i

print(paste0("Mats and params file path is: ",mats_and_params_file))
print(paste0("Using output directory: ",output_dir))
print(paste0("Using output file name: ",output_file_name))
print(paste0("Using SNV calls folder: ",SNV_calls_folder))
print(paste0("Using INDEL calls folder: ",INDEL_calls_folder))

##Define the REFERENCE GERMLINE MUTATIONS set for each foetus
#Get reference set of SNPs/ INDELs for the 18 week foetus.  Get from the post-LCM-artefact-filtering mutation sets, but must not be the "reduced" versions (i.e. SNPs still retained)
print("Loading up mats_and_params file")
load(mats_and_params_file) #Path to full output file of the "Mutation_filtering_get_parameters.R" script

print("Define set of germline SNVs and INDELs")
germline_SNVs = COMB_mats$mat$mut_ref[COMB_mats$mat$Mut_type == "SNV" & log10(filter_params$germline_pval) > -6 & filter_params$bb_rhoval < 0.05 & rowSums(COMB_mats$NR) > 500]
germline_INDELs = COMB_mats$mat$mut_ref[COMB_mats$mat$Mut_type == "INDEL" & log10(filter_params$germline_pval) > -6 & filter_params$bb_rhoval < 0.05  & rowSums(COMB_mats$NR) > 500]

#Function to calculate sensitivity from a vcf file of passed mutations for an individual sample, and a reference list of germline mut_refs
get_proportion_of_germline_detected = function(file, reference_germline_muts, skip = "#CHROM") { #By default will skip vcf header. If no header, set skip to 0.
  require(data.table)
  if(file.info(file)$size == 0) {
    pick_up_rate <- 0
  } else {
    pass_vcf <- fread(file, skip = skip, header = FALSE) #Import the vcf file
    pass_vcf[,2] <- sapply(pass_vcf[,2], function(x) gsub("^\\s+|\\s+$","",x)) #strip the white space from pos column
    mut_refs = apply(pass_vcf[,c(1,2,4,5)],1,paste,collapse = "-") #Convert the vcf into vector mutations in format "Chrom-Pos-Ref-Alt". Note that columns must be in the correct vcf order
    pick_up_rate = sum(reference_germline_muts %in% mut_refs)/length(reference_germline_muts) #Calculate what proportion of germline SNVs have been called
  }
  names(pick_up_rate) = file
  print(paste(file," completed"))
  return(pick_up_rate)
}

#Function to create sensitivity dataframe for INDELs & SNVs for each sample
#NB. May need to alter the "skip" argument depending on the exact format of the individual sample vcf files (i.e. if have a header or not)

create_sensitivity_df=function(germline_SNV_set, #A vector of germline SNVs in the format "Chrom-Pos-Ref-Alt" 
                               SNV_calls_folder, #The folder where all the vcf files of called SNVs are for each sample 
                               SNV_calls_file_ext="_complete_final_retained_2.vcf", #The extension for the SNV vcf files (the extension will be removed to get sample names)
                               germline_INDEL_set=NULL, #A vector of germline INDELs in the format "Chrom-Pos-Ref-Alt" 
                               INDEL_calls_folder=NULL, #The folder where all the vcf files of called INDELs are for each sample 
                               INDEL_calls_file_ext=".pindel.annot.vcf_pass_flags") { #The extension for the INDEL vcf files (the extension will be removed to get sample names - these should then match the SNV sample names)
  
  #Get sensitivity for germline SNVs for each sample
  pass_SNV_files = list.files(SNV_calls_folder,pattern = paste0("\\",SNV_calls_file_ext,"$"),full.names=TRUE) #get list of all the files with the appropriate extension in the specified folder
  SNV_sensitivity <- sapply(pass_SNV_files, get_proportion_of_germline_detected, reference_germline_muts = germline_SNV_set) #iterate through these files to get the sensitivity
  names(SNV_sensitivity) <- gsub(pattern=paste0(SNV_calls_folder,"/"),replacement="",x=gsub(pattern = SNV_calls_file_ext,replacement =  "", x =pass_SNV_files)) #name the resulting vector with sample names
  
  #Get sensitivity for germline INDELs for each sample
  if(!is.null(germline_INDEL_set)) {
    pass_INDEL_files = list.files(INDEL_calls_folder,pattern = paste0("\\",INDEL_calls_file_ext,"$"),full.names=TRUE) #get list of all the files with the appropriate extension in the specified folder
    INDEL_sensitivity <- sapply(pass_INDEL_files, get_proportion_of_germline_detected, skip = 0, reference_germline_muts = germline_INDEL_set) #iterate through these files to get the sensitivity
    names(INDEL_sensitivity) <- gsub(pattern=paste0(INDEL_calls_folder,"/"),replacement="",gsub(pattern = INDEL_calls_file_ext,replacement =  "", x = pass_INDEL_files)) #name the resulting vector with sample names
    
    #Collate into data frame
    results_df=data.frame(Sample = names(SNV_sensitivity), SNV_sensitivity = SNV_sensitivity, INDEL_sensitivity = INDEL_sensitivity)
    rownames(results_df)<-NULL #no need to keep rownames as sample name is in the Sample column
    return(results_df)
  } else {
    results_df=data.frame(Sample = names(SNV_sensitivity), SNV_sensitivity = SNV_sensitivity)
    rownames(results_df)<-NULL #no need to keep rownames as sample name is in the Sample column
    return(results_df)
  }
}

##RUN THESE FUNCTIONS TO GET SENSITIVITY DATAFRAMES FOR BOTH SNVs and INDELs
print("Run the function to create the sensitivity data frame")
sensitivity_analysis=create_sensitivity_df(germline_SNV_set=germline_SNVs,
                                           SNV_calls_folder=SNV_calls_folder,
                                           SNV_calls_file_ext="_complete_final_retained_3.vcf",
                                           germline_INDEL_set=germline_INDELs,
                                           INDEL_calls_folder = INDEL_calls_folder)


##SAVE THE DATAFRAMES
print("Saving the dataframe")
write.table(sensitivity_analysis, file = paste0(output_dir,"/",output_file_name), row.names = FALSE)