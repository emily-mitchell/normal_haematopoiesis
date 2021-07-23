
####################################################################################################
#
# We have a large sample of simulated data sets (from each conditional distribution);
# This is a (stratified) sample from the posterior predictive distribution, for each summary statistic;
#
# Furtheremore, for each summary statistic, we have a sample from the posterior predictive distribution of
# the standardised residual of the specified summary statistic;
# In addition to these "simulated" standardised residuals, 
# we also have "observed" standardised residuals;
#
# This sample is the output from the scripts:
# resim_conditionals_drivers.R;
# compute_conditionals_drivers.R;
#
# This sample is the input for the present script (compute_ppp_value_drivers.R);
#
# The present script (compute_ppp_value_drivers.R) reads in a vector of summary statistic names ("stat_name_vect"), associated with a specific ID ("stat_set_ID");
# For this specified vector of summary statistic names ("stat_name_vect")
# we select the corresponding columns of ("simulated" and "observed") standardised residuals,
# from which we compute the ("simulated" and "observed") chi-squared discrepency variables;
#
# From the sample of ("simulated" and "observed") chi-squared discrepency variables,
# we estimat the posterior probability that the "observed discrepency" value exceeds the "simulated discrepency" value.
#
# This posterior probability is the posterior predictive p-value.
#
# If the p-value is close to 0.5, this is evidence that the observed data is compatible with the posterior distribution
# with respect to all the features of the data specified by the vector of summary statistic names ("stat_name_vect").
#
# If the p-value is close to zero, this is evidence that the observed data is incompatible with the posterior distribution,
# in the sense that the observed data is exceptional with respect to at least one feature of the data specified by the vector of summary statistic names ("stat_name_vect").
#
# Furtheremore, if the p-value is close to zero, this is evidence that the "range of models" (the parameter space) on which we have computed the posterior distribution is too small,
# in the sense that all these fail to explain some features of the data specified by the vector of summary statistic names ("stat_name_vect").
#
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

##SET FILEPATHS
source_directory = toString(args[1])
pooled_resim_file = toString(args[2])
obs_stats_file = toString(args[3])
output_directory = toString(args[4])
n_resims_per_post_obs = as.integer(args[5])
patient_ID = toString(args[6])
stat_set_ID = toString(args[7])
n_stat_set=as.integer(args[8])

begin_index=9
end_index=begin_index+n_stat_set-1
stat_name_vect=args[begin_index:end_index]

cat( "0.0. compute_ppp_value_drivers.R main: n_stat_set = ", n_stat_set, "\n", sep = "\t" )
cat( "0.0. compute_ppp_value_drivers.R main: length( stat_name_vect ) = ", length( stat_name_vect ) , "\n", sep = "\t" )
cat( "0.0. compute_ppp_value_drivers.R main: stat_name_vect = ", stat_name_vect , "\n", sep = "\t" )
# readLines("stdin",n=1)

library(tidyr)
library(dplyr)
# library(seqinr)

####################################################################################################
#
# Read in a large stratified sample from the posterior predictive distribution;
#
# This sample is the output from the scripts:
# resim_conditionals_drivers.R;
# compute_conditionals_drivers.R;
#

stratified_resim_table = read.table( file = pooled_resim_file, sep = "\t", header=TRUE, row.names = NULL, stringsAsFactors = FALSE )

cat( "1.0. compute_ppp_value_drivers.R main: dim( stratified_resim_table ) = ", dim( stratified_resim_table ), "\n", sep = "\t" )
cat( "1.0. compute_ppp_value_drivers.R main: names( stratified_resim_table ) = ", names( stratified_resim_table ), "\n", sep = "\t" )
# readLines("stdin",n=1)

#
# Create vectors of column names:
#
# resim_col_names_vect = c( "post_obs_index", "n_resims", cols_params, hidden_variable_col_names_vect, cols_stats, conditional_mean_col_names_vect, conditional_SD_col_names_vect, sim_residual_col_names_vect, obs_residual_col_names_vect )
#

cols_params= c("number_drivers_per_year","gamma_shape","gamma_rate")
# number_stats=32
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")

cat( "1.0. compute_ppp_value_drivers.R main: length( cols_stats ) = ", length( cols_stats ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	  
n_top_drivers = 5
top_drivers_col_names_vect = sapply( 1:n_top_drivers, function(x){ paste( "driver_fitness_", x, sep="" ) } )

n_top_clades = 5
top_clades_col_names_vect = sapply( 1:n_top_clades, function(x){ paste( "clade_fitness_", x, sep="" ) } )

driver_count_col_names_vect = c( "shared_drivers", "private_drivers", "total_drivers" )

hidden_variable_col_names_vect = c( top_drivers_col_names_vect, top_clades_col_names_vect, driver_count_col_names_vect )


conditional_mean_col_names_vect = sapply( cols_stats, function(x){ paste( "mean_", x, sep="" ) } )
conditional_SD_col_names_vect = sapply( cols_stats, function(x){ paste( "SD_", x, sep="" ) } )

sim_residual_col_names_vect = sapply( cols_stats, function(x){ paste( "sim_std_resid_", x, sep="" ) } )
obs_residual_col_names_vect = sapply( cols_stats, function(x){ paste( "obs_std_resid_", x, sep="" ) } )

####################################################################################################
#
# Read summary statistics computed from observed data from file "obs_stats_file";
#

obs_stats_table = read.table( file = obs_stats_file, sep = "\t", header=TRUE, row.names = NULL, stringsAsFactors = FALSE )

obs_stats_vect = unlist( obs_stats_table[ 1, cols_stats ] )

####################################################################################################
#
# Perform a Bayesian Chi-squared calculation;
#
####################################################################################################
#
# Select summary statitics to be used in the Chi-squared calculation;
# Select observations from the posterior predictive distribution to be used in the Chi-squared calculation;
#

# stats_included_col_names_vect = c(ltt_stats) # These are the summary stats used in the ABC computation;
stats_included_col_names_vect = stat_name_vect

stats_included_index_vect = which( cols_stats %in% stats_included_col_names_vect )
# stats_included_index_vect = 2:length( cols_stats ) # These are ALL the summary stats contained in "obs_stats_table" excluding column 1 ("mean_muts_data");
# stats_included_index_vect = 1:length( cols_stats )

sim_residual_included_col_names_vect = sim_residual_col_names_vect[ stats_included_index_vect ]
obs_residual_included_col_names_vect = obs_residual_col_names_vect[ stats_included_index_vect ]

residual_included_col_names_vect = c( sim_residual_included_col_names_vect, obs_residual_included_col_names_vect )

####################################################################################################
# 
# Create data frame "ppc_sample_table";
# 

ppc_sample_col_names_vect = c( "post_obs_index", "n_resims", "chisq_sim", "chisq_obs", "diff_chisq", residual_included_col_names_vect )
ppc_sample_table = as.data.frame( matrix( NA, nrow=nrow( stratified_resim_table ), ncol=length( ppc_sample_col_names_vect ) ) )
names( ppc_sample_table ) <- ppc_sample_col_names_vect

ppc_sample_table[ , c( "post_obs_index", "n_resims", residual_included_col_names_vect ) ] = stratified_resim_table[ , c( "post_obs_index", "n_resims", residual_included_col_names_vect ) ]

complete_rows_index_vect = complete.cases( ppc_sample_table[ , residual_included_col_names_vect ] )
ppc_sample_table = ppc_sample_table[ complete_rows_index_vect, ]

good_est_index_vect = which( ppc_sample_table[[ "n_resims" ]] >= n_resims_per_post_obs )
ppc_sample_table = ppc_sample_table[ good_est_index_vect, ]

# 
# Compute squares of standardised residuals;
#

# temp_table <- ppc_sample_table %>% mutate_at( residual_included_col_names_vect, function(x){ x^2 } )
temp_table <- ppc_sample_table[ , residual_included_col_names_vect ] %>% mutate_all( function(x){ x^2 } )

#
# Perform a Bayesian Chi-squared calculations;
#
ppc_sample_table[[ "chisq_sim" ]] = rowSums( temp_table[ , sim_residual_included_col_names_vect ] )
ppc_sample_table[[ "chisq_obs" ]]  = rowSums( temp_table[ , obs_residual_included_col_names_vect ] )
ppc_sample_table[[ "diff_chisq" ]] = ppc_sample_table[[ "chisq_sim" ]] - ppc_sample_table[[ "chisq_obs" ]]

####################################################################################################
# 
# Create data frame "posterior_obs_info_table";
# 

posterior_obs_col_vect <- unique( stratified_resim_table[[ "post_obs_index" ]] )

posterior_obs_info_col_names_vect = c( "post_obs_index", "n_sims", "n_complete_sims" )
posterior_obs_info_table = as.data.frame( matrix( NA, nrow=length( posterior_obs_col_vect ), ncol=length( posterior_obs_info_col_names_vect) ) )
names( posterior_obs_info_table ) <- posterior_obs_info_col_names_vect

posterior_obs_info_table[[ "post_obs_index" ]] = posterior_obs_col_vect

# 
# Count number of resims obtained for each "post_obs_index";
#

index_vect_list = lapply( as.list( posterior_obs_col_vect ), function(x){ which( stratified_resim_table[[ "post_obs_index" ]] == x ) } )

sim_count_list = lapply( index_vect_list, length )

complete_sim_count_list = lapply( index_vect_list, function(x){ length( which( complete.cases( stratified_resim_table[ x, stat_name_vect ] ) ) ) } )

posterior_obs_info_table[[ "n_sims" ]] = unlist( sim_count_list )
posterior_obs_info_table[[ "n_complete_sims" ]] = unlist( complete_sim_count_list )

sim_count_index_vect = order( posterior_obs_info_table[[ "n_complete_sims" ]], decreasing = TRUE )
posterior_obs_info_table = posterior_obs_info_table[ sim_count_index_vect, ]

# 
# Write data frame "posterior_obs_info_table";
#

posterior_obs_info_outfile = paste( output_directory, "/", "posterior_obs_info_table", "_stat_set_", stat_set_ID, ".txt", sep = "" )
write.table( posterior_obs_info_table, file = posterior_obs_info_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

####################################################################################################
# 
# Create data frame "single_stat_tests_table";
# 

single_stat_tests_col_names_vect = c( "stat_name", "n_sims", "n_reject", "p_reject", "p_formatted" )
single_stat_tests_table = as.data.frame( matrix( NA, nrow=length( stat_name_vect ), ncol=length( single_stat_tests_col_names_vect) ) )
names( single_stat_tests_table ) <- single_stat_tests_col_names_vect

single_stat_tests_table[[ "stat_name" ]] = stat_name_vect

# 
# Compute squares of standardised residuals;
#

# temp_table <- ppc_sample_table %>% mutate_at( residual_included_col_names_vect, function(x){ x^2 } )
temp_table <- ppc_sample_table[ , residual_included_col_names_vect ] %>% mutate_all( function(x){ x^2 } )

# ppc_sample_table[[ "chisq_sim" ]] = rowSums( temp_table[ , sim_residual_included_col_names_vect ] )
# ppc_sample_table[[ "chisq_obs" ]]  = rowSums( temp_table[ , obs_residual_included_col_names_vect ] )
# ppc_sample_table[[ "diff_chisq" ]] = ppc_sample_table[[ "chisq_sim" ]] - ppc_sample_table[[ "chisq_obs" ]]

diff_table = as.data.frame( matrix( NA, nrow=nrow( ppc_sample_table ), ncol=length( stat_name_vect ) ) )
names( diff_table ) <- stat_name_vect

diff_table[ , stat_name_vect ] = temp_table[ , sim_residual_included_col_names_vect ] - temp_table[ , obs_residual_included_col_names_vect ]


n_sims_finite_vect = sapply( stat_name_vect, function(x){ length( which( is.finite( diff_table[[ x ]] ) ) ) } )

n_reject_vect = sapply( stat_name_vect, function(x){ length( which( diff_table[[ x ]][ is.finite( diff_table[[ x ]] ) ] > 0 ) ) } )

single_stat_tests_table[[ "n_sims" ]] = n_sims_finite_vect
single_stat_tests_table[[ "n_reject" ]] = n_reject_vect
single_stat_tests_table[[ "p_reject" ]] = single_stat_tests_table[[ "n_reject" ]] / single_stat_tests_table[[ "n_sims" ]]
single_stat_tests_table[[ "p_formatted" ]] = round( single_stat_tests_table[[ "p_reject" ]], digits = 4 )

# 
# Write data frame "single_stat_tests_table";
#

single_stat_tests_outfile = paste( output_directory, "/", "single_stat_tests_table", "_stat_set_", stat_set_ID, ".txt", sep = "" )
write.table( single_stat_tests_table, file = single_stat_tests_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

####################################################################################################
# 
# Create data frame "ppc_result_table";
# 

ppc_result_col_names_vect = c( "n_post_obs_represented", "min_resims_per_post_obs", "n_sims", "n_reject", "p_reject", "p_formatted", stats_included_col_names_vect )
ppc_result_table = as.data.frame( matrix( NA, nrow=1, ncol=length( ppc_result_col_names_vect) ) )
names( ppc_result_table ) <- ppc_result_col_names_vect

ppc_result_table[ , stats_included_col_names_vect ] = obs_stats_table[ , stats_included_col_names_vect ]

complete_rows_row_index_vect = which( is.finite( ppc_sample_table[[ "chisq_obs" ]] ) )
ppc_sample_table = ppc_sample_table[ complete_rows_index_vect, ]

n_post_obs_good_est = length( unique( ppc_sample_table[[ "post_obs_index" ]] ) )

#Estimate (Bayesian) p-value:
n_sims = nrow( ppc_sample_table )
n_reject = length( which( ppc_sample_table[[ "diff_chisq" ]] > 0 ) )
p_reject = n_reject / n_sims
p_formatted = round( p_reject, digits = 4 )

ppc_result_table[[ "n_post_obs_represented" ]][ 1 ] = n_post_obs_good_est
ppc_result_table[[ "min_resims_per_post_obs" ]][ 1 ] = n_resims_per_post_obs
ppc_result_table[[ "n_sims" ]][ 1 ] = n_sims
ppc_result_table[[ "n_reject" ]][ 1 ] = n_reject
ppc_result_table[[ "p_reject" ]][ 1 ] = p_reject
ppc_result_table[[ "p_formatted" ]][ 1 ] = p_formatted

# 
# Write data frame "ppc_result_table";
#

# pp_test_outfile = paste( output_directory, "/", "test_result_table", ".txt", sep = "" )
ppc_result_outfile = paste( output_directory, "/", "ppc_result_table", "_stat_set_", stat_set_ID, ".txt", sep = "" )
write.table( ppc_result_table, file = ppc_result_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

####################################################################################################
# 
# Write data frame "ppc_sample_table" to outfile;
#

ppc_sample_outfile = paste( output_directory, "/", "ppc_sample_table", "_stat_set_", stat_set_ID, ".txt", sep = "" )
write.table( ppc_sample_table, file = ppc_sample_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

####################################################################################################
#
# Plot distributions;
#
setwd( output_directory )
system("mkdir pdfs")

#Plot posterior distribution of Chisq (= chisq_obs);
plot_title="Distribution of chisq_obs"
# pdf(paste0("pdfs/",ID, "_chisq_obs.pdf"))
pdf(paste0("pdfs/", "patient_", patient_ID, "_stat_set_", stat_set_ID, "_chisq_obs.pdf"))
hist( ppc_sample_table$chisq_obs, main = plot_title, xlab = "chisq_obs", breaks = 20)
dev.off()

#Plot posterior predictive distribution of Chisq (= chisq_sim);
plot_title="Distribution of chisq_sim"
# pdf(paste0("pdfs/",ID, "_chisq_sim.pdf"))
pdf(paste0("pdfs/", "patient_", patient_ID, "_stat_set_", stat_set_ID, "_chisq_sim.pdf"))
hist( ppc_sample_table$chisq_sim, main = plot_title, xlab = "chisq_sim", breaks = 20)
dev.off()

#Plot posterior distribution of diff_chisq (= chisq_sim - chisq_obs);
plot_title=paste0("Distribution of chisq_sim - chisq_obs (P value = ", p_formatted, ")")
# pdf(paste0("pdfs/",ID, "_diff_chisq.pdf"))
pdf(paste0("pdfs/", "patient_", patient_ID, "_stat_set_", stat_set_ID, "_diff_chisq.pdf"))
hist( ppc_sample_table$diff_chisq, main = plot_title, xlab = "chisq_sim - chisq_obs", breaks = 20)
abline(v = 0,col="red")
dev.off()

####################################################################################################

cat( "done!", "\n" )


