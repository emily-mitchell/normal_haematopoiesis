## Run on farm5
args = commandArgs(TRUE)

##SET param_ID AND FILEPATHS
# RUNID = as.character(args[1])
input_file_path = toString(args[1])
resim_directory = toString(args[2])
output_directory = toString(args[3])
source_directory = toString(args[4])

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

#Set working directory
ID = "KX003"
Iteration = "KX003_5"
stat_time=c(0.2,0.4,0.6,0.8)
my_working_directory = output_directory

setwd(my_working_directory)

load(input_file_path)

pp_param_table <- as.data.frame(KX003.parameters.ppc)
pp_summary_stats_table <- as.data.frame(KX003.summary_stats.ppc)
obs_summary_stats_table <- as.data.frame(KX003.summary_stats.data)

cat( "0.0. pp_chi_squared_test: names( pp_param_table ) = ", names( pp_param_table ), "\n", sep = "\t" )
cat( "0.0. pp_chi_squared_test: names( obs_summary_stats_table ) = ", names( obs_summary_stats_table ), "\n", sep = "\t" )

cols_params=colnames(pp_param_table)
cols_stats=colnames(obs_summary_stats_table)

names(pp_summary_stats_table) <- names(obs_summary_stats_table)

##Set up data frame "pp_test_table"
cols_stats_conditional_means=sapply( cols_stats , function(x){ paste0( "mean_conditional_", x ) } )
cols_stats_conditional_vars=sapply( cols_stats , function(x){ paste0( "var_conditional_", x ) } )
cols_pp_test=c( "pp_sim_index", cols_params, cols_stats, cols_stats_conditional_means, cols_stats_conditional_vars, "n_resims", "chisq_sim", "chisq_obs", "diff_chisq" )
# pp_test_table=matrix(NA,nrow=nrow(pp_param_table),ncol=length(cols_pp_test))
pp_test_table=as.data.frame( matrix(NA,nrow=nrow(pp_param_table),ncol=length(cols_pp_test)) )
names(pp_test_table) <- cols_pp_test

##Fill data frame "pp_test_table"
pp_test_table[[ "pp_sim_index" ]] = 1:nrow(pp_param_table)
pp_test_table[[ "n_resims" ]] = rep( 0, nrow(pp_param_table) )
pp_test_table[ , cols_params ] = pp_param_table
pp_test_table[ , cols_stats ] = pp_summary_stats_table

##Estimate moments of conditional distribution
for(pp_sim_index in 1:nrow(pp_param_table)) {
  
cat( "1.0. pp_chi_squared_test: pp_sim_index = ", pp_sim_index, "\n", sep = "\t" )
# readLines("stdin",n=1)
  
  #Define list of simulation tree and parameter files
  stats_resim_pattern=paste0("summary_stats_pp_sim_", pp_sim_index, "_")
  stats_resim_files=list.files(resim_directory,pattern=stats_resim_pattern)
  
  stats_resim_files_vect = sapply( stats_resim_files, function(x){ paste0( resim_directory, x ) } )
  stats_resim_table_list = list()
  if( length( stats_resim_files ) > 0 ) {
	  for( i in 1:length( stats_resim_files ) ) {
		get_result = try( stats_resim_table <- get( load( stats_resim_files_vect[ i ] ) ) )
		if( class( get_result ) != "try-error" ) {
			stats_resim_table_list = append( stats_resim_table_list, list( stats_resim_table ) )
		}
	  }
  }
  
cat( "1.1. pp_chi_squared_test: length( stats_resim_table_list ) = ", length( stats_resim_table_list ), "\n", sep = "\t" )

  if( length( stats_resim_table_list ) > 0 ) {

	  if( length( stats_resim_table_list ) == 1 ) {
		  summary_stats.resim == stats_resim_table_list[[ 1 ]]
	  }else{
		  summary_stats.resim <- do.call("rbind", stats_resim_table_list)
	  }
	  
	  summary_stats.resim = as.data.frame( summary_stats.resim )
	  
cat( "1.2. pp_chi_squared_test: dim( summary_stats.resim ) = ", dim( summary_stats.resim ), "\n", sep = "\t" )
# readLines("stdin",n=1)
  
	  if( nrow( summary_stats.resim ) > 0 )
	  {
		pp_test_table[[ "n_resims" ]][ pp_sim_index ] = nrow( summary_stats.resim )
		
		conditional_mean_vect = sapply( summary_stats.resim, mean )
		conditional_var_vect = sapply( summary_stats.resim, var )
		
		pp_test_table[ pp_sim_index, cols_stats_conditional_means ] = conditional_mean_vect
		pp_test_table[ pp_sim_index, cols_stats_conditional_vars ] = conditional_var_vect
		
	  }
  } # if( length( stats_resim_files ) > 0 )
  
} # for(pp_sim_index in 1:nrow(pp_param_table))

cat( "1.1. pp_chi_squared_test: dim( pp_param_table ) = ", dim( pp_param_table ), "\n", sep = "\t" )
cat( "1.1. pp_chi_squared_test: names( pp_param_table ) = ", names( pp_param_table ), "\n", sep = "\t" )
print( head( pp_param_table ) )

cat( "1.2. pp_chi_squared_test: dim( pp_summary_stats_table ) = ", dim( pp_summary_stats_table ), "\n", sep = "\t" )
cat( "1.2. pp_chi_squared_test: names( pp_summary_stats_table ) = ", names( pp_summary_stats_table ), "\n", sep = "\t" )
print( head( pp_summary_stats_table ) )

cat( "1.3. pp_chi_squared_test: dim( obs_summary_stats_table ) = ", dim( obs_summary_stats_table ), "\n", sep = "\t" )
cat( "1.3. pp_chi_squared_test: names( obs_summary_stats_table ) = ", names( obs_summary_stats_table ), "\n", sep = "\t" )
print( head( obs_summary_stats_table ) )

complete_rows_index_vect = complete.cases( pp_test_table[ , c( cols_params, cols_stats, cols_stats_conditional_means, cols_stats_conditional_vars ) ] )
pp_test_table = pp_test_table[ complete_rows_index_vect, ]

cat( "1.4. pp_chi_squared_test: dim( pp_test_table ) = ", dim( pp_test_table ), "\n", sep = "\t" )
cat( "1.4. pp_chi_squared_test: names( pp_test_table ) = ", names( pp_test_table ), "\n", sep = "\t" )
print( head( pp_test_table ) )

##Fill data frame "pp_test_table" columns: "chisq_sim", "chisq_obs", "diff_chisq";

ltt_stats = grep("ltt",colnames(pp_summary_stats_table))
coals_stats = grep("coals",colnames(pp_summary_stats_table))
clade_stats = grep("clade",colnames(pp_summary_stats_table))
singleton_stats = grep("singletons",colnames(pp_summary_stats_table))

sumstats_to_include = c(ltt_stats) # These are the summary stats used in the ABC computation;

sumstats_to_include_index_vect = 2:length(cols_stats)

included_cols_stats = cols_stats[ sumstats_to_include_index_vect ]
included_cols_stats_conditional_means = cols_stats_conditional_means[ sumstats_to_include_index_vect ]
included_cols_stats_conditional_vars = cols_stats_conditional_vars[ sumstats_to_include_index_vect ]

pp_summary_stats_vect_list = as.list( as.data.frame( t( as.matrix( pp_test_table[ , included_cols_stats ] ) ) ) )
pp_summary_stats_conditional_means_vect_list = as.list( as.data.frame( t( as.matrix( pp_test_table[ , included_cols_stats_conditional_means ] ) ) ) )
pp_summary_stats_conditional_vars_vect_list = as.list( as.data.frame( t( as.matrix( pp_test_table[ , included_cols_stats_conditional_vars ] ) ) ) )

obs_summary_stats_vect = c( as.matrix( obs_summary_stats_table[ , included_cols_stats ] ) )

####################################################################################################

calc_chisq_stat <- function(y,m,v)
{
 	
	vect = (y-m)*(y-m) / v
	
	stat = sum( vect )
	
	return( stat )
	
}

####################################################################################################

chisq_sim_list = Map( function(y,m,v){ calc_chisq_stat(y,m,v) }, pp_summary_stats_vect_list, pp_summary_stats_conditional_means_vect_list, pp_summary_stats_conditional_vars_vect_list )
chisq_sim_vect = unlist( chisq_sim_list )

cat( "3.1. pp_chi_squared_test: chisq_sim_vect = ", chisq_sim_vect, "\n", sep = "\t" )

chisq_obs_list = Map( function(m,v){ calc_chisq_stat( obs_summary_stats_vect, m, v ) }, pp_summary_stats_conditional_means_vect_list, pp_summary_stats_conditional_vars_vect_list )
chisq_obs_vect = unlist( chisq_obs_list )

cat( "3.1. pp_chi_squared_test: chisq_obs_vect = ", chisq_obs_vect, "\n", sep = "\t" )

pp_test_table[[ "chisq_sim" ]] = chisq_sim_vect
pp_test_table[[ "chisq_obs" ]] = chisq_obs_vect
pp_test_table[[ "diff_chisq" ]] = pp_test_table[[ "chisq_sim" ]] - pp_test_table[[ "chisq_obs" ]]

test_table_outfile = paste( output_directory, "/", "test_table", ".txt", sep = "" )
write.table( pp_test_table, file = test_table_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )


complete_rows_row_index_vect = which( is.numeric( pp_test_table[[ "chisq_obs" ]] ) )
pp_test_table = pp_test_table[ complete_rows_index_vect, ]

#Estimate (Bayesian) p-value:
n_pp_sims = nrow( pp_test_table )
n_reject = length( which( pp_test_table[[ "diff_chisq" ]] > 0 ) )
p_reject = n_reject / n_pp_sims
p_formatted = round( p_reject, digits = 4 )

pp_test_result_table=data.frame(p_reject=p_reject, n_sims=n_pp_sims, n_reject=n_reject)

pp_test_result_table_outfile = paste( output_directory, "/", "test_result_table", ".txt", sep = "" )
write.table( pp_test_result_table, file = pp_test_result_table_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

#Assign the key objects to objects incorporating the model name
assign(paste0(ID,".parameters.ppc"),pp_param_table)
assign(paste0(ID,".summary_stats.ppc"),pp_summary_stats_table)
assign(paste0(ID,".test_table.ppc"),pp_test_table)
assign(paste0(ID,".summary_stats.data"),obs_summary_stats_table)

save("KX003.parameters.ppc",
     "KX003.summary_stats.ppc",
     "KX003.test_table.ppc",
     "KX003.summary_stats.data",
     file = paste0(ID,"_tables"))

assign(paste0(ID,".result_table.data"),pp_test_result_table)

save("KX003.test_table.ppc",
     "KX003.result_table.data",
     file = paste0(ID,"_result"))
	

#Plot distributions;
system("mkdir pdfs")

#Plot posterior distribution of Chisq (= chisq_obs);
plot_title="Distribution of chisq_obs"
pdf(paste0("pdfs/",ID, "_chisq_obs.pdf"))
hist( pp_test_table$chisq_obs, main = plot_title, xlab = "chisq_obs", breaks = 20)
dev.off()

#Plot posterior predictive distribution of Chisq (= chisq_sim);
plot_title="Distribution of chisq_sim"
pdf(paste0("pdfs/",ID, "_chisq_sim.pdf"))
hist( pp_test_table$chisq_sim, main = plot_title, xlab = "chisq_sim", breaks = 20)
dev.off()

#Plot posterior distribution of diff_chisq (= chisq_sim - chisq_obs);
plot_title=paste0("Distribution of chisq_sim - chisq_obs (P value = ", p_formatted, ")")
pdf(paste0("pdfs/",ID, "_diff_chisq.pdf"))
hist( pp_test_table$diff_chisq, main = plot_title, xlab = "chisq_sim - chisq_obs", breaks = 20)
abline(v = 0,col="red")
dev.off()
