
####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC sampling.
# 
# At the start of "round 1", we generate a large sample of pairs (parameter vector, simulated data set)
# by sampling parameter vectors from the "prior" distribution,
# and for each parameter vector we run a simulation of the model;
# by the running script: sim_prior_drivers.R (which uses the R package rsimpop);
#
# This sample (the output from the script: sim_prior_drivers.R)
# is the input for the present script (collect_prior_sims_drivers.R);
#
# In fact, the script: sim_prior_drivers.R;
# is run many times within an "job array" on the computing cluster (farm5),
# so that may separate output files are generated;
# The present script (collect_prior_sims_drivers.R) collects and pools together all these output files;
#
# This pooled sample is then the input for the script: run_ABC_drivers.R;
#
# We can also run the present script (sim_prior_drivers.R)
# in order to count the number of completed simulations.
# 
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

source_directory = toString(args[1])
sim_directory = toString(args[2])
pooled_sim_directory = toString(args[3])
pooled_sim_file_name = toString(args[4])
n_sims_per_job = as.integer(args[5])
n_jobs = as.integer(args[6])


sim_file_list = list.files( path = sim_directory, pattern = "sim_" )

setwd( sim_directory )
sim_table_list = list()
for( i in 1:length( sim_file_list ) )
{
	sim_file = sim_file_list[[ i ]]
	
	if( file_test( "-f", sim_file ) )
	{
		sim_result = try( sim_table <- read.table( file = sim_file, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE ) )
		if( class( sim_result ) != "try-error" )
		{
			sim_table_list = append( sim_table_list, list( sim_table ) )
			
		}
		
	}
	
}

if( length( sim_table_list ) > 0 ) 
{
	if( length( sim_table_list ) == 1 ) 
	{
		pooled_sim_table = sim_table_list[[ 1 ]]
		
	}else
	{
		pooled_sim_table <- do.call( "rbind", sim_table_list )
		
	}
	
} # if( length( sim_files ) > 0 )


pooled_sim_outfile = paste( pooled_sim_directory, "/", pooled_sim_file_name, sep = "" )
write.table( pooled_sim_table, file = pooled_sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

###############################################################################################
#
# Count number of sims collected;
# 
# Write data outfile "sim_count_file.txt";
#

n_complete_sims = length( which( complete.cases( pooled_sim_table ) ) )

sim_count_outfile = paste( pooled_sim_directory, "/", "sim_count_file", ".txt", sep = "" )
sink( file = sim_count_outfile )
cat( "number of output files = ", length( sim_table_list ), "\n", sep = "\t" )
cat("\n")
cat( "number of sims collected = ", nrow( pooled_sim_table ), "\n", sep = "\t" )
cat("\n")
cat( "number of sims completed = ", n_complete_sims, "\n", sep = "\t" )
sink()

###############################################################################################

cat( "done!", "\n" )

