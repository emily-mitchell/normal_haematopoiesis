#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
PP_SIMS_FILE=${3}
RESIM_PATH=${4}
OUTPUT_PATH=${5}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${PP_SIMS_FILE} ${RESIM_PATH} ${OUTPUT_PATH} ${SOURCE_PATH}"

R_script_name="KX003_pp_chi_squared_test.R"
Rout_filename="KX003_pp_chi_squared_test"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
#export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
# Rscript "--no-save --no-restore --verbose ./KX003_pp_chi_squared_test.R ${R_PARAMS} > KX003_pp_chi_squared_test.Rout"
# Rscript "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}.Rout"
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}.Rout"
