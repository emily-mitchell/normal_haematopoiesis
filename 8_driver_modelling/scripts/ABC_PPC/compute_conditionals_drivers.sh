#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
CONDITIONAL_RESIM_PATH=${3}
POST_SAMPLE_FILE=${4}
OBS_STATS_FILE=${5}
POOLED_CONDITIONAL_RESIM_PATH=${6}
POOLED_CONDITIONAL_RESIM_FILE_NAME=${7}
N_RESIMS_PER_POST_OBS=${8}
PATIENT_ID=${9}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${SOURCE_PATH} ${CONDITIONAL_RESIM_PATH} ${POST_SAMPLE_FILE} ${OBS_STATS_FILE} ${POOLED_CONDITIONAL_RESIM_PATH} ${POOLED_CONDITIONAL_RESIM_FILE_NAME} ${N_RESIMS_PER_POST_OBS} ${PATIENT_ID}"


R_script_name="compute_conditionals_drivers.R"
Rout_filename="compute_conditionals_drivers"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}.Rout"

