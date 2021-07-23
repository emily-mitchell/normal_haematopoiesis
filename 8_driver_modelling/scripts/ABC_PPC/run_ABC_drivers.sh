#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
PATIENT_ID=${3}
TREE_FILE=${4}
OBS_STATS_PATH=${5}
OBS_STATS_FILE_NAME=${6}
POOLED_SIM_FILE=${7}
ABC_OUTPUT_PATH=${8}
POST_SAMPLE_FILE_NAME=${9}
N_SIMS_ACCEPT=${10}
N_SIMS_MAX=${11}
ABC_STAT_SET_ID=${12}
LENGTH_ABC_STAT_NAME_VECT=${13}
# ABC_STAT_NAME_VECT=${14}

ARG_ARRAY=( "$@" ) # ARG_ARRAY is indexed from 0! (while $ is indexed from 1!)

ABC_STAT_NAME_VECT=${ARG_ARRAY[@]:13:${LENGTH_ABC_STAT_NAME_VECT}} # ARG_ARRAY is indexed from 0!

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

# RUNID='KX004'

R_PARAMS="${SOURCE_PATH} ${TREE_FILE} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${POOLED_SIM_FILE} ${ABC_OUTPUT_PATH} ${POST_SAMPLE_FILE_NAME} ${N_SIMS_ACCEPT} ${N_SIMS_MAX} ${PATIENT_ID} ${ABC_STAT_SET_ID} ${LENGTH_ABC_STAT_NAME_VECT} ${ABC_STAT_NAME_VECT[@]}"

echo "R_PARAMS = ${R_PARAMS}"


R_script_name="run_ABC_drivers.R"
Rout_filename="run_ABC_drivers"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_PATIENT_${PATIENT_ID}_JOB_${JOB_ID}.Rout"

