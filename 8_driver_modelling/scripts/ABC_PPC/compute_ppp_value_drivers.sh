#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
POOLED_CONDITIONAL_RESIM_FILE=${3}
OBS_STATS_FILE=${4}
PPC_OUTPUT_PATH=${5}
N_RESIMS_PER_POST_OBS=${6}
PATIENT_ID=${7}
PPC_STAT_SET_ID=${8}
LENGTH_PPC_STAT_NAME_VECT=${9}
# PPC_STAT_NAME_VECT=${10}

ARG_ARRAY=( "$@" ) # ARG_ARRAY is indexed from 0! (while $ is indexed from 1!)

PPC_STAT_NAME_VECT=${ARG_ARRAY[@]:9:${LENGTH_PPC_STAT_NAME_VECT}} # ARG_ARRAY is indexed from 0!

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${SOURCE_PATH} ${POOLED_CONDITIONAL_RESIM_FILE} ${OBS_STATS_FILE} ${PPC_OUTPUT_PATH} ${N_RESIMS_PER_POST_OBS} ${PATIENT_ID} ${PPC_STAT_SET_ID} ${LENGTH_PPC_STAT_NAME_VECT} ${PPC_STAT_NAME_VECT[@]}"

R_script_name="compute_ppp_value_drivers.R"
Rout_filename="compute_ppp_value_drivers"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}.Rout"








