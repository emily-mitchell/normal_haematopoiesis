#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
PATIENT_ID=${3}
TREE_FILE=${4}
OBS_STATS_PATH=${5}
OBS_STATS_FILE_NAME=${6}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${SOURCE_PATH} ${TREE_FILE} ${OBS_STATS_PATH} ${OBS_STATS_FILE_NAME} ${PATIENT_ID}"

echo "R_PARAMS = ${R_PARAMS}"


R_script_name="calc_obs_stats_test_drivers.R"
Rout_filename="calc_obs_stats_test_drivers"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_PATIENT_${PATIENT_ID}_JOB_${JOB_ID}.Rout"

