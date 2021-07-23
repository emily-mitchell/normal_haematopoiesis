#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
OUTPUT_PATH=${3}
POST_SAMPLE_FILE=${4}
CONDITIONAL_RESIM_DIR_NAME=${5}
BEGIN_POST_OBS=${6}
END_POST_OBS=${7}
N_RUNS_PER_POST_OBS=${8}
N_RESIMS_PER_RUN=${9}
N_PATIENTS=${10}
# PATIENT_ARRAY=${11}

ARG_ARRAY=( "$@" ) # ARG_ARRAY is indexed from 0! (while $ is indexed from 1!)

BEGIN_INDEX=10 # ARG_ARRAY is indexed from 0!
PATIENT_ARRAY=${ARG_ARRAY[@]:${BEGIN_INDEX}:${N_PATIENTS}} # ARG_ARRAY is indexed from 0!

BEGIN_INDEX=${BEGIN_INDEX}+${N_PATIENTS} # ARG_ARRAY is indexed from 0!
PATIENT_AGE_ARRAY=${ARG_ARRAY[@]:${BEGIN_INDEX}:${N_PATIENTS}} # ARG_ARRAY is indexed from 0!

BEGIN_INDEX=${BEGIN_INDEX}+${N_PATIENTS} # ARG_ARRAY is indexed from 0!
PATIENT_N_WGS_ARRAY=${ARG_ARRAY[@]:${BEGIN_INDEX}:${N_PATIENTS}} # ARG_ARRAY is indexed from 0!


JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

# RUNID='KX004'

TEMP_JOB_INDEX=$((${JOB_INDEX}-1))
TEMP_OBS_COUNT=$((${TEMP_JOB_INDEX}/${N_RUNS_PER_POST_OBS}))
TEMP_RESIM_RUN_INDEX=$((${TEMP_JOB_INDEX}%${N_RUNS_PER_POST_OBS}))

POST_OBS_INDEX=$((${TEMP_OBS_COUNT}+${BEGIN_POST_OBS}))
RESIM_RUN_INDEX=$((${TEMP_RESIM_RUN_INDEX}+1))
echo "POST_OBS_INDEX = ${POST_OBS_INDEX}"
echo "RESIM_RUN_INDEX = ${RESIM_RUN_INDEX}"

R_PARAMS="${SOURCE_PATH} ${OUTPUT_PATH} ${POST_SAMPLE_FILE} ${CONDITIONAL_RESIM_DIR_NAME} ${N_RESIMS_PER_RUN} ${POST_OBS_INDEX} ${RESIM_RUN_INDEX} ${N_PATIENTS} ${PATIENT_ARRAY[@]} ${PATIENT_AGE_ARRAY[@]} ${PATIENT_N_WGS_ARRAY[@]}"
echo "R_PARAMS = ${R_PARAMS}"


R_script_name="resim_conditionals_drivers.R"
Rout_filename="resim_conditionals"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}_POST_OBS_${POST_OBS_INDEX}_RUN_INDEX_${RESIM_RUN_INDEX}.Rout"
