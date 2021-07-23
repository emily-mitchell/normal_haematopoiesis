#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
PATIENT_ID=${3}
PATIENT_AGE=${4}
PATIENT_N_WGS=${5}
SIM_PATH=${6}
PRIOR_SAMPLE_FILE=${7}
USE_PRIOR_SAMPLE_BOOL=${8}
N_SIMS_PER_RUN=${9}

JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

# RUNID='KX004'

SIM_RUN_INDEX=${JOB_INDEX}
echo "SIM_RUN_INDEX = ${SIM_RUN_INDEX}"

RAND_SEED=${JOB_INDEX} # Not used!

R_PARAMS="${SOURCE_PATH} ${SIM_PATH} ${PRIOR_SAMPLE_FILE} ${USE_PRIOR_SAMPLE_BOOL} ${PATIENT_AGE} ${PATIENT_N_WGS} ${N_SIMS_PER_RUN} ${SIM_RUN_INDEX} ${RAND_SEED}"
echo "R_PARAMS = ${R_PARAMS}"


R_script_name="sim_prior_drivers.R"
Rout_filename="sim_prior_drivers"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_PATIENT_${PATIENT_ID}_JOB_${JOB_ID}_RUN_INDEX_${SIM_RUN_INDEX}.Rout"

