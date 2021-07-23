#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
SIM_PATH=${3}
POOLED_SIM_PATH=${4}
POOLED_SIM_FILE_NAME=${5}
N_SIMS_PER_RUN=${6}
N_RUNS=${7}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

R_PARAMS="${SOURCE_PATH} ${SIM_PATH} ${POOLED_SIM_PATH} ${POOLED_SIM_FILE_NAME} ${N_SIMS_PER_RUN} ${N_RUNS}"

R_script_name="collect_prior_sims_drivers.R"
Rout_filename="collect_prior_sims_drivers"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
# export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}.Rout"

