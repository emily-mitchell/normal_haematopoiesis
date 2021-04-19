#!/bin/bash
SOURCE_PATH=${1}
LOGS_PATH=${2}
PP_SIMS_FILE=${3}
RESIM_PATH=${4}
N_RUNS_PER_PP_SIM=${5}
N_RESIMS_PER_RUN=${6}
JOB_INDEX=${7}

# JOB_INDEX=${LSB_JOBINDEX}
JOB_ID=${LSB_JOBID}

RUNID='KX003'

TEMP_JOB_INDEX=$((${JOB_INDEX}-1))
TEMP_PP_SIM_INDEX=$((${TEMP_JOB_INDEX}/${N_RUNS_PER_PP_SIM}))
TEMP_RESIM_RUN_INDEX=$((${TEMP_JOB_INDEX}%${N_RUNS_PER_PP_SIM}))

PP_SIM_INDEX=$((${TEMP_PP_SIM_INDEX}+1))
RESIM_RUN_INDEX=$((${TEMP_RESIM_RUN_INDEX}+1))
echo "PP_SIM_INDEX = ${PP_SIM_INDEX}"
echo "RESIM_RUN_INDEX = ${RESIM_RUN_INDEX}"

R_PARAMS="${PP_SIMS_FILE} ${RESIM_PATH} ${PP_SIM_INDEX} ${RESIM_RUN_INDEX} ${N_RESIMS_PER_RUN} ${SOURCE_PATH}"
echo "R_PARAMS = ${R_PARAMS}"


R_script_name="KX003_resim_conditionals.R"
Rout_filename="KX003_resim_conditionals"

export PATH=/software/R-3.6.1/bin:${PATH}
export R_HOME=$(R RHOME)
export R_LIBS_USER="/nfs/users/nfs_e/em16/R/x86_64-pc-linux-gnu-library/3.6"
#export R_LIBS_USER="/nfs/users/nfs_k/kd7/R/x86_64-pc-linux-gnu-library/3.6"

# CMD="/software/bin/R-dev CMD BATCH"
CMD="R CMD BATCH"

cd ${SOURCE_PATH}

# Rscript "--no-save --no-restore --verbose ./KX003_resim_conditionals.R ${R_PARAMS} > KX003_resim_conditionals.Rout"
# Rscript "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}_PP_SIM_${PP_SIM_INDEX}_RUN_INDEX_${RESIM_RUN_INDEX}.Rout"
${CMD} "--no-save --no-restore --args ${R_PARAMS}" ${SOURCE_PATH}${R_script_name} "${LOGS_PATH}${Rout_filename}_ID_${JOB_ID}_PP_SIM_${PP_SIM_INDEX}_RUN_INDEX_${RESIM_RUN_INDEX}.Rout"


