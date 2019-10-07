#!/bin/bash
#
#SBATCH -p wmglab
#SBATCH --job-name=combD
#SBATCH -c 36
#SBATCH --output=job_out/comb_data.out
#SBATCH --error=job_out/comb_data.err
#SBATCH --mem=5G
#SBATCH --mail-type=END,FAIL
uname -n 1>&1 # print node name
# Model parameters
nrn_model_ver="maxHlin"
tms_mode=1 # MagProx100 monophasic
lsn=1 # layer set number
nrn_pop_name="nrn_pop1"
Efield_solution="M1_PA_MCB70"
model_prefix="tms_${nrn_model_ver}_w${tms_mode}_ls_${lsn}_E_${Efield_solution}_P_${nrn_pop_name}" 
cell_ids='[6:10]'
# Go to MATLAB directory and launch function
main_fold=".." # or replace with absolute path to TMSsim_Aberra2019/
cd ${main_fold}/mat

echo "Launching matlab"
mat_args="'${model_prefix}',${cell_ids}"
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "addPaths; combineLayerThreshData($mat_args);"

echo "Finished"