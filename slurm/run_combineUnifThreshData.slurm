#!/bin/bash
#
#SBATCH -p wmglab
#SBATCH --job-name=combD
#SBATCH -c 1
#SBATCH --output=job_out/comb_data.out
#SBATCH --error=job_out/comb_data.err
#SBATCH --mem=5G
#SBATCH --mail-type=END,FAIL
uname -n 1>&1 # print node name
# Model parameters
model_prefix="utms_maxHlin_w1_rax5"
cell_ids='[6:10]'
# Go to MATLAB directory and launch function
main_fold=".." # or replace with absolute path to TMSsim_Aberra2019/
cd ${main_fold}/mat

echo "Launching matlab"
mat_args="'${model_prefix}',${cell_ids}"
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "addPaths; combineUnifThreshData($mat_args);"

echo "Finished"