#!/bin/bash
#
#SBATCH -p wmglab, common
#SBATCH --mail-type=FAIL
uname -n 1>&1 # print node name
# Model parameters
cell_id=$SLURM_ARRAY_TASK_ID
# passed in from submit_simTMSlayerThresh.sh:
# nrn_model_ver,tms_mode,layer_set_num,nrn_pop_name,
# Efield_solution, model_prefix, and over_write 
if [ -z "$replace_axon" ]; then # if not exported from 'submit' script 
    if [ $cell_id -gt 15 ]; then # cell_ids 16-25 correspond to L5 and L6 PCs
        replace_axon=2 # disable axon terminal    
    else
        replace_axon=0 # leave axon intact
    fi
fi
echo "replace_axon = ${replace_axon}" 
if [ $SLURM_ARRAY_TASK_ID -eq $cell_id_start ]
then
    echo "Saving SLURM script"
    # Copy SLURM script to mat/nrn_sim_data/<model_prefix> folder
    mkdir -p ${main_fold}/mat/nrn_sim_data/${model_prefix}/ && cp -R -u -p ./run_simTMSlayerThresh.slurm ${main_fold}/mat/nrn_sim_data/${model_prefix}
    # Copy SLURM script to $job_fold folder
    cp -R -u -p ./run_simTMSlayerThresh.slurm ${job_fold}    
fi     

# make special executable
chmod +x ${main_fold}/nrn/special

# Go to MATLAB directory and launch threshold loop
cd ${main_fold}/mat

echo "Launching matlab"
mat_args="${cell_id},'${model_prefix}','nrn_model_ver','${nrn_model_ver}','layer_set_num',${layer_set_num},\
         'Efield_solution','${Efield_solution}','nrn_pop_name','${nrn_pop_name}',\
         'mode',${tms_mode},'replace_axon',${replace_axon},'over_write',${over_write}"
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "addPaths;; simTMSlayerThresh($mat_args);"

echo "Finished"