#!/bin/bash
#
#SBATCH -p wmglab, common
#SBATCH --mail-type=FAIL
uname -n 1>&1 # print node name
# Model parameters
cell_id=$SLURM_ARRAY_TASK_ID
# model_prefix and over_write passed in from submit_simUnifEthreshMap.sh
if [ -z "$nrn_model_ver" ]; then 
    nrn_model_ver="maxH" # default neuron model version
fi
tms_mode=1 # MagProX100 - monophasic
dtheta=15 # (deg) polar angle step size
dphi=10 # (deg) azimuthal angle step size
if [ -z "$replace_axon" ]; then # if not exported from 'submit' script 
    if [ $cell_id -gt 15 ]; then # cell_ids 16-25 correspond to L5 and L6 PCs
        replace_axon=2 # disable axon terminal    
    else
        replace_axon=0 # leave axon intact
    fi
fi
echo "replace_axon = ${replace_axon}" 
if [ ${cell_id} -eq ${cell_id_start} ]; then
    echo "Saving SLURM script"
    # Copy SLURM script to mat/nrn_sim_data/<model_prefix> folder
    mkdir -p ${main_fold}/mat/nrn_sim_data/${model_prefix}/${nrn_model_ver} && cp -R -u -p ./run_simUnifEthreshMap.slurm ${main_fold}/mat/nrn_sim_data/${model_prefix}
    # Copy SLURM script to $job_fold folder
    cp -R -u -p ./run_simUnifEthreshMap.slurm ${job_fold}    
fi     

# make special executable
chmod +x ${main_fold}/nrn/special

# Go to MATLAB directory and launch threshold loop
cd ${main_fold}/mat

echo "Launching matlab"
mat_args="${cell_id},'${model_prefix}','nrn_model_ver','${nrn_model_ver}',\
         'mode',${tms_mode},'dtheta',${dtheta},'dphi',${dphi},'over_write',\
         ${over_write},'replace_axon',${replace_axon}"
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "addPaths; simUnifEthreshMap($mat_args);"

echo "Finished"