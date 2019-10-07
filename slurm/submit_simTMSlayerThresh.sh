#!/bin/bash
# Submit simTMSlayerThresh jobs for all cells
export main_fold=".." # or replace with absolute path to TMSsim_Aberra2019/
numCPUs=9
total_mem=30G
job_name="tms"
slurm_file="run_simTMSlayerThresh.slurm" # tms_thresh.slurm or tms_thresh_ids.slurm
export nrn_model_ver="maxH"
export tms_mode=1 # MagProx100 monophasic
export layer_set_num=1 # layer set number
export nrn_pop_name="nrn_pop1"
export Efield_solution="M1_PA_MCB70"
export model_prefix="tms_${nrn_model_ver}_w${tms_mode}_ls_${layer_set_num}_E_${Efield_solution}_P_${nrn_pop_name}" 
export cell_id_start=1
export cell_id_end=25
export over_write=0
# make job output folder  
export job_fold="job_out/${model_prefix}"
mkdir -p $job_fold
echo sbatch --array=${cell_id_start}-${cell_id_end} \
    --output=${job_fold}/cell_%a.out \
    --error=${job_fold}/cell_%a.err \
    --job-name=${job_name} --cpus-per-task=$numCPUs --mem=${total_mem} ${slurm_file}
# submit job array
sbatch --array=${cell_id_start}-${cell_id_end} \
    --output=${job_fold}/cell_%a.out \
    --error=${job_fold}/cell_%a.err \
    --job-name=${job_name} --cpus-per-task=${numCPUs} --mem=${total_mem} ${slurm_file}

# run run_combineLayerThreshData.slurm to compile data into one .mat file