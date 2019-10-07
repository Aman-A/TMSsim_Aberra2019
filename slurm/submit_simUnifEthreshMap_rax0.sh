#!/bin/bash
# Submit simUnifEthreshMap_rax0 jobs for L5/6 PCs
export main_fold=".." # or replace with absolute path to TMSsim_Aberra2019/
numCPUs=36
total_mem=30G
job_name="utms"
slurm_file="run_simUnifEthreshMap.slurm" # tms_thresh.slurm or tms_thresh_ids.slurm
export nrn_model_ver="maxH"
export cell_id_list=16-25
export replace_axon=0 # leave axon intact, instead of disabling terminal
export cell_id_start=16
export over_write=0
export model_prefix="utms_${nrn_model_ver}_w1_rax${replace_axon}"
# make job output folder  
export job_fold="job_out/${model_prefix}"
mkdir -p $job_fold
echo sbatch --array=${cell_id_list} \
    --output=${job_fold}/cell_%a.out \
    --error=${job_fold}/cell_%a.err \
    --job-name=${job_name} --cpus-per-task=$numCPUs --mem=${total_mem} ${slurm_file}
# submit job array
sbatch --array=${cell_id_list} \
    --output=${job_fold}/cell_%a.out \
    --error=${job_fold}/cell_%a.err \
    --job-name=${job_name} --cpus-per-task=${numCPUs} --mem=${total_mem} ${slurm_file}

# run run_combineUnifThreshData.slurm to compile data into one .mat file