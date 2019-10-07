#!/bin/bash
#
#SBATCH -p wmglab
#SBATCH --job-name=ctmsTC
#SBATCH -c 36
#SBATCH --output=job_out/ctmsTC.out
#SBATCH --error=job_out/ctmsTC.err
#SBATCH --mem=30G
#SBATCH --mail-type=END,FAIL
uname -n 1>&1 # print node name
# Model parameters
nrn_model_ver="maxH"
tms_mode=4 # cTMS1 waveforms
durs=(30 60 120)
cTMS_durs="[${durs[*]}]"
lsn=1 # layer set number
nrn_pop_name="nrn_pop1-nrn_pop6_all"
Efield_solution="M1_PA_Magstim70mm"
data_files=()
for i in "${durs[@]}"; do
    data_files+=("'tms_${nrn_model_ver}_w${tms_mode}_${i}us_ls_${lsn}_E_${Efield_solution}_P_${nrn_pop_name}'")
done
data_files_str="{${data_files[*]}}"
data_output_file="tms_${nrn_model_ver}_w${tms_mode}_${durs[0]}-${durs[-1]}us_ls_${lsn}_E_${Efield_solution}_P_${nrn_pop_name}_timeconstants" 
# Go to MATLAB directory and launch function
main_fold="/dscrhome/asa45/TMSsim_Aberra2019"
cd ${main_fold}/mat

echo "Launching matlab"
mat_args="${cTMS_durs},${data_files_str},'${data_output_file}'"
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "addPaths; ctmsTimeConstants(${mat_args});"

echo "Finished"