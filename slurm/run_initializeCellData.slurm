#!/bin/bash
#
#SBATCH -p wmglab
#SBATCH -o job_out/init_data.out
#SBATCH -e job_out/init_data.err
#SBATCH --job-name=init_tms
#SBATCH -c 36
#SBATCH --mail-type=END,FAIL

# run initialize_cell_data()
# generates cell model data and 6 neuron populations

cd ../mat # go to matlab code directory
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "initializeCellData()"
