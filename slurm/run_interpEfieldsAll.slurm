#!/bin/bash
#
#SBATCH -p wmglab
#SBATCH -o job_out/interpE.out
#SBATCH -e job_out/interpE.err
#SBATCH --job-name=interpE
#SBATCH -c 36
#SBATCH --mail-type=END,FAIL

# run interpEfieldsAll()
# calls interpEfield for all 6 neuron populations (nrn_pop1-6) and interpolates E-field for both
# E-field solutions in used in Aberra et al. 2019

cd ../mat # go to matlab code directory
/opt/apps/matlabR2016a/bin/matlab -nodisplay -nodesktop -nosplash -r "interpEfieldsAll()"
