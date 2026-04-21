#!/bin/bash -l

###
# SLURM OPTIONS
###

# Name of job
#SBATCH --job-name=ps-oct_slices

# Partition
#SBATCH --partition=agsmall

# Timing
#SBATCH --time=00:30:00

# Mem per node request
# In testing, used max of 40G
#SBATCH --mem=40G

# Request a specific number of cores
# per slice aka task
#SBATCH --cpus-per-task=10

# %A is job number and %a is array index
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

# Must set mail-type to ARRAY_TASKS to get notified per array job and not entire set
## SBATCH --mail-type=ARRAY_TASKS
## SBATCH --mail-type=END,FAIL
## SBATCH --mail-user=huxfo013@umn.edu

# Set to the slice numbers you want to analyze
# Can give as 1-3 for range e.g. 1,2,3
# OR give as 1,5,7 e.g. for particular slices
#SBATCH --array=1

# Scratch Space request
# Tune for slice range listed above
# The max space required=500GB (max size for 100 tiles)
#SBATCH --tmp=500G


###
# Processing
###

# Fetch relevant code from github
git clone https://github.com/rhuxfo/midb_cmc_ps-oct.git /tmp/midb_cmc_ps-oct
cp /tmp/midb_cmc_ps-oct/main_codes/* /tmp/
cp /tmp/midb_cmc_ps-oct/slurm/Test.py /tmp/

# Actually copy data to local scratch
module load conda
source activate EnvZarr_py312

# Launch the matlab code per slice
export MATLABPATH=/tmp/
module load matlab
matlab -nodisplay -nodesktop -nosplash -r "run(pyrunfile("Test.py")); exit;"
