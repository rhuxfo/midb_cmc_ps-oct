#!/bin/bash -l        

# Name of job
#SBATCH --job-name=cmc_ps-oct_slice_analysis

# Must set mail-type to ARRAY_TASKS to get notified per array job and not entire set
#SBATCH --mail-type=ARRAY_TASKS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=huxfo013@umn.edu

# Set to the slice numbers you want to analyze
#SBATCH --array=156-157

# %A is job number and %a is array index
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

# Actual call for each array job
# 1) use the slice number to copy appropriate raw data
# 2) Given a txt file and a slice number, write out the wrapper function.
mydate=$(awk -F, -e '$2=='${SLURM_ARRAY_TASK_ID}' { print $1 }' <Kumquat_Right_Hemisphere_Slices_Sheet2.csv )
dir_date=$(date -d "$mydate" +%m%d%Y)
module load s5cmd
s5cmd sync 's3://midb-cmc-nonhuman/PS-OCT/KQRH/Raw/'${dir_date}'/Slice_'${SLURM_ARRAY_TASK_ID}'_Tile_*_840_*.dat\' /scratch.local/ 

python3 pre-analysis_script.py ${SLURM_ARRAY_TASK_ID} 
s3cmd get --skip-existing s3://midb-cmc-nonhuman/PS-OCT/Code/Current\ Processing\ code/ComTom_W_Ch1_shifted.dat /scratch.local/
s3cmd get --skip-existing s3://midb-cmc-nonhuman/PS-OCT/Code/Current\ Processing\ code/ComTom_W_Ch2_shifted.dat /scratch.local/

#3) Actually launch the matlab code with the parfor loop
module load matlab/R2019a
matlab -nodisplay -nodesktop -nosplash -r "run('/scratch.local/slice_${SLURM_ARRAY_TASK_ID}_wrapper.m'); exit;"

# 4) Write it back to the S3 bucket
s3cmd sync /scratch.local/Enface/ s3://midb-cmc-nonhuman/PS-OCT/KQRH/Enface/
