#!/bin/bash -l        

###
# SLURM OPTIONS
###

# Name of job
#SBATCH --job-name=ps-oct_slices

# Partition
# SBATCH --partition=agsmall

# Mem per node request
# In testing, used max of 40G
#SBATCH --mem=40G

# Request a specific number of cores
# per slice aka task
#SBATCH --cpus-per-task=6

# Scratch Space request=500GB (max size for 100 tiles)
#SBATCH --tmp=500000

# Must set mail-type to ARRAY_TASKS to get notified per array job and not entire set
#SBATCH --mail-type=ARRAY_TASKS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=huxfo013@umn.edu

# Set to the slice numbers you want to analyze
#SBATCH --array=156-157

# %A is job number and %a is array index
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

###
# Processing
###

# Fetch the write date from the csv sheet
mydate=$(awk -F, -e '$2=='${SLURM_ARRAY_TASK_ID}' { print $1 }' <Kumquat_Right_Hemisphere_Slices_Sheet2.csv )
dir_date=$(date -d "$mydate" +%m%d%Y)

# Actually copy data to local scratch
module load s5cmd
s5cmd sync 's3://midb-cmc-nonhuman/PS-OCT/KQRH/Raw/'${dir_date}'/Slice_'${SLURM_ARRAY_TASK_ID}'_Tile_*_840_*.dat' /scratch.local/ 

# Write out wrapper functions for a given slice
python3 pre-analysis_script.py ${SLURM_ARRAY_TASK_ID} 
s3cmd get --skip-existing s3://midb-cmc-nonhuman/PS-OCT/Code/Current\ Processing\ code/ComTom_W_Ch1_shifted.dat /scratch.local/
s3cmd get --skip-existing s3://midb-cmc-nonhuman/PS-OCT/Code/Current\ Processing\ code/ComTom_W_Ch2_shifted.dat /scratch.local/

# Launch the matlab code per slice
module load matlab/R2019a
matlab -nodisplay -nodesktop -nosplash -r "run('/scratch.local/slice_${SLURM_ARRAY_TASK_ID}_wrapper.m'); exit;"

# 4) Write it back to the S3 bucket following bucket structure
# Bucket structure is different than how the data is saved to scratch.
# Do not want Orientation dir, or CDP, or A1A2 dirs.
s5cmd sync /scratch.local/Stitched/AbsoOri/ 's3://midb-cmc-nonhuman/PS-OCT/KQRH/Enface/Orientation/'
s5cmd sync /scratch.local/Stitched/Cross/ 's3://midb-cmc-nonhuman/PS-OCT/KQRH/Enface/Cross/'
s5cmd sync /scratch.local/Stitched/Reflectivity/ 's3://midb-cmc-nonhuman/PS-OCT/KQRH/Enface/Reflectivity/'
s5cmd sync /scratch.local/Stitched/Retardance/ 's3://midb-cmc-nonhuman/PS-OCT/KQRH/Enface/Retardance/'

