#!/bin/bash -l

###
# SLURM OPTIONS
###

# Name of job
#SBATCH --job-name=ps-oct_slices

# Partition
#SBATCH --partition=agsmall

# Timing
#SBATCH --time=7:30:00

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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yeatt002@umn.edu

# Set to the slice numbers you want to analyze
# Can give as 1-3 for range e.g. 1,2,3
# OR give as 1,5,7 e.g. for particular slices
#SBATCH --array=1

# Scratch Space request
#SBATCH --tmp=10G


###
# Processing
###

# Fetch the write date from the csv sheet
mydate=$(awk -F, -e '$2=='${SLURM_ARRAY_TASK_ID}' { print $1 }' <Kumquat_Right_Hemisphere_Slices_Sheet2.csv )
DIR_DATE=$(date -d "$mydate" +%m%d%Y)
echo $DIR_DATE

# Actually copy data to local scratch
module load rclone
MOUNT_PATH=/tmp/cmc-s3-bucket
mkdir $MOUNT_PATH
rclone mount "ceph:midb-cmc-nonhuman/PS-OCT/KQRH/Raw/${DIR_DATE}/" $MOUNT_PATH &
sleep 5 # Takes rclone a second to actually mount

# Write out wrapper functions for a given slice
python3 pre-analysis_script.py ${SLURM_ARRAY_TASK_ID} 
git clone https://github.com/rhuxfo/midb_cmc_ps-oct.git /tmp/midb_cmc_ps-oct
cp /tmp/midb_cmc_ps-oct/main_codes/* /tmp/

# Launch the matlab code per slice
export MATLABPATH=/tmp/
module load matlab/R2019a
matlab -nodisplay -nodesktop -nosplash -r "run('/tmp/slice_${SLURM_ARRAY_TASK_ID}_wrapper.m'); exit;"

# 4) Write it back to the S3 bucket following bucket structure
# Bucket structure is different than how the data is saved to scratch.
# Do not want Orientation dir, or CDP, or A1A2 dirs.
module load s5cmd
s5cmd sync /tmp/Stitched/AbsoOri/ 's3://midb-cmc-nonhuman/PS-OCT/Moe/Enface/Orientation/'
s5cmd sync /tmp/Stitched/Cross/ 's3://midb-cmc-nonhuman/PS-OCT/Moe/Enface/Cross/'
s5cmd sync /tmp/Stitched/Reflectivity/ 's3://midb-cmc-nonhuman/PS-OCT/Moe/Enface/Reflectivity/'
s5cmd sync /tmp/Stitched/Retardance/ 's3://midb-cmc-nonhuman/PS-OCT/Moe/Enface/Retardance/'

kill %1
fusermount3 -u /tmp/cmc-s3-bucket
