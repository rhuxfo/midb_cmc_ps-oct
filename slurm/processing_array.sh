#!/bin/bash -l

###
# SLURM OPTIONS
###

# Name of job
#SBATCH --job-name=ps-oct_slices

# Partition
#SBATCH --partition=agsmall

# Timing
#SBATCH --time=5:30:00

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
#SBATCH --tmp=300G


###
# Processing
###

# Input options
# Options must be provided in this order from the command line
# e.g. sbatch example_array.sh Moe.csv NAME_OF_RCLONE Moe
CSV_FILE=$1
RCLONE_NAME=$2
SUBJECT_NAME=$3

# Fetch relevant code from github
git clone https://github.com/rhuxfo/midb_cmc_ps-oct.git /tmp/midb_cmc_ps-oct
cp /tmp/midb_cmc_ps-oct/main_codes/* /tmp/
cp /tmp/midb_cmc_ps-oct/slurm/${CSV_FILE} ./
cp /tmp/midb_cmc_ps-oct/slurm/pre-analysis_script.py ./

# Fetch the write date from the csv sheet
mydate=$(awk -F, -e '$2=='${SLURM_ARRAY_TASK_ID}' { print $1 }' <${CSV_FILE} )
DIR_DATE=$(date -d "$mydate" +%m%d%Y)
echo $DIR_DATE

# Actually copy data to local scratch
module load rclone
MOUNT_PATH=/tmp/cmc-s3-bucket
mkdir $MOUNT_PATH
rclone mount "${RCLONE_NAME}:midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/Raw/${DIR_DATE}/" $MOUNT_PATH &
sleep 5 # Takes rclone a second to actually mount

# Write out wrapper functions for a given slice
python3 pre-analysis_script.py --csvfile ${CSV_FILE} --slicenum ${SLURM_ARRAY_TASK_ID} 

# Launch the matlab code per slice
export MATLABPATH=/tmp/
module load matlab/R2019a
matlab -nodisplay -nodesktop -nosplash -r "run('/tmp/slice_${SLURM_ARRAY_TASK_ID}_wrapper.m'); exit;"

# 4) Write it back to the S3 bucket following bucket structure
# Bucket structure is different than how the data is saved to scratch.
# Do not want Orientation dir, or CDP, or A1A2 dirs.
SAVE_PATH=/scratch.global/PSOCT
module load s5cmd
s5cmd sync $SAVE_PATH/Stitched/AbsoOri/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/Enface/Orientation/"
s5cmd sync $SAVE_PATH/Stitched/Cross/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/Enface/Cross/"
s5cmd sync $SAVE_PATH/Stitched/Reflectivity/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/Enface/Reflectivity/"
s5cmd sync $SAVE_PATH/Stitched/Retardance/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/Enface/Retardance/"

s5cmd sync $SAVE_PATH/jpegs/Cross/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Cross/"
s5cmd sync $SAVE_PATH/jpegs/Reflectivity/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Reflectivity/"
s5cmd sync $SAVE_PATH/jpegs/Retardance/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Retardance/"

s5cmd sync $SAVE_PATH/TComp/Cross/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/3DTiles/Cross/"
s5cmd sync $SAVE_PATH/TComp/Reflectivity/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/3DTiles/Reflectivity/"
s5cmd sync $SAVE_PATH/TComp/Retardance/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/3DTiles/Retardance/"
s5cmd sync $SAVE_PATH/TComp/AbsoOri/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/3DTiles/Orientation/"

kill %1
fusermount3 -u /tmp/cmc-s3-bucket
