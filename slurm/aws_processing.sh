#!/bin/bash -l

###
# SLURM OPTIONS
###

# Name of job
#SBATCH --job-name=ps-oct_slices

# Partition
#SBATCH --partition=agsmall

# Timing
#SBATCH --time=15:30:00

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

# Input options
# Options must be provided in this order from the command line
# e.g. sbatch example_array.sh NAME_OF_RCLONE Moe 3dtile 30

# Name of your rclone configuration for the bucket
RCLONE_NAME=$1
# Name of the monkey
SUBJECT_NAME=$2
# Can only be 'enface' or '3dtile'
# Whether you want to generate enface or 3d tiles
ENFACE_VS_3DTILE=$3
# The tile number you want to generate IFF generating 3D tiles 
3D_TILE_NUM=$4

CSV_FILE=${SUBJECT_NAME}.csv

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
#rclone mount "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Raw/${SUBJECT_NAME}/PS-OCT/${DIR_DATE}/" $MOUNT_PATH &
rclone mount "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/PS-OCT/${SUBJECT_NAME}/Raw/${DIR_DATE}/" $MOUNT_PATH &
sleep 5 # Takes rclone a second to actually mount

# Write out wrapper functions for a given slice
python3 pre-analysis_script.py --csvfile ${CSV_FILE} --slicenum ${SLURM_ARRAY_TASK_ID} --enface_vs_3dtile ${ENFACE_VS_3DTILE} 

# Launch the matlab code per slice
export MATLABPATH=/tmp/
module load matlab/R2019a
matlab -nodisplay -nodesktop -nosplash -r "run('/tmp/slice_${SLURM_ARRAY_TASK_ID}_wrapper.m'); exit;"

# 4) Write it back to the S3 bucket following bucket structure
# Bucket structure is different than how the data is saved to scratch.
# Do not want Orientation dir, or CDP, or A1A2 dirs.
SAVE_PATH=/scratch.local/PSOCT
module load rclone
rclone sync $SAVE_PATH/Stitched/AbsoOri/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/Enface/Orientation/"
rclone sync $SAVE_PATH/Stitched/Cross/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/Enface/Cross/"
rclone sync $SAVE_PATH/Stitched/Reflectivity/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/Enface/Reflectivity/"
rclone sync $SAVE_PATH/Stitched/Retardance/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/Enface/Retardance/"

rclone sync $SAVE_PATH/jpegs/Cross/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/jpegs/Cross/"
rclone sync $SAVE_PATH/jpegs/Reflectivity/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/jpegs/Reflectivity/"
rclone sync $SAVE_PATH/jpegs/Retardance/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/jpegs/Retardance/"

rclone sync $SAVE_PATH/TComp/Cross/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/Cross/"
rclone sync $SAVE_PATH/TComp/Reflectivity/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/Reflectivity/"
rclone sync $SAVE_PATH/TComp/Retardance/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/Retardance/"
rclone sync $SAVE_PATH/TComp/AbsoOri/ "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/Orientation/"

kill %1
fusermount3 -u /tmp/cmc-s3-bucket
