#!/bin/bash -l

###
# SLURM OPTIONS
###

# Name of job
#SBATCH --job-name=ps-oct_check

# Partition
#SBATCH --partition=agsmall

# Timing
#SBATCH --time=1:30:00

# Mem per node request
# In testing, used max of 40G
#SBATCH --mem=10G

# Request a specific number of cores
# per slice aka task
#SBATCH --cpus-per-task=1

# %A is job number and %a is array index
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

# Must set mail-type to ARRAY_TASKS to get notified per array job and not entire set
#SBATCH --mail-type=ARRAY_TASKS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=huxfo013@umn.edu

# Set to the slice numbers you want to analyze
# Can give as 1-3 for range e.g. 1,2,3
# OR give as 1,5,7 e.g. for particular slices
#SBATCH --array=1

# Scratch Space request
# Tune for slice range listed above
# The max space required=500GB (max size for 100 tiles)
#SBATCH --tmp=10G


###
# Processing
###

SLST=$1
SLEND=$2
RCLONE_NAME=cmcs3
SUBJECT_NAME=Zebel

# Fetch relevant code from github
git clone https://github.com/rhuxfo/midb_cmc_ps-oct.git /tmp/midb_cmc_ps-oct
cp /tmp/midb_cmc_ps-oct/main_codes/* /tmp/

# Actually copy data to local scratch
module load rclone
MOUNT_PATH=/tmp/cmc-s3-bucket
mkdir $MOUNT_PATH
rclone mount "${RCLONE_NAME}:midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/Enface/Reflectivity/" $MOUNT_PATH &
sleep 5 # Takes rclone a second to actually mount

mkdir /tmp/Reflectivity/

# Launch the matlab code per slice
export MATLABPATH=/tmp/
module load matlab/R2019a
matlab -nodisplay -nodesktop -nosplash -r "Dir= '${MOUNT_PATH}'; Sdir = '/tmp/Reflectivity'; SL =${SLST}:${SLEND}; gifStack(Dir,Sdir,SL);  exit;"

# 4) Write it back to the S3 bucket following bucket structure
# Bucket structure is different than how the data is saved to scratch.
# Do not want Orientation dir, or CDP, or A1A2 dirs.
module load s5cmd
#s5cmd sync /tmp/Stitched/AbsoOri/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Orientation/"
#s5cmd sync /tmp/Stitched/Cross/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Cross/"
s5cmd sync /tmp/Reflectivity/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Reflectivity/"
#s5cmd sync /tmp/Stitched/Retardance/ "s3://midb-cmc-nonhuman/PS-OCT/${SUBJECT_NAME}/jpegs/Retardance/"

kill %1
fusermount3 -u /tmp/cmc-s3-bucket
