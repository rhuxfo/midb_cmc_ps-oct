#!/bin/bash -l
###
# SLURM OPTIONS
###
# Name of job
#SBATCH --job-name=lincconvert
# Partition
#SBATCH --partition=msismall
# Timing
#SBATCH --time=20:30:00

# Mem per node request
# FIXME Need to figure out how much it needs
#SBATCH --mem=100G

# Request a specific number of cores
# per slice aka task
#SBATCH --cpus-per-task=32

# %A is job number and %a is array index
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

# Set to the slice numbers you want to analyze

# Can give as 1-3 for range e.g. 1,2,3
# OR give as 1,5,7 e.g. for particular slices

#SBATCH --array=102-109

# Scratch Space request
# Tune for slice range listed above
# The max space required=500GB (max size for 100 tiles)
# FIXME How big is one zarr file?

#SBATCH --tmp=250G

###

# Processing
###
# Input options
# Options must be provided in this order from the command line
# e.g. sbatch LC.sh NAME_OF_RCLONE Moe EnvZarr_py312 Cross

###
# INPUTS
###
# Name of your rclone configuration for the bucket
RCLONE_NAME=$1

# Name of the monkey
SUBJECT_NAME=$2

# Conda environment
ENV=$3

# Contrast to be stitched
CONTRAST=$4

SUB_CON=$5

# CSV filename
CSV_FILE=${SUBJECT_NAME}.csv

# Out directory
PROJECT_SPACE=$MSIPROJECT/shared/lincZarr

###
# SETUP
###
# Fetch the write date from the csv sheet

git clone https://github.com/rhuxfo/midb_cmc_ps-oct.git /tmp/midb_cmc_ps-oct
cp /tmp/midb_cmc_ps-oct/slurm/${CSV_FILE} ./

tile_x=$(awk -F, -e '$2=='${SLURM_ARRAY_TASK_ID}' { print $3 }' <${CSV_FILE} )
tile_y=$(awk -F, -e '$2=='${SLURM_ARRAY_TASK_ID}' { print $4 }' <${CSV_FILE} )

# Load linc-convert environment

module load conda
source activate ${ENV}
echo ${ENV}
# Actually copy data to local scratch

module load rclone
MOUNT_PATH=/tmp/cmc-s3-bucket

mkdir $MOUNT_PATH

# FIXME not sure what path this should be

rclone mount "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/${CONTRAST}/" $MOUNT_PATH --read-only &
sleep 5 # Takes rclone a second to actually mount

echo "Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/${CONTRAST}/"
###
# LINC-CONVERT
###
# Name the yaml file
YAML_FILE=${PROJECT_SPACE}/${SUBJECT_NAME}_stitch.yaml

# Generate the tile config file
linc-convert psoct generate_tile_config \
    --columns ${tile_y} \
    --rows ${tile_x} \
    --tile-size 1000 \
    --base-dir $MOUNT_PATH \
    --naming-format slice_${SLURM_ARRAY_TASK_ID}_tile_{tile_number:03d}_${SUB_CON}.mat \
    --output ${YAML_FILE} \
    --grid-type row-by-row\
    --order right-down \
    --overlap-percentage 0.1

# Actually do the stitching

linc-convert psoct mosaic ${YAML_FILE} \
    --out ${PROJECT_SPACE}/${SUBJECT_NAME}_slice_${SLURM_ARRAY_TASK_ID}_${CONTRAST}.zarr \
    --zarr_version 3 \
    --shard 1024 \
    --chunk 256 \
    --driver tensorstore \
    --tile-overlap 0.1 \
    --overwrite \
    --voxel-size 5 --voxel-size 5 --voxel-size 5 --verbose

# Copy to bucket

rclone copy ${PROJECT_SPACE}/${SUBJECT_NAME}_slice_${SLURM_ARRAY_TASK_ID}_${CONTRAST}.zarr "${RCLONE_NAME}:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/Zarr/${CONTRAST}/${SUBJECT_NAME}_slice_${SLURM_ARRAY_TASK_ID}_${CONTRAST}/" --s3-no-check-bucket --transfers 32 --checkers 32 --progress
###
# CLEANUP
###
kill %1
fusermount3 -u /tmp/cmc-s3-bucket#!/bin/bash -l

