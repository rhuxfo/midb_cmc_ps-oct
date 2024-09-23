#!/usr/bin/python3

#
# Read a given number of slices from a given monkey
# Store each slice as its own group
# Write out zarr file
#
# Note: No current PS-OCT code to stitch slices
#       but the "origin" is stored in csv files for each slice
#

import scipy.io as io
import numpy as np
import ome_zarr
from ome_zarr.writer import write_image
from ome_zarr.reader import Reader
from ome_zarr.io import parse_url
import ome_zarr.reader
import zarr

# Path to mounted s3 bucket
# must mount with rclone mount ceph:midb-cmc-nonhuman/PS-OCT/ /scratch.local/cmc-s3-bucket & before running this program.
monkey_name = 'KQRH'
path_to_bucket = f'/scratch.local/cmc-s3-bucket/{monkey_name}/Enface/'
contrasts = {'Cross': 'EnCr', 'Orientation': 'EnAO', 'Reflectivity': 'EnRef', 'Retardance': 'EnR'}
slice_range = range(1, 10)

# Initialize zarr file
zs = zarr.storage.FSStore(f'{monkey_name}_{min(slice_range)}_{max(slice_range)}.zarr')

for slice_num in slice_range:
    slice_name = 'Slice_' + str(slice_num)
    this_idx = len(contrasts)-1

    for contrast in contrasts:
        this_con_dir = ''.join([path_to_bucket, contrast])
        try:
            this_contrast = io.loadmat('/'.join([this_con_dir, slice_name+'_'+contrasts[contrast]+'.mat']))
        except FileNotFoundError:
            print(f'The contrast {contrast} does not exist for Slice {slice_num}')
            continue

        # Each Enface slice is a dict with the data entry being TEnCr, TEnAO, etc.
        this_contrast = this_contrast['T' + contrasts[contrast]]

        # Build the cxy array
        try:
            loaded_contrasts
        except NameError:
            loaded_contrasts = np.zeros((len(contrasts), this_contrast.shape[0], this_contrast.shape[1]))

        loaded_contrasts[this_idx, :, :] = this_contrast
        this_idx -= 1

    # Create Zarr group
    zg1 = zarr.group(store=zs, path=slice_name, overwrite=True)

    # Add Zarr group metadata
    axes_metadata = [{"name": "c", "type": "channel"}, {"name": "y", "type": "space", "unit": "pixel"}, {"name": "x", "type": "space", "unit": "pixel"}]

    # Write OME-Zarr image
    write_image(image=loaded_contrasts, group=zg1, axes=axes_metadata)

    # Slices are not the same size. Free loaded_contrasts to get right size for next slice.
    del loaded_contrasts
    print(f'{slice_name} stored as group')

print(f'{monkey_name} zarr file written. Reading back in ...')

###
# Read it back in
###

url = f'{monkey_name}_{min(slice_range)}_{max(slice_range)}.zarr/Slice_1'

# It is finding the file, otherwise it fails on this line.
reader = Reader(parse_url(url))

# nodes may include images, labels etc
nodes = list(reader())
print(ome_zarr.reader.Multiscales(nodes[0]))
print(nodes)
print(nodes[0])
# first node will be the image pixel data
image_node = nodes[0]
print(image_node.metadata)

dask_data = image_node.data
print(dask_data)
