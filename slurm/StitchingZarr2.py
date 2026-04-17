import numpy as np
import zarr
import scipy.io as sio
import h5py
#from matplotlib import pyplot as plt
#import fsspec
import argparse
from rclone_python import rclone

parser = argparse.ArgumentParser()
parser.add_argument("filename",help='Full save path of Tile to convert')
parser.add_argument("Contrast",help='contrast of tile for naming')
parser.add_argument("Tile_name",help='Tile_cross, Tile_R, etc')
args = parser.parse_args()

load_path = args.filename+'.mat'
#fs = fsspec.filesystem('s3',)
save_path = args.filename+'.zarr'

mat_contents = h5py.File(load_path,'r')

tile_data = mat_contents.get(args.Tile_name)
#tile_data = np.swapaxes(mat_contents.get('Tile_cross'), 0,1)
tile_np = np.array(tile_data)

#plt.imshow(tile_np[::-1,:,65], cmap='gray', vmin=40); plt.show()
#create zarr file and then save data to that file
tile_zarr = zarr.open(save_path, mode='w', shape=tile_data.shape, chunks=tile_data.shape, dtype=tile_data.dtype)
tile_zarr[:,:,:] = tile_np[:,:,:]
rclone.copy(save_path,"AWScmc:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/${SUBJECT_NAME}/PS-OCT/3DTiles/{args.Contrast}/", args=["--s3-no-check-bucket"])
