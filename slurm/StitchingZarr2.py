import numpy as np
import zarr
from rclone_python import rclone

save_path = filename+'.zarr'

tile_np = np.array(tile_data)

#create zarr file and then save data to that file
tile_zarr = zarr.open(save_path, mode='w', shape=tile_data.shape, chunks=tile_data.shape, dtype=tile_data.dtype)
tile_zarr[:,:,:] = tile_np[:,:,:]
rclone.copy(save_path,"AWScmc:cmc-msi-accesspoint-2-254319122668/CMC/Derivatives/Moe/PS-OCT/3DTiles/{Contrast}/", args=["--s3-no-check-bucket"])
