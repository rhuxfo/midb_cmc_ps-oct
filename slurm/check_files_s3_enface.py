#!/usr/bin/env python3
import sys
import os
# fetch the contents of that directory
this_list = os.listdir(f'/tmp/cmc-s3-bucket/')
# Check that all files for that slice exist                        
for slice_num in range(1, 328):
  filename = f'Slice_{slice_num}_EnR.mat'

# if file not in dir, then print to outfile
  if filename not in this_list:
    print(filename, flush = True)
