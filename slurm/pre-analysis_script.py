#!/usr/bin/env python3
import sys
import csv
import datetime

# Copy the appropriate raw data
with open('Kumquat_Right_Hemisphere_Slices_Sheet2.csv') as csvfile:
        slice_reader = csv.reader(csvfile)
        for row in slice_reader:
            if row[1] != sys.argv[1]:
                    continue
            else:
                print(row)
                date = row[0]
                dir_date = datetime.datetime.strptime(date, "%m/%d/%Y").strftime("%m%d%Y")
                print(dir_date)

                slice_num, tile_x, tile_y = row[1:4]
                max_tile = int(tile_x)*int(tile_y)

# Write the wrapper function
with open(f'/tmp/slice_{slice_num}_wrapper.m', 'w') as filename:
       filename.write(f"""
P.dir = '/tmp/cmc-s3-bucket/';
P.Sdir ='/tmp/';
P.autofolder=1;
P.DCf1 = '/scratch.local/ComTom_W_Ch1_shifted.dat';
P.DCf2 = '/scratch.local/ComTom_W_Ch2_shifted.dat';
P.Slices = {slice_num}:{slice_num};
P.tiles = 1:{max_tile};
P.buffers = 1:100;
P.baseN = 'Slice_';
P.tileN = '_Tile_';
P.XTiles = {tile_x};
P.YTiles = {tile_y};
P.overlap = 10;
P.depthstart = 20;
P.depthcut = 300;
P.Ch_dB_limit = 55;
P.Rline = 264;
P.disper = 1;
P.wind = 1;
P.BGremoval = 1;
P.CDP = 0;
P.CH12 = 1;
P.Flect = 1;
P.Retar = 0;
P.Cr = 0;
P.Orio = 0;
P.AbOrio = 0;
P.En = 1;
P.StitchOnly = 0;
P.TCsv = 0;
P.Ensv = 1;
P.Stsv = 1;
Status = PMSDOCT_2024_FCN(P);
""")
       filename.close()

