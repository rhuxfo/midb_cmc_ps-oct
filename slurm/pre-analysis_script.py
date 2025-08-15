#!/usr/bin/env python3
import sys
import csv
import datetime
import argparse

# Inputs
parser = argparse.ArgumentParser()
parser.add_argument('--csvfile', help='Provide name of the relevant csv file containing tile sizes here. e.g. Moe.csv')
parser.add_argument('--slicenum', help='Please provide the slice number to be analyzed here e.g. 1')
parser.add_argument('--enface_vs_3dtile', help='Please provide "enface" or "3dtile" representing which you would like to analyze.')
parser.add_argument('--enface_vs_3dtile', help='Please provide "enface" or "3dtile" representing which you would like to analyze.')
parser.add_argument('--3dtile_num', help='Please provide the number of the 3d tile you wish to analyze. Only provide if 3dtile is chosen for --enface_vs_3dtile.')
args = parser.parse_args()

# Sanity checks on inputs
if not args.enface_vs_3dtile:
    raise ValueError("Please provide 'enface' or '3dtile'")

if args.enface_vs_3dtile not in ['3dtile','enface']:
    raise ValueError("--enface_vs_3dtile option can only be 'enface' or '3dtile'.")

# Copy the appropriate raw data
with open(args.csvfile) as csvfile:
        slice_reader = csv.reader(csvfile)
        for row in slice_reader:
            if row[1] != args.slicenum:
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
        if args.enface_vs_3dtile == 'enface':
            filename.write(f"""
    P.dir = '/tmp/cmc-s3-bucket/';
    P.Sdir ='/scratch.local/PSOCT/';
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
    P.depthstart = 10;
    P.depthcut = 200;
    P.NoiseCut = 55;
    P.disper = 1;
    P.wind = 1;
    P.BGremoval = 1;
    P.Flect = 1;
    P.Retar = 1;
    P.Cr = 1;
    P.Orio = 0;
    P.AbOrio = 1;
    P.En = 1;
    P.StitchOnly = 0;
    P.Flip = 0;
    P.TCsv = 0;
    P.Ensv = 1;
    P.Stsv = 1;
    P.img = 0;
    Status = PMSDOCT_2024_FCN(P);
    """)

        elif args.enface_vs_3dtile == '3dtile':
            filename.write(f"""
    P.dir = '/tmp/cmc-s3-bucket/';
    P.Sdir ='/scratch.local/PSOCT/';
    P.autofolder=1;
    P.DCf1 = '/scratch.local/ComTom_W_Ch1_shifted.dat';
    P.DCf2 = '/scratch.local/ComTom_W_Ch2_shifted.dat';
    P.Slices = {slice_num}:{slice_num};
    P.tiles = {args.3dtile_num};
    P.buffers = 1:100;
    P.baseN = 'Slice_';
    P.tileN = '_Tile_';
    P.XTiles = {tile_x};
    P.YTiles = {tile_y};
    P.overlap = 10;
    P.depthstart = 10;
    P.depthcut = 200;
    P.NoiseCut = 55;
    P.disper = 1;
    P.wind = 1;
    P.BGremoval = 1;
    P.Flect = 1;
    P.Retar = 1;
    P.Cr = 1;
    P.Orio = 0;
    P.AbOrio = 1;
    P.En = 0;
    P.StitchOnly = 0;
    P.Flip = 1;
    P.TCsv = 1;
    P.Ensv = 0;
    P.Stsv = 0;
    P.img = 0;
    Status = PMSDOCT_2024_FCN(P);
    """)
       
        filename.close()

