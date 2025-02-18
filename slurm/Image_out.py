#!/usr/bin/env python3
import sys
import argparse

# Inputs
parser = argparse.ArgumentParser()
parser.add_argument('--slicenum', help='Please provide the slice number to be analyzed here e.g. 1')
args = parser.parse_args()

slice_num=args.slicenum
# Write the wrapper function
with open(f'/tmp/slice_{slice_num}_in.m', 'w') as filename:
       filename.write(f"""
D = '/tmp/cmc-s3-bucket/';
Sv ='/tmp/';
Slices = 1:{slice_num};

Img = gifStack(D,Sv,Slices);
""")
       filename.close()
