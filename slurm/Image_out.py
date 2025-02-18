#!/usr/bin/env python3
import sys
import csv
import datetime
import argparse

# Inputs
parser = argparse.ArgumentParser()
parser.add_argument('--slicenum', help='Please provide the slice number to be analyzed here e.g. 1')
args = parser.parse_args()

# Write the wrapper function
with open(f'/tmp/slice_{slice_num}_in.m', 'w') as filename:
       filename.write(f"""
D = '/tmp/cmc-s3-bucket/';
Sv ='/tmp/';
Slices = {slice_num}:{slice_num};

Img = gifStack(D,Sv,Slices);
""")
       filename.close()
