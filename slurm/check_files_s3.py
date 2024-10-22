#!/usr/bin/env python3
import argparse
import csv
import datetime
import sys
import os

# Inputs
parser = argparse.ArgumentParser()
parser.add_argument('--csvfile', help='Provide name of the relevant csv file containing tile sizes here. e.g. Moe.csv')
args = parser.parse_args()

# Open the csv
with open(args.csvfile) as csvfile:
        slice_reader = csv.reader(csvfile)
        count=0
        # Go slice by slice
        for row in slice_reader:
                # extract the date
                date = row[0]

                try:
                    dir_date = datetime.datetime.strptime(date, "%m/%d/%Y").strftime("%m%d%Y")
                except ValueError:
                    continue

                # fetch the contents of that directory
                this_date_list = os.listdir(f'/tmp/cmc-s3-bucket/{dir_date}')

                # Check that all files for that slice exist
                # Row[2] = # of x tiles, row[3] = # of Y tiles
                for tile_num in range(1, (int(row[2])*int(row[3])) + 1):
                        for buffer_num in range(1, 101):
                            filename = f'Slice_{row[1]}_Tile_{tile_num}_840_{buffer_num}.dat'

                            # if file not in dir, then print to outfile
                            if filename not in this_date_list:
                                print(filename, flush = True)

