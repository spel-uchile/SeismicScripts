#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from os import listdir
from os.path import isfile, join

parser = argparse.ArgumentParser(description='Obspy wrapper: Apply \"change extension\" operation for infolder')
parser.add_argument('directory', help='directory to use', action='store')
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--extin', action='store', help='extin to filter', required=True)
parser.add_argument('--extout', action='store', help='extout to filter', required=True)
#parser.add_argument('--outfile', action='store', help='output file')
#parser.add_argument('--indiv', action='store_true', help='output individual plots for every file')
args = parser.parse_args()

# 1) Make sure user inputs are correct
extin = args.extin
print(extin)
extout = args.extout
print(extout)
# Convert to real (no symlink) and full path
infolder_path = args.directory
infolder_path = os.path.normcase(infolder_path)
infolder_path = os.path.normpath(infolder_path)
infolder_path = os.path.realpath(infolder_path)
print(infolder_path)
# Get all files in folder that contain "extin" and "extout"
onlyfiles = [f for f in listdir(infolder_path) if isfile(join(infolder_path, f))]
onlyfiles.sort()
sel_files = []
for file_i in onlyfiles:
    filename, file_extension = os.path.splitext(file_i)
    file_extension = file_extension[1:]
    #print(file_extension)
    if extin == file_extension:
        sel_files.append(file_i)
sel_files.sort()
infolder_path = sel_files
print(infolder_path)
#exit(0)

# 2) Change names
for file_i in infolder_path:
    filename, file_extension = os.path.splitext(file_i)
    new_file_i = filename + "." + extout
    os.rename(file_i, new_file_i)