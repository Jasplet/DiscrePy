#! /usr/bin/env python
## This is a short script to handle eigm files and extract the important lag information to create a summary text file.
import sys
import splitwavepy as sw

file = sys.argv[1]

ef = []
with open(file,'r') as reader:
    eigen_files = reader.readlines()
    for eigen in eigen_files:
        e = eigen.strip().split("/n")[-1]
        ef.append(e)
    print(ef)
#
# for file in files
#     sw.load(file)
