# !/usr/bin/env python

# script generates bed file that scpecifies regions with coverage below defined level
# usage: samtools depth <file.bam> | bed_from_areas_covered_above_N.py N

import sys
import re
import pandas as pd
import subprocess
import fileinput

if(len(sys.argv) < 1):
	print("usage: samtools depth <file.bam> | bed_from_areas_covered_above_N.py N")
	exit(0)

depth_in = open(sys.argv[1], "r+")

N = int(sys.argv[2])

while(1): # read only first few lines, until you get first position covered below required level
	line = depth_in.readline()
	if not line:
		break
	x = line.rstrip().split("\t")
	chrom=x[0]
	start_pos = int(x[1])
	last_pos = start_pos
	num = int(x[2])
	if(num <= N): # read until first position below required depth
		break

while (1):
	line = depth_in.readline()
	if not line:
		break
	x = line.rstrip().split("\t")
	if(int(x[2]) > N): # depth is more than we want
		continue
	this_pos = int(x[1])
	if(last_pos != this_pos-1):
		print(chrom+"\t"+str(start_pos)+"\t"+str(last_pos))
		chrom = x[0]
		start_pos = int(x[1])
		last_pos = int(x[1])
	last_pos = this_pos

