#! /usr/bin/python3

from multiprocessing import Pool
from collections import Counter
import pysam
import sys
import csv

samfile=pysam.AlignmentFile(sys.argv[1], "rb")
idxfile=sys.argv[2]
outDir=sys.argv[3]
nproc=sys.argv[4]
cpu=int(nproc)

idxstats=open(idxfile, 'r')
chrs = idxstats.read().split()
idxstats.close()

def splitBigBam(chr):
	Itr = samfile.fetch(str(chr), multiple_iterators=True)
	splitBam = pysam.AlignmentFile(outDir + "/" + str(chr) + ".bam", "wb", template=samfile)
	for read in Itr:
		spl = sorted(read.query_name.split("_"), key=len)
		bc = spl[0]
		rd = spl[1]
		read.set_tag("XU", bc)
		read.query_name = rd
		splitBam.write(read)
	splitBam.close()

pool = Pool(processes=cpu)
toy_out = pool.map(splitBigBam, chrs)
pool.close()

