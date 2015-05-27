#!/usr/bin/env python
# Fredrik Boulund 2015
# Convert output from EMBOSS digest in SRS format to FASTA.

from sys import argv, exit
from textwrap import TextWrapper

if len(argv)<2:
    print "usage: convert_srs_to_fasta.py FILE"
    exit()

with open(argv[1]) as srs:
    tw = TextWrapper()
    for line in srs:
        if line.startswith("# Sequence:"):
            header = line.split()[2]
        elif line.startswith("Feature:"):
            feature = line.split()[1]
        elif line.startswith("Length:"):
            seqlen = line.split()[1]
            print ">{h}{f}_{length}".format(h=header, f=feature, length=seqlen)
        elif line.startswith("Sequence:"):
            seq = line.split()[1]
            for seqline in tw.wrap(seq):
                print seqline
