#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Sample sequences from a FASTA file 

from read_fasta import read_fasta
from sys import argv, exit, maxint
import argparse

from numpy.random import choice


def parse_args(argv):
    """Parse commandline arguments.
    """

    desc = """Sample sequences from FASTA files with replacement. Fredrik Boulund 2015"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FASTA", 
            help="FASTA file to sample from.")
    parser.add_argument("-n", metavar="N", required=True, type=int,
            help="Number of sequences to sample from FASTA file [%(default)s].")
    parser.add_argument("--maxlength", metavar="M", type=int,
            default=0,
            help="Maximum length of sequences to sample from, 0 means no limit [%(default)s]")
    parser.add_argument("--minlength", metavar="m", type=int,
            default=0,
            help="Minimum length of sequences to sample from, 0 means no limit [%(default)s].")
    parser.add_argument("--replace", dest="replace", action="store_true",
            default=False,
            help="Sample sequences with replacement [%(default)s].")
    parser.add_argument("-o", "--outfile", metavar="FILE", dest="outfile",
            default="",
            help="Write output to FILE instead of STDOUT.")

    if len(argv)<2:
        parser.print_help()
        exit()
    
    options = parser.parse_args()
    return options


def sample_fasta(fastafile, n, replacement=False, outfile="", maxlength=0, minlength=0):
    """Sample sequences from FASTA.

    Will write to STDOUT if outfile evaluates to False.
    """

    seqs = []
    for header, seq in read_fasta(fastafile):
        seqlen = len(seq)
        if not maxlength:
            maxlength = maxint
        if seqlen >= minlength and seqlen <= maxlength:
            seqs.append((header,seq))

    chosen_samples = choice(len(seqs), size=n, replace=replacement)
    if outfile:
        with open(outfile, 'w') as f:
            for index in chosen_samples:
                header, seq = seqs[index]
                f.write(">"+header+"\n")
                f.write(seq+"\n")
    else:
        for index in chosen_samples:
            header, seq = seqs[index]
            print ">"+header
            print seq


if __name__ == "__main__":
    options = parse_args(argv)
    sample_fasta(options.FASTA, 
                 options.N,
                 options.replace,
                 options.outfile, 
                 options.maxlength,
                 options.minlength)
