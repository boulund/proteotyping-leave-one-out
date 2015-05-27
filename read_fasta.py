#!/usr/bin/env python
# Fredrik Boulund 2015
# Yield sequences from FASTA file

def read_fasta(filename, keep_formatting=True):
    """Read sequence entries from FASTA file
    NOTE: This is a generator, it yields after each completed sequence.
    Usage example:
    for header, seq in read_fasta(filename):
        print ">"+header
        print seq
    """

    with open(filename) as fasta:
        line = fasta.readline().rstrip()
        if not line.startswith(">"):
            raise IOError("Not FASTA format? First line didn't start with '>'")
        if keep_formatting:
            sep = "\n"
        else:
            sep = ""
        first = True
        seq = []
        header = ""
        while fasta:
            if line == "": #EOF
                yield header, sep.join(seq)
                break
            elif line.startswith(">") and not first:
                yield header, sep.join(seq)
                header = line[1:]
                seq = []
            elif line.startswith(">") and first:
                header = line[1:]
                first = False
            else:
                seq.append(line)
            line = fasta.readline().rstrip()
