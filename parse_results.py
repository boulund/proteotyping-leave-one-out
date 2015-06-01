#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Pretty print leave-one-out results

from sys import argv, exit
import os


def parse_results(filename):
    """ Parse hits to Escherichia coli (562) and total hits from proteotyping output.
    """
    with open(filename) as f:
        lines = f.readlines()[0:15]
        for l in lines:
            if "562" in l:
                correct = l.split()[1]
            if "Total" in l:
                total = l.split()[1]
        sample = filename.split("/", 1)[1].split(".")[0]
        if "blacklist" in filename:
            blacklist = "blacklisted"
        else:
            blacklist = "not blacklisted"
        return sample, blacklist, total, correct


results = []
for fn in os.listdir("results"):
    results.append(parse_results("results/"+fn))

for res in sorted(results, key=lambda x: x[0]):
    print "{:<50} {:>15} {:>10} {:>10}".format(*res)

