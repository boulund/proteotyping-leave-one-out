#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Pretty print leave-one-out results

from __future__ import division
from sys import argv, exit
import os


def parse_results(filename):
    """ Parse hits to Escherichia coli (562) and total hits from proteotyping output.
    """
    with open(filename) as f:
        lines = f.readlines()[0:15]
        for l in lines:
            if "562" in l:
                correct = int(l.split()[1])
            if "Total" in l:
                total = int(l.split()[1])
        sample = filename.split("/", 1)[1].split(".")[0]
        if "blacklist" in filename:
            blacklist = True
        else:
            blacklist = False 
        return sample, blacklist, total, correct


results = []
for fn in os.listdir("results"):
    results.append(parse_results("results/"+fn))

blacklisted = {a[0]:a[2:] for a in results if a[1]}
non_blacklisted = {a[0]:a[2:] for a in results if not a[1]}
completed = set(blacklisted.keys()).intersection(set(non_blacklisted.keys()))

print "{:<50} {:<10}".format("Strain", "Sensitivity difference")
outputfmt = "{:<50} {:<1.9f}" 
for strain in sorted(completed):
    blacklisted_tot_disc = blacklisted[strain][0]
    non_blacklisted_tot_disc = non_blacklisted[strain][0]
    blacklisted_correct_disc = blacklisted[strain][1]
    non_blacklisted_correct_disc = non_blacklisted[strain][1]

    sensitivity_diff = non_blacklisted_correct_disc/non_blacklisted_tot_disc - blacklisted_correct_disc/blacklisted_tot_disc
    print outputfmt.format(strain, sensitivity_diff)

blacklisted_tot_mean = sum([a[0] for a in blacklisted.values()])/len(completed)
blacklisted_cor_mean = sum([a[1] for a in blacklisted.values()])/len(completed)
non_blacklised_tot_mean = sum([a[0] for a in non_blacklisted.values()])/len(completed)
non_blacklised_cor_mean = sum([a[1] for a in non_blacklisted.values()])/len(completed)

print "Average number of discriminative fragments blacklisted"
print blacklisted_tot_mean
print blacklisted_cor_mean 
print "Average number of discriminative fragments non-blacklisted"
print non_blacklised_tot_mean 
print non_blacklised_cor_mean 


# Print raw data
print "{:<50} {:>15} {:>10} {:>10}".format("Strain", "Blacklisted", "Tot disc.", "Cor. disc.")
for result in sorted(results, key=lambda x: x[0]):
    print "{:<50} {:>15} {:>10} {:>10}".format(*result)
