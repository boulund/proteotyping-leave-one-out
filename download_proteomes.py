#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Download proteomes from Uniprot


import requests

baseurl = "http://www.uniprot.org/uniprot/"
outputdir = "proteomes/"
organism_IDs = [a.split("	")[2] for a in open("proteomes.txt") if not a.split("	")[2].startswith("Organism")]
for organism in organism_IDs:
    payload = {"query": "organism:{}".format(organism),
               "format": "fasta",
               "include": "yes"}
    r = requests.get(baseurl, params=payload)
    with open(outputdir+organism+".fasta", "w") as f:
        f.write(r.text)
    print "Downloaded proteome for", organism


