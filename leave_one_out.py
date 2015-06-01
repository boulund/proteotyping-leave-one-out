#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Proteotyping leave-one-out analysis

from sys import argv, exit
from collections import defaultdict
from numpy.random import choice, seed
from subprocess import call
from textwrap import TextWrapper
from fnmatch import fnmatch
from multiprocessing import Pool, current_process
import os 


from read_fasta import read_fasta
from sample_fasta import sample_fasta


def find_files(directory, pattern, folder_pattern):
    """ Generator that yields files recursively using glob pattern.
    """

    for root, subfolders, files in os.walk(directory, followlinks=True):
        if folder_pattern in root:
            genome_list = []
            for basename in files:
                if fnmatch(basename, pattern):
                    filename = os.path.join(root, basename)
                    genome_list.append(filename)
            yield (root, genome_list)


def read_proteome_genome_paths(protdir, genomedir):
    """ Construct a dictionary with proteome:genome mappings.

    Disregards any proteomes without corresponding genome sequences.
    """
    paths = defaultdict(list)
    accnos = defaultdict(list)
    proteome_paths = [p for p in find_files(protdir, "*.faa", "Escherichia_coli")]
    genome_paths = [p for p in find_files(genomedir, "*.fna", "Escherichia_coli")]

    for proteomedir, proteomelist in proteome_paths:
        strain = proteomedir.rsplit("/", 1)[1]
        for genomedir, genomelist in genome_paths:
            if strain in genomedir:
                paths[strain].extend(proteomelist)
                for filename in genomelist:
                    with open(filename) as genome_file:
                        header_line = genome_file.readline()
                        accno = header_line.split("ref|")[1].split("|", 1)[0]
                        accnos[strain].append(accno)
    return paths, accnos


def digest_proteome(proteome_path):
    """ Digest a proteome using EMBOSS digest.
    """
    
    digest_call = "digest -menu 1 -mono N -outfile {output} -rformat2 srs {input}"
    outfile = proteome_path+".srs"
    retcode = call(digest_call.format(input=proteome_path, output=outfile), shell=True)
    return outfile


def convert_srs_2_fasta(srsfile, fastafile):
    """ Convert EMBOSS digest's SRS2 format to FASTA.
    """
    
    with open(srsfile) as srs:
        with open(fastafile, "w") as fasta:
            tw = TextWrapper()
            for line in srs:
                if line.startswith("# Sequence:"):
                    header = line.split()[2]
                elif line.startswith("Feature:"):
                    feature = line.split()[1]
                elif line.startswith("Length:"):
                    seqlen = line.split()[1]
                    fasta.write(">{h}{f}_{length}\n".format(h=header, f=feature, length=seqlen))
                elif line.startswith("Sequence:"):
                    seq = line.split()[1]
                    for seqline in tw.wrap(seq):
                        fasta.write(seqline+"\n")


def run_workflow(strain):
    """ Run the complete workflow on a single strain.
    """
    workername = current_process().name

    print workername+": Concatenating all proteome files for", strain
    call("cat {} > ./concat_proteomes/{}.fasta".format(" ".join(paths[strain]), strain), shell=True)
    concat_proteome = "./concat_proteomes/{}.fasta".format(strain)

    print workername+": Digesting", concat_proteome
    outfile = digest_proteome(concat_proteome)
    
    print workername+": Converting SRS to FASTA..."
    fastafile = outfile+".fasta"
    convert_srs_2_fasta(outfile, fastafile)

    sampled_fasta = fastafile+".sampled"
    num_fragments = 10000
    print workername+": Sampling", num_fragments, "fragments from", fastafile
    sample_fasta(fastafile, 
            num_fragments, 
            replacement=False, 
            outfile=sampled_fasta, 
            minlength=6, 
            maxlength=45)
    
    print workername+": Running BLAT on", sampled_fasta
    mappings = sampled_fasta+".blast8"
    blat_call = "blat -t=dnax -q=prot -tileSize=5 -stepSize=5 -minScore=10 -out=blast8 -minIdentity=85 {database} {query} {output}"
    db = "/shared/genomes/NCBI/bacterial/bacterial_genomes.fasta"
    call(blat_call.format(database=db, query=sampled_fasta, output=mappings), shell=True)

    print workername+": Running proteotyping without blacklisting on", mappings
    resultsfile = "./results/"+os.path.basename(mappings)+".results"
    proteotyping_call = "proteotyping.py --accno_annotation_pickle accno_annotation.pkl --taxtree_pickle taxtree.pkl --gene_info /shared/db/NCBI/gene/gene_info --blacklist_accnos blacklist.txt --logfile {logfile} --output {output} {input}"
    call(proteotyping_call.format(input=mappings, 
                                  output=resultsfile, 
                                  logfile=strain+".log"), shell=True)

    print workername+": Running proteotyping WITH blacklisting on", mappings
    resultsfile = "./results/"+os.path.basename(mappings)+".blacklist.results"
    proteotyping_call = "proteotyping.py --accno_annotation_pickle accno_annotation.pkl --taxtree_pickle taxtree.pkl --gene_info /shared/db/NCBI/gene/gene_info --leave-out {blacklist} --blacklist_accnos blacklist.txt --logfile {logfile} --output {output} {input}"
    call(proteotyping_call.format(input=mappings, 
                                  output=resultsfile, 
                                  logfile=strain+"blacklisted.log",
                                  blacklist=",".join(accnos[strain])), shell=True)




if __name__ == "__main__":
    seed(1337) # Seed the RNG for replication purposes

    proteome_dir = "/shared/genomes/NCBI/bacterial/20150525_proteomes/"
    genome_dir = "/shared/genomes/NCBI/bacterial/20150210/"
    paths, accnos = read_proteome_genome_paths(proteome_dir, genome_dir)
    print "Loaded paths to", len(paths), "proteomes and", len(accnos), "genomes from disk."

    print "Starting worker pool with 6 processess..."
    pool = Pool(processes=8)
    #strains = choice(paths.keys(), size=12, replace=True)
    strains = paths.keys()
    pool.map(run_workflow, strains)
    print "Finished."
