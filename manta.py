#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import argparse
import re
import subprocess
import logging

#manta_bin = "/home/faraut/dynawork/CNVPipeline/work/vargoats/cobalt2.0/manta-1.6.0.centos6_x86_64/bin"
manta_bin = "/ccc/cont007/home/fg0047/farautth/work/vargoats/manta-1.6.0.centos6_x86_64/bin"
configManta = manta_bin + "/configManta.py"

def get_bamfiles(sample_file):
    with open(sample_file) as fin:
        bams = [os.path.abspath(line.rstrip()) for line in fin]
    return bams


def extract_valid_regions(reference, chromregex='(chr)?([1-9][0-9]?|[XY])'):
    fai = reference + ".fai"
    region_file = "regions.bed"
    p = re.compile(chromregex)
    with open(fai) as fin, open("regions.bed", "w") as fout:
        for line in fin:
            fields = line.split("\t")
            if p.match(fields[0]):
                fout.write("%s\t%d\t%d\n" % (fields[0], 0, int(fields[1])))
    command = "bgzip -f %s; tabix %s.gz" % (region_file, region_file)
    subprocess.call(command, shell=True)
    return "%s.gz" % region_file


def config_manta(bamlist_file, reference, regions, workdir):
    """ Run manta configuration script """
    bams = get_bamfiles(bamlist_file)
    command = configManta
    for bam in bams:
        command += " --bam %s " % bam
    command += " --referenceFasta %s " % reference
    command += " --callRegions %s " % regions
    command += " --runDir %s " % workdir
    try:
        print(command)
        subprocess.call(command, shell=True)
    except subprocess.CalledProcessError:
        print("Manta configuration failed")


def run_manta(workdir, threads=1, memory=12):
    manta_workflow = "%s/runWorkflow.py" % workdir
    manta_workflow += " --jobs %d " % threads
    manta_workflow += " -g %d " % memory
    manta_workflow += " --quiet "
    subprocess.call(manta_workflow, shell=True)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="manta.py",
        description="manta wrapper ")
    parser.add_argument("-b", "--bamlist", required=True,
                        help="A file with a list of BAM files")
    parser.add_argument("-r", "--reference", required=True,
                        help="The reference genome (associated fai required)")
    parser.add_argument("--regions",
                        help="Restrict manta detection to the specified "
                             "regions (compressed tabix bed file)")
    parser.add_argument("-w", "--work", required=True,
                        help="working dir")
    parser.add_argument("-t", "--threads",
                        default=1, type=int,
                        help="number of threads")
    parser.add_argument("-m", "--mem",
                        default=12, type=int,
                        help="memory available")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="increase verbosity")
    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()
    bamlist = args.bamlist
    reference = args.reference
    regions = args.regions
    workdir = args.work
    threads = args.threads
    memory =  args.mem

    logging.basicConfig(format="manta log: %(message)s",level=logging.INFO)

    logging.info('Configuring manta')
    if regions is None:
        regions = extract_valid_regions(reference)
    config_manta(bamlist, reference, regions, workdir)

    logging.info('Running manta')
    run_manta(workdir, threads, memory)


if __name__ == "__main__":
    sys.exit(main())
