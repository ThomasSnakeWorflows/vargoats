#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import sys

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pysam import VariantFile
from Bio.SeqUtils import lcc


def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


def svlen(record):
    return max(len(record.alts[0]), len(record.ref)) - 1


def variant_complexity(record):
    alt = record.alts[0]
    ref = record.ref
    if len(ref) > len(alt):
        return lcc.lcc_simp(ref)
    else:
        return lcc.lcc_simp(alt)


def valid_record(record):
    if ("<" not in record.alts[0] and
            len(record.alts) == 1 and
            record.alts[0][0] == record.ref[0]):
        return True
    return False


def valid_chrom(chrom, valid_chromosomes):
    if not valid_chromosomes:
        return True
    elif chrom in valid_chromosomes:
        return True
    return False


def svtype(record):
    if len(record.alts[0]) > len(record.ref):
        return "INS"
    else:
        return "DEL"


def get_valid_chromosomes(chromosomes):
    if not chromosomes:
        return []
    else:
        return set(chromosomes.split(","))


def main(vcffile, outvcf, prefix):
    vcf_in = VariantFile(vcffile)
    vcf_out = VariantFile(outvcf, 'w', header=vcf_in.header)
    for record in vcf_in:
        id = record.id
        record.id = "%s_%s" % (prefix, id)
        vcf_out.write(record)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract SV allele sequences")
    parser.add_argument('-v', '--vcf', required=True,
                        help='the vcf file file')
    parser.add_argument('-o', '--out', required=True,
                        help='the vcf file file')
    parser.add_argument('-p', '--prefix', required=True,
                        help='the output vcf file file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.vcf, args.out, args.prefix)
