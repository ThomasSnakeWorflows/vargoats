
import os
import shutil
from collections import defaultdict

DEFAULT_THREADS = 3


def tool_exists(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    from shutil import which
    return which(name) is not None


def get_threads(rule, default=DEFAULT_THREADS):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "threads" in cluster_config[rule]:
        return cluster_config[rule]["threads"]
    if "default" in cluster_config and "threads" in cluster_config["default"]:
        return cluster_config["default"]["threads"]
    return default


def get_mem(rule):
    threads = get_threads(rule)
    print(threads)
    return threads * 4


def set_sample_batches(sample_file, size_batches=10):
    """
       Construct batches os samples of size size_batches
    """
    samples = defaultdict()
    samples_batch = []
    start_batch_nb = 1
    with open(sample_file) as fin:
        for line in fin:
            bam = os.path.abspath(line.rstrip())
            samples_batch.append(bam)
            if len(samples_batch) == size_batches:
                samples["batch%03d" % start_batch_nb] = samples_batch
                samples_batch = []
                start_batch_nb += 1
            if samples_batch:
                samples["batch%03d" % start_batch_nb] = samples_batch
    return samples


def get_chromosomes_old(wildcards, reference, chromregex='(chr)?([1-9][0-9]?|[XY])'):
    referencefai = reference + ".fai"
    p = re.compile(chromregex)
    chromosomes = []
    with open(referencefai) as fin:
        contigs = [line.split("\t")[0] for line in fin]
    for contig in contigs:
        if p.match(contig):
            chromosomes.append(contig)
    return chromosomes


def get_chromosomes(reference, chromregex='(chr)?([1-9][0-9]?|[XY])'):
    referencefai = reference + ".fai"
    p = re.compile(chromregex)
    chromosomes = []
    with open(referencefai) as fin:
        contigs = [line.split("\t")[0] for line in fin]
    for contig in contigs:
        if p.match(contig):
            chromosomes.append(contig)
    return chromosomes


def get_bam_files(wildcards):
    return [bam for bam in sample_batches[wildcards.batch]]


if not tool_exists("tabix"):
    print("tabix is not available")
    exit(1)


if not tool_exists("bgzip"):
    print("bgzip is not available")
    exit(1)


bamlistfile = os.path.abspath(config['samples'])
sample_batches = set_sample_batches(bamlistfile)
reference = os.path.abspath(config['reference'])

if not os.path.isfile("%s.fai" % config['reference']):
    print("reference genome must be fai indexed")


chromosomes = get_chromosomes(reference, config['chromregex'])

workdir: config["workdir"]

localrules: chromregion, configchrommanta, concat

# Wildcard constraints
wildcard_constraints:
    batch="|".join(sample_batches.keys()),
    chrom="|".join(chromosomes)

#
# rule all:
#     input:
#         expand("{batch}/{chrom}/rundir/workflow.exitcode.txt",
#                batch=sample_batches.keys(),
#                chrom=chromosomes)


rule all:
    input:
        expand("{batch}/diploidSV.vcf.gz.tbi", batch=sample_batches.keys())


rule regions:
    input:
        reference = reference,
        fai = reference+".fai"
    output:
        "regions.bed.gz"
    params:
        chromregex = config['chromregex']
    shell:
        "cat {input.fai} | egrep \'^{params.chromregex}\' "
        " | awk -v OFS='\\t' '{{ print $1,0,$2}}' | bgzip > {output}; "
        " tabix {output.vcf} "


rule chromregion:
    input:
        reference = reference,
        fai = reference+".fai"
    output:
        "regions/{chrom}.bed.gz"
    params:
        chrom = lambda w: w.chrom
    shell:
        "cat {input.fai} | grep -w  \'^{params.chrom}\' "
        " | awk -v OFS='\\t' '{{ print $1,0,$2}}' | bgzip > {output}; "
        " tabix {output} "


rule configchrommanta:
    input:
        "%s.fai" % reference,
        reference=reference,
        bams=get_bam_files,
        regions="regions/{chrom}.bed.gz"
    output:
        "{batch}/{chrom}/rundir/runWorkflow.py"
    threads:
        1
    params:
        rundir = "{batch}/{chrom}/rundir"
    log:
        stdout = "logs/{batch}/{chrom}_config.o",
        stderr = "logs/{batch}/{chrom}_config.e"
    run:
        import subprocess
        command = "configManta.py "
        for bam in input.bams:
            command += "--bam %s " % bam
        command += "--referenceFasta %s " % input.reference
        command += "--runDir %s " % params.rundir
        command += "--callRegions %s " % input.regions
        try:
            print(command)
            subprocess.call(command, shell=True)
        except subprocess.CalledProcessError:
            print("Manta configuration failed")


rule runchrommanta:
    input:
        mantarunner = "{batch}/{chrom}/rundir/runWorkflow.py"
    output:
        "{batch}/{chrom}/rundir/workflow.exitcode.txt",
        "{batch}/{chrom}/rundir/results/variants/diploidSV.vcf.gz"
    threads:
        get_threads("runchrommanta",16)
    params:
        mem=get_mem("runchrommanta")
    log:
        stdout = "logs/{batch}/{chrom}_run.o",
        stderr = "logs/{batch}/{chrom}_run.e"
    shell:
        """
        python2 {input} -j {threads} -g {params.mem} --quiet \
         1>{log.stdout} 2>{log.stderr}
        """

rule concat:
    input:
        expand("{{batch}}/{chrom}/rundir/results/variants/diploidSV.vcf.gz",
                  chrom=chromosomes)
    output:
        "{batch}/diploidSV.vcf.gz.tbi",
        vcf="{batch}/diploidSV.vcf.gz"

    log:
        stdout = "logs/{batch}/concat.o",
        stderr = "logs/{batch}/concat.e"
    shell:
        """
        bcftools concat -Oz -o {output.vcf} {input}
        tabix {output.vcf}
        """




# """
# fg_sar start manta_{wildcards.batch}_{wildcards.chrom}
# python2 {input} -j {threads} -g {params.mem} --quiet \
#  1>{log.stdout} 2>{log.stderr}
# fg_sar stop manta_{wildcards.batch}_{wildcards.chrom}
# """
