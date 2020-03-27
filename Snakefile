
import os
import shutil
from collections import defaultdict

DEFAULT_THREADS = 28


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


def get_chromosomes(wildcards, reference, chromregex='(chr)?([1-9][0-9]?|[XY])'):
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
    return [bam for bam in sample_batches[wildcards.batch]


if not tool_exists("tabix"):
    print("tabix is not available")
    exit(1)


if not tool_exists("bgzip"):
    print("bgzip is not available")
    exit(1)


if not os.path.isfile("%s.fai" % config['reference']):
    print("reference genome must be fai indexed")


bamlistfile = os.path.abspath(config['samples'])
sample_batches = set_sample_batches(bamlistfile)
reference = os.path.abspath(config['reference'])

svtypes = ["DEL", "INS"]

workdir: config["workdir"]

localrules: bamlist, cleanmanta, regions

# Wildcard constraints
wildcard_constraints:
    batch="|".join(sample_batches.keys()),
    svtype="|".join(svtypes),


rule all:
    input:
        expand("{batch}.tar.gz", batch=sample_batches.keys()),
        expand("{svtype}.tar.gz", svtype=svtypes)


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
        " tabix {output} "

rule configmanta:
    input:
        reference=reference,
        bams=get_bam_files,
        regions="regions.bed.gz"
    output:
        "{batch}/rundir/workflow.exitcode.txt",
        "{batch}/rundir/results/variants/diploidSV.vcf.gz"
    threads:
        1
    params:
        rundir = "{batch}/rundir"
    log:
        stdout = "logs/{batch}/run.o",
        stderr = "logs/{batch}/run.e"
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


rule runmanta:
    input:
        mantarunner = "{batch}rundir/runWorkflow.py"
    output:
        "{batch}/rundir/workflow.exitcode.txt",
        "{batch}/rundir/results/variants/diploidSV.vcf.gz"
    threads:
        get_threads("runmanta")
    params:
        mem=get_mem("runmanta")
    log:
        stdout = "logs/run.o",
        stderr = "logs/run.e"
    shell:
        "python2 {input} -j {threads} -g {params.mem} --quiet "
        " 1>{log.stdout} 2>{log.stderr} "


rule cleanmanta:
    input:
        expand("{svtype}/{{batch}}.diploidSV.vcf.gz", svtype=svtypes),
        exitcode="{batch}/rundir/workflow.exitcode.txt"
    output:
        tar="{batch}.tar.gz"
    priority: 10
    shell:
        """
        exitcode=`cat {wildcards.batch}/rundir/workflow.exitcode.txt`
        if [ $exitcode == 0 ];
        then
            rm -fr {wildcards.batch}/rundir/workspace
            tar cvzf {wildcards.batch}.tar.gz {wildcards.batch}
            rm -fr {wildcards.batch}
        fi
        """

rule splitvcf:
    input:
        exitcode="{batch}/rundir/workflow.exitcode.txt",
        diploidsv="{batch}/rundir/results/variants/diploidSV.vcf.gz"
    output:
        "{svtype}/{batch}.diploidSV.vcf.gz"
    priority: 20
    shell:
        """
        exitcode=`cat {wildcards.batch}/rundir/workflow.exitcode.txt`
        if [ $exitcode == 0 ];
        then
           bcftools view -i 'INFO/SVTYPE=="{wildcards.svtype}"' {input.diploidsv} -Oz -o {output}
          tabix {output}
        fi
        """

rule tarsvtype:
    input:
        expand("{{svtype}}/{batch}.diploidSV.vcf.gz", batch=sample_batches.keys()),
        expand("{batch}.tar.gz", batch=sample_batches.keys())
    output:
        "{svtype}.tar.gz"
    shell:
        """
        tar cvzf {wildcards.svtype}.tar.gz {wildcards.svtype}
        rm -fr {wildcards.svtype}
        """
