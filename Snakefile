
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


def set_sample_batches(sample_file, size_batches=4):
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


rule bamlist:
    input:
        bamlistfile
    output:
        "{batch}/bamlist.txt"
    run:
        with open(output[0], "w") as fout:
            for bam in sample_batches[wildcards.batch]:
                fout.write("%s\n" % bam)


rule runmanta:
    input:
        reference = reference,
        bamlist = "{batch}/bamlist.txt",
        regions = "regions.bed.gz"
    output:
        "{batch}/rundir/workflow.exitcode.txt",
        "{batch}/rundir/results/variants/diploidSV.vcf.gz"
    threads:
        get_threads("runmanta")
    params:
        mem = get_mem("runmanta"),
        rundir = "{batch}/rundir"
    log:
        stdout = "logs/{batch}/run.o",
        stderr = "logs/{batch}/run.e"
    shell:
        """
	    python2 ../manta.py -b {input.bamlist} -r {input.reference} \
         --regions {input.regions} -w {params.rundir} \
         --threads {threads} --mem {params.mem} 1>{log.stdout} 2>{log.stderr}
        """
        #fg_sar start manta
        #python2 ../manta.py -b {input.bamlist} -r {input.reference} \
        # --regions {input.regions} -w {params.rundir} \
        # --threads {threads} --mem {params.mem} 1>{log.stdout} 2>{log.stderr}
        #fg_sar stop manta
        #"""


rule cleanmanta:
    input:
        expand("{svtype}_merged.vcf.gz", svtype=svtypes),
        exitcode="{batch}/rundir/workflow.exitcode.txt"
    output:
        tar="{batch}.tar.gz"
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
        expand("{{svtype}}/{batch}.diploidSV.vcf.gz", batch=sample_batches.keys())
    output:
        "{svtype}.tar.gz"
    shell:
        """
        tar cvzf {wildcards.svtype}.tar.gz {wildcards.svtype}
        """


# rule mergevcf:
#     input:
#         expand("{{svtype}}/{batch}.diploidSV.vcf.gz", batch=sample_batches.keys())
#     output:
#         merged="{svtype}_merged.vcf.gz",
#         tbi="{svtype}_merged.vcf.gz.tbi"
#     params:
#         chromosomes = lambda w: get_chromosomes(w, reference)
#     shadow: "shallow"
#     shell:
#         """
#         ls {wildcards.svtype}/batch*.diploidSV.vcf.gz > sample_files.txt
#         svimmer --max_distance 10 --max_size_difference 10 --ids \
#            sample_files.txt {params.chromosomes} > sample_merged.vcf
#         bcftools sort sample_merged.vcf -Oz -o {output.merged}
#         tabix {output.merged}
#         """
