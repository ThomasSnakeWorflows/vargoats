# vargoats

### CNVpipeline for Vargoats

**Environment and runing on genotoul**

```bash
module load bioinfo/samtools-1.9
module load bioinfo/bcftools-1.9
module load bioinfo/snakemake-4.8.0

manta=/home/faraut/dynawork/CNVPipeline/work/vargoats/cobalt2.0/manta-1.6.0.centos6_x86_64/bin
here=`pwd`
export PATH=$manta:$here:$PATH
```

Running the CNVPipeline

```bash
snakemake --configfile config.yaml \
          --cluster-config cluster.yaml \
          --drmaa " --mem={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" \
          --jobs 4 -p -n
```

**Environment and runing on cobablt**

```bash
module load extenv/fg
module load fgtools
fg_project fg0047
module load bcftools/1.6
module load tabix/0.2.6
module load snakemake/4.8.0
```

```bash
snakemake --configfile config.yaml \
          --cluster-config cluster.yaml \
          --cluster "ccc_msub -A fg0047 -q broadwell -Q normal -n 1 -c {threads} -T {cluster.time} -o manta_job_%j.out -e manta_job_%j.err " \
          --jobs 95 -p  > snake.log
```
