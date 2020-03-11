
# Thursday december 5 2019
# Copying 10 bams on genotoul for testing purpose

module load bioinfo/samtools-1.9
module load bioinfo/bcftools-1.9
module load bioinfo/snakemake-4.8.0

manta=/home/faraut/dynawork/CNVPipeline/work/vargoats/cobalt2.0/manta-1.6.0.centos6_x86_64/bin

here=`pwd`
export PATH=$manta:$here:$PATH

snakemake --configfile config.yaml \
          --cluster-config cluster.yaml \
          --cluster "srun --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1" \
          --jobs 4 -p -n


