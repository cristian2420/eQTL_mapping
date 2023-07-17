#!/bin/bash
#SBATCH --job-name=eQTL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10g
#SBATCH --time=20:00:00
#SBATCH --output=pcs.out
#SBATCH --error=pcs.err

cd /path/to/pipeline/

WORKDIR=/path/to/working/directory/
log_path=${WORKDIR}/logs/

start=`date +%s`
date
echo 'Running snakemake' 

snakemake --jobs 100 --latency-wait 60 --snakefile Snakefile --configfile snake_conf.yaml --cluster-config cluster_peer.json --cluster "sbatch --time={cluster.walltime} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem={cluster.memory} -e {rule}.{jobid}.{wildcards}.err -o {rule}.{jobid}.{wildcards}.out --export ALL --parsable" --stats $log_path/snakemake.stats >& $log_path/snakemake.log  --rerun-incomplete --use-conda 

end=`date +%s`
runtime=$((end-start))
echo 'Running time ' ${runtime} ' seconds'

