# Bulk_RNA-Seq_Analysis

# 3) Instructions to run on Slurm managed HPC
3A. Download version controlled repository
```
git clone https://github.com/kevinboyd76/Bulk_RNA-Seq_Analysis.git
```
3B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
3C. Modify samples and config file
```
vim samples.csv
vim config.yml
```
3D. Dry Run
```
snakemake -npr
```
3E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 20 --use-envmodules --rerun-incomplete --latency-wait 300 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output} --job-name {cluster.name}'"
```
