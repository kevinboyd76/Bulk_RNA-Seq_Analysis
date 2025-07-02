# Bulk_RNA-Seq_Analysis

![RNA-seq](/images/Bulk_RNA-seq.png)  
- OpenAI. (2025). Scientific data visualization: RNA-seq pipeline schematic [AI-generated image]. DALL-E. Retrieved from ChatGPT interface.
  
---

## 1) Project Description

**Bulk RNA-seq Analysis** is a reproducible Snakemake pipeline that automates the complete RNA-seq data processing workflow: **quality control**, **adapter trimming**, **alignment**, and **read quantification** with merged counts. Instead of running each step manually, the pipeline executes all tasks in a robust and modular fashion, from **raw FASTQ** inputs to **merged gene-level count matrices**, ready for downstream differential expression analysis. This pipeline takes raw FASTQs and outputs individual sample counts files and a single merged counts file. 

### Key Features

+ **Comprehensive QC**  
  + Runs **FastQC** on raw and trimmed FASTQs  
  + Aggregates results with **MultiQC**

+ **Efficient Adapter Trimming**  
  + Uses **fastp** to trim paired-end reads, generating JSON and HTML QC reports

+ **Accurate Alignment**  
  + Uses **HISAT2** to align reads to a user-specified genome index

+ **Gene-Level Quantification**  
  + Uses **featureCounts** to quantify aligned reads by gene annotation

+ **Merged Count Matrix**  
  + Merges individual count files into a single `merged_counts.csv` for direct input to DESeq2, EdgeR, or other differential expression tools

+ **Parallel Execution**  
  + Supports HPC environments via Snakemake's job management, re-running only missing or outdated outputs

---

## 2) Intended Use Case

This workflow is designed for **bulk RNA-seq** experiments where you want to:

+ Start with raw paired-end FASTQ files  
+ Perform thorough QC before and after trimming  
+ Align to a reference genome using **HISAT2**  
+ Count gene-level expression with **featureCounts**  
+ Generate a single count matrix for downstream DE analysis

---

## 3) Dependencies and Configuration

All parameters and module versions are specified in `config/config.yml`.

**Key fields include:**  
+ `Hisat2_Index`: path to your HISAT2 genome index  
+ `Annotation_GTF`: path to your reference annotation GTF file  
+ `fastqc`, `multiqc`, `fastp`, `hisat2`, `samtools`, `subread`: HPC modules or conda environments used for each step

---

## 4) Tools & Modules

The workflow uses the following tools:

+ **FastQC** — quality checks on raw and trimmed reads
+ **MultiQC** — aggregates all FastQC reports
+ **fastp** — adapter trimming with visual and JSON reports
+ **HISAT2** — alignment of paired-end reads
+ **Samtools** — sorting, indexing BAM files
+ **featureCounts** — quantification of gene-level read counts
+ **Python** — merges individual count files into a single CSV matrix

---

## 5) Example Data

Include a small test dataset in a `resources/` folder (optional). Edit `samples.csv` to test the workflow with this data before analyzing your actual samples.

---

## 6) Explanation of `samples.csv`

Your `config/samples.csv` should define each sample's paired-end FASTQ files.

| sample     | fastq1                          | fastq2                         |
|------------|---------------------------------|--------------------------------|
| **Sample_A** | /path/to/Sample_A_R1.fastq.gz | /path/to/Sample_A_R2.fastq.gz  |
| **Sample_B** | /path/to/Sample_B_R1.fastq.gz | /path/to/Sample_B_R2.fastq.gz  |
| **Sample_C** | /path/to/Sample_C_R1.fastq.gz | /path/to/Sample_C_R2.fastq.gz  |

+ **sample**: unique sample name used for output naming  
+ **fastq1/fastq2**: absolute or relative paths to the paired FASTQ files

---

## 7) Examples of Output

**1. Raw FASTQ QC**  
+ FastQC HTML reports in `results/qc/raw/fastqc/`  
+ MultiQC report in `results/qc/raw/multiqc/`

**2. Trimmed FASTQs**  
+ Adapter-trimmed FASTQs in `results/trimmed/`  
+ `fastp` HTML & JSON reports per sample

**3. QC on Trimmed FASTQs**  
+ FastQC HTML reports in `results/qc/trimmed/fastqc/`  
+ MultiQC report in `results/qc/trimmed/multiqc/`

**4. Aligned Files**  
+ Sorted BAMs in `results/aligned/`  
+ BAM index files (`.bam.bai`)

**5. Quantification Files**  
+ Individual `*.counts.txt` files per sample in `results/counts/`  
+ `merged_counts.csv` in `results/counts/` — final matrix for DESeq2 or EdgeR

---

# 8) Instructions to run on Slurm managed HPC
8A. Download version controlled repository
```
git clone https://github.com/kevinboyd76/Bulk_RNA-Seq_Analysis.git
```
8B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
8C. Modify samples and config file
```
vim samples.csv
vim config.yml
```
8D. Dry Run
```
snakemake -npr
```
8E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 20 --use-envmodules --rerun-incomplete --latency-wait 300 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output} --job-name {cluster.name}'"
```
