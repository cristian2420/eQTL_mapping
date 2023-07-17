# eQTL Mapping Pipeline

------
* Cristian Gonzalez-Colin (cgonzalez@lji.org)
* Vijayanand Lab (https://www.lji.org/labs/vijayanand/)
* La Jolla Institute for Immunology (LJI)
* La Jolla, CA USA
* Current version: (03/02/2023)
------

## About it

This pipeline was developd for eQTL calling for the DICE Tissue project (unpublished). It was implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/) v7.14.0 workflow manager. Cluster configuration file (cluster.json) needs to be modified according to the cluster/cloud enviroment to work properly.

Linear models are fitted using [MatrixeQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) and two different multiple correction methods are included: [eigenMT](https://github.com/joed3/eigenMT) and a permutation-based method. The last one requires more computational resources and time to evaluate 1,000 (default) permutations.

## Pipeline setup

In order to properly run the pipeline three environment are provided. Snakemake would automatically set the environments to run specific steps.

* DLCP.yaml
* pyEigenMT.yaml
* bcftools.yaml


## Data preparation

### Config file:

Configuration file ```snake_conf.yaml``` has to be in the same folder as the ```Snakefile```. Make proper changes to it. Files needed are explained below.

### Covariates file:

Tab separated file with donors in the columns and covariates in the rows. This is a general covariates file, PEER Factors and PCs will be added within the pipeline.

### Snpfiles file:

Tab separated file with the location of genotype files in matrix eQTL format:
|CHR|SNP|SNP_LOC|
|---|---|-------|
|1|SNPFILE.txt|SNPLOC.txt|
|2|SNPFILE.txt|SNPLOC.txt|

Both files need to be prepared followed matrix eQTL [toy dataset](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt). [SNPFILE.txt](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt) & [SNPLOC.txt](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/snpsloc.txt)

### DATASETS file:

Pipeline runs multiple clusters per cell in parallel. To do that a data set .csv file is needed with the following columns.  
|cell|tissue|subset|expFile|donorFile|
|----|------|------|-------|---------|
|CD4|Lung|0|file1.txt|donors1.txt|
|CD4|Lung|1|file2.txt|donors2.txt|
|CD8|Lung|0|file3.txt|donors3.txt|

* expFile: Tab separated file with expression data to use for a given cluster in a cell type. Donors in columns and genes in rows. Name of donors and genes are needed.
* donorFile: List of donors to use in the analysis. Needs to be a set of all donors.

### Gene annotation file [geneannot]:

File was prepared following matrix eQTL [toy dataset](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/geneloc.txt).

### NOTES

The ```bin/``` folder needs to be copy to the working directory.

The ```pbs_submit.sh``` file shows and example of how the pipeline was run.

Please cite the following manuscript if you are using this repository:

## Contact

Please email Cristian Gonzalez-Colin (cgonzalez@lji.org)  and/or Vijayanand Pandurangan (vijay@lji.org).
