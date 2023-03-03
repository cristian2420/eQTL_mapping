# UNDER CONSTRUCTION: eQTL Mapping Pipeline

------
* Cristian Gonzalez-Colin (cgonzalez@lji.org)
* Vijayanand Lab (https://www.lji.org/labs/vijayanand/)
* La Jolla Institute for Immunology (LJI)
* La Jolla, CA USA
* Current version: (03/02/2023)
------

## About it

This pipeline was developd for eQTL calling for the DICE Tissue project (unpublished). It was implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/) v7.14.0 workflow manager. Cluster configuration file (cluster.json) needs to be modified according to the cluster/cloud enviroment to work properly.

Linear models are fitted using MatrixeQTL and two different multiple correction methods are included: eigenMT and a permutation-based method. The last one requires more computational resources and time to evaluate 1,000 (default) permutations.

## Pipeline setup

The following tools/versions have been used in the pipeline. Different versions of these tools may results in different results. Version used is specified in parentheses.

*
*
*
*



## Data preparation

### Config file:

### Running file:



## Contact

Please email Cristian Gonzalez-Colin (cgonzalez@lji.org)  and/or Vijayanand Pandurangan (vijay@lji.org).
