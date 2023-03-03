import pandas as pd
import os
import re

configfile: "snake_conf.yaml"
#configfile: "snake_conf_10x.yaml"
workdir: config['config']['workdir']
#report: "report/workflow.rst"

class Table:
    def __init__(self, file):
        self.data = pd.read_csv(file)
    def extract(self, wildcards):
        #sel = datasets[ datasets["cell"] == wildcards.cell ]
        sel = self.data[ self.data["cell"] == wildcards.cell ]
        sel = sel[ sel["subset"].astype(str) == wildcards.subset ]
        sel = sel[ sel["tissue"] == wildcards.tissue ]
        return sel.to_dict(orient="list")

    def donor(self, wildcards):
        dic = self.extract(wildcards)
        return {"donorFile": dic["donorFile"]}

    def expression(self, wildcards):
        dic = self.extract(wildcards)
        return {"expFile": dic["expFile"]}

def snpfiles(wildcards):
    dic = SNP_FILES[ SNP_FILES['CHR'] == wildcards.chrom ].to_dict(orient = 'list')
    return {"snps": dic["SNP"], "snpsloc": dic["SNP_LOC"]}


#### Load dataset file
dataset_file = config["inpFiles"]["datasets"]
datasets = Table(dataset_file)

TISSUES = datasets.data.to_dict(orient='list')['tissue']
CELLS = datasets.data.to_dict(orient='list')['cell']
SUBSETS = datasets.data.to_dict(orient='list')['subset']
PERMUTATIONS = int(config["parameters"]["permutations"])
#######
SNP_FILES = pd.read_table(config["inpFiles"]["snpfiles"])
CHROMS =  SNP_FILES['CHR'].to_list()
###
covariates = config['parameters']['covariates']
covariates = '-'.join(covariates)

run_var_genes = config['parameters']['run_var_genes']
n_var_genes = config['parameters']['n_var_genes']
print(run_var_genes)
print(n_var_genes)
mst_var_gns = '_VarGenes' + str(run_var_genes) + '_Ngenes_' + str(n_var_genes)
print(mst_var_gns)

if run_var_genes:
    if type(n_var_genes) != int or n_var_genes < 1:
        sys.exit("ERROR: if 'run_var_genes' == True, 'n_var_genes' int > 0 & < number genes in the matrix")


GPREFIX = "MAF_" + str(config['parameters']['MAF']) + "_covariates" + covariates + "_genPCs" + str(config['parameters']['PCS']) + "_expPEER" + str(config['parameters']['Factors']) + mst_var_gns
######
rule all:
    input:
        expand("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig.tsv", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS)
        #expand("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig.tsv", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS)

    #expand(expand("results/matrix_eqtl/eQTL/{{chrom}}/{tissue}/{cell}/{subset}/Output_All_cis_sig.tsv", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS ), chrom=CHROMS)


rule expression_preprocess:
    input:
        unpack(datasets.donor),
        unpack(datasets.expression)
    output:
        exp_file = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.txt")
    params:
        script = "bin/general_run/process_expression.R",
        peerfactors_n = config["parameters"]["Factors"]
    shell:
        "Rscript {params.script} -e {input.expFile} -d {input.donorFile} -q -o {output.exp_file} "

rule peer_calc:
    input:
        expression = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.txt"
    output:
        "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/factors.txt"
    params:
        script = "bin/general_run/run_peer_analysis.py",
        peerfactors_n = config["parameters"]["Factors"],
        odir = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/",
        run_var_genes = run_var_genes,
        n_genes = n_var_genes
    shell:
        "~/miniconda3/envs/peer/bin/python {params.script} {input.expression} {params.peerfactors_n} {params.odir} {params.run_var_genes} {params.n_genes} "

rule subset_geno:
    input:
        unpack(datasets.donor),
        genotype = config["inpFiles"]["vcf_genotype"]
    output:
        temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.vcf.gz")
    shell:
        "bcftools view -S {input.donorFile} -Oz -o {output}  {input.genotype}"

rule LD_prune:
    input:
        vcf = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.vcf.gz"
    output:
        prunein = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.prune.in",
        pruneout = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.prune.out")
    params:
        "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype"
    shell:
        "plink --vcf {input.vcf} --indep-pairwise 200 100 0.1 --out {params}"

rule pca:
    input:
        genotype = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.vcf.gz",
        variants = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.prune.in"
    output:
        eigenvec = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenvec"),
        eigenval = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenval")
    params:
        "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype"
    shell:
        "plink --vcf {input.genotype} --extract {input.variants} --pca header tabs --out {params}"

rule merge_cov:
    input:
        covariates = config["inpFiles"]["covariates"],
        peer_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/factors.txt",
        eigenvec = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenvec"
    output:
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt"
    params:
        script = "bin/general_run/merge_covariates.R",
        pcs = config["parameters"]["PCS"],
        factors = config["parameters"]["Factors"],
        covar = covariates
    run:
        if not params.covar:
            shell("Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --factors {params.factors} --pcs {params.pcs} --output {output.covariates} ")
        else:
            shell("Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --covariates {params.covar} --factors {params.factors} --pcs {params.pcs} --output {output.covariates} ")
    # shell:
    #     "Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --covariates {params.covar} --factors {params.factors} --pcs {params.pcs} --output {output.covariates} "


rule matrix_eQTL:
    input:
        # snps = config["inpFiles"]["snps"],
        # snpsloc = config["inpFiles"]["snpsloc"],
        unpack(snpfiles),
        expression = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.txt",
        geneloc = config["inpFiles"]["geneannot"],
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt"
    output:
        cis = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_cis.txt",
        trans = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_tra.txt",
        cis_all = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_all_cis.txt",
        qqplot = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_qqplot.png"
    params:
        script = "bin/general_run/Matrix_eQTL.R",
        null = "FALSE",
        MAF = config["parameters"]["MAF"],
        outputdir = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/",
        prefix = "Output"
    shell:
        "Rscript {params.script}  --nullDist {params.null} --snpfile {input.snps} \
         --snplocation {input.snpsloc} --expressionfile {input.expression} \
         --covariates {input.covariates} --genelocation {input.geneloc} --MAF {params.MAF} \
         --output {params.prefix} --outdir {params.outputdir} "

rule eigenMT:
    input:
        unpack(datasets.donor),
        unpack(snpfiles),
        qtl = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_all_cis.txt",
        geneloc = config["inpFiles"]["geneannot"]
    output:
        eigen = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_eigen_cis_sig.tsv"
    params:
        script = "bin/general_run/eigenMT.py",
        chromosome = "{chrom}"
    shell:
        "python {params.script} --QTL {input.qtl} --GEN {input.snps} --GENPOS {input.snpsloc} --PHEPOS {input.geneloc} --OUT {output.eigen} --sample_list {input.donorFile} --CHROM {params.chromosome} "

rule null_distribution:
    input:
        # snps = config["inpFiles"]["snps"],
        # snpsloc = config["inpFiles"]["snpsloc"],
        unpack(snpfiles),
        expression = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.txt",
        geneloc = config["inpFiles"]["geneannot"],
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt"
    output:
        cis = temp("temp/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_{permutation}_cis.txt"),
        qqplot = temp("temp/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_{permutation}_qqplot.png")
    params:
        script = "bin/general_run/Matrix_eQTL.R",
        null = "TRUE",
        MAF = config["parameters"]["MAF"],
        outputdir = "temp/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/",
        prefix = "Output_{permutation}"
    shell:
        "Rscript {params.script}  --nullDist {params.null} --snpfile {input.snps} "
        " --snplocation {input.snpsloc} --expressionfile {input.expression} "
        " --covariates {input.covariates} --genelocation {input.geneloc} --MAF {params.MAF} "
        " --output {params.prefix} --outdir {params.outputdir} "


rule merge_null:
    input:
        null = expand("temp/matrix_eqtl/null_distr/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_{permutation}_cis.txt", chrom = CHROMS, permutation = list(range(1,PERMUTATIONS + 1)))
    output:
        null = temp("results/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/null_distr.tsv")
    params:
        Rfile = "bin/general_run/Merge_Results_null.R"
    shell:
        "Rscript {params.Rfile} {input.null} {output.null}"


rule calculate_FDR:
    input:
        cis = expand("results/matrix_eqtl/eQTL/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_cis.txt", chrom = CHROMS),
        tra = expand("results/matrix_eqtl/eQTL/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_tra.txt", chrom = CHROMS),
        null = "results/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/null_distr.tsv",
        eigen = expand("results/matrix_eqtl/eQTL/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_eigen_cis_sig.tsv", chrom = CHROMS)
    output:
        cis = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig.tsv",
        tra = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_tra_sig.tsv",
        qqcis = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig_QQplot.png"
    params:
        Rfile = "bin/general_run/Calculate_FDR_merge.R",
        prefix = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output",
        cis_results = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/",
        cis_name = "Output_cis.txt",
        tra_results = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/",
        tra_name = "Output_tra.txt",
        eigen_results = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/",
        eigen_name = "Output_eigen_cis_sig.tsv"
    shell:
        "Rscript {params.Rfile} -c {params.cis_results} -t {params.tra_results} -i {params.cis_name} -r {params.tra_name} -n {input.null} -o {params.prefix} -m {params.eigen_results} -e {params.eigen_name} "
