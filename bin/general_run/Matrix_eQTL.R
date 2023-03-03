# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.
library(parallel)
library(MatrixEQTL)
library(data.table)

# Get command line arguments
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option("--nullDist", action="store_true", default=FALSE,
              help="Should create a null distribution?  [default %default]"),
  make_option("--snpfile", action="store", default=NA, type='character',
              help="File of SNPs data based on matrix eQTL format."),
  make_option("--snplocation", action="store", default=NA, type='character',
              help="File of SNPs location data based on matrix eQTL format."),
  make_option("--expressionfile", action="store", default=NA, type='character',
              help="Gene Expression data."),
  make_option("--covariates", action="store", default=NA, type='character',
              help="File covariates based on matrix eQTL format."),
  make_option("--genelocation", action="store", default=NA, type='character',
              help="Gene location data based on matrix eQTL format."),
  make_option("--MAF", action="store", default=0.05, type='double',
              help="Minor allele frquency (MAF) filter to use [default %default]"),
  make_option("--output", action="store", default="genotype", type='character',
              help="Output prefix used to save the data. [default %default]"),
  make_option("--outdir", action="store", default="", type='character',
              help="Directory where the output data will be store. [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

###Debug
if(FALSE){
  opt <- list(nullDist = FALSE,
              snpfile = "/mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/genotyping/matrix_eQTL/snps/MAF_0.05/snp_22.txt",
              snplocation = "/mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/genotyping/matrix_eQTL/snpsloc/MAF_0.05/snppos_22.txt",
              expressionfile = "/mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/results/covariates/tumor/CD8/cell/GE.txt",
              covariates = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/covariates/LungInfo_covariates.txt",
              genelocation = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/geneAnnotation/DLCP_geneAnnot_hg38.matrixeQTL.txt",
              MAF = 0.01, #default 0.05
              output = "genotype", #default ./genotype
              outdir = "/mnt/beegfs/cgonzalez/R24/") #default NULL
}

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Building null distribution
nullDist <- opt$nullDist

# Genotype file name
SNP_file_name <- opt$snpfile
snps_location_file_name <- opt$snplocation

# Gene expression file name
expression_file_name <- opt$expressionfile
gene_location_file_name <- opt$genelocation

# Covariates file name
# Set to character() for no covariates
# covariates_file_name = character();
covariates_file_name <- opt$covariates
#####
cat("########################################################\n")
cat("###############      ARGUMENTS      ####################\n")
cat("CALCULATE NULL DISTRIBUTION: ", opt$nullDist,"\n")
cat("SNP FILE: ", opt$snpfile,"\n")
cat("SNP LOCATION FILE: ", opt$snplocation,"\n")
cat("EXPRESSION FILE: ", opt$expressionfile,"\n")
cat("GENE LOCATION FILE: ", opt$genelocation,"\n")
cat("COVARIATES FILE: ", opt$covariates,"\n")
cat("MAF: ", opt$MAF,"\n")
cat("OUTPUT PREFIX: ", opt$output,"\n")
cat("OUTPUT DIRECTORY: ", opt$outdir,"\n")
cat("########################################################\n")
####
# Output file name
dir.create(opt$outdir, recursive = T, showWarnings = F)
output_file_name_cis <- file.path(opt$outdir, paste0(opt$output, "_cis.txt"))
output_file_name_tra <- file.path(opt$outdir, paste0(opt$output, "_tra.txt"))
output_file_name_qqplot <- file.path(opt$outdir, paste0(opt$output, "_qqplot.png"))
output_file_name_all_cis <- file.path(opt$outdir, paste0(opt$output, "_all_cis.txt"))

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1e-8;
if(nullDist) {
    pvOutputThreshold_tra = 0;
}

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

# Minimum MAF
minMAF = opt$MAF;

# Maximum missingness
maxMiss = 0.05;

## Load genotype data
if(nrow(fread(SNP_file_name, header = TRUE, nrows = 10)) == 0) {
	# Touch the output file and exit if there are no suitable snps
	try(system(paste("touch", output_file_name_cis),
		intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_tra),
		intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_qqplot),
		intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_all_cis),
		intern = TRUE, ignore.stderr = TRUE))
	quit(save = "no")
}


snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# Subset and reorder snps columns using the gene columns as reference
cols <- sapply(gene$columnNames, function(x) which(x == snps$columnNames))
snps$ColumnSubsample(cols)

## ## Load annotation, adding sex covariate
## Load covariates
covariates_file_name <- read.table(covariates_file_name, header = T, check.names = F)

if(nrow(covariates_file_name)>0) {
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(covariates_file_name));

  colscvrt <- sapply(gene$columnNames, function(x) which(x == cvrt$columnNames))
  cvrt$ColumnSubsample(colscvrt)

}else{
  cvrt <- SlicedData$new()
}
## Permutate samples if we are calculating the null distribution
if(nullDist) {
	random <- sample(1:snps$nCols(), size = snps$nCols());
	snps$ColumnSubsample(random);
}

## Filter SNPs with high missingness
miss.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
	slice = snps[[sl]];
	miss.list[[sl]] = unlist(apply(slice, 1,
		function(x) sum(is.na(x))))/ncol(slice);
}
miss = unlist(miss.list)
cat('SNPs before missingness filtering:',nrow(snps),'\n')
snps$RowReorder(miss < maxMiss);
cat('SNPs after missingness filtering:',nrow(snps),'\n')

## Filter SNPs with low MAF

maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
	slice = snps[[sl]];
	maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
	maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

## Look at the distribution of MAF
# hist(maf[maf<0.1],seq(0,0.1,length.out=100))
cat('SNPs before MAF filtering:',nrow(snps),'\n')
if(any(maf > minMAF)) {
	snps$RowReorder(maf>minMAF);
	cat('SNPs after MAF filtering:',nrow(snps),'\n')
} else {
	# Touch the output file and exit if there are no suitable snps
	try(system(paste("touch", output_file_name_cis),
		intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_tra),
		intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_qqplot),
		intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_all_cis),
		intern = TRUE, ignore.stderr = TRUE))
	quit(save = "no")
}

# Quantile normalization of the gene expression values
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);

## Run the analysis

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);


me = Matrix_eQTL_main(
		snps = snps,
		gene = gene,
		cvrt = cvrt,
		output_file_name = NULL,
		pvOutputThreshold = pvOutputThreshold_tra,
		useModel = useModel,
		errorCovariance = errorCovariance,
		verbose = TRUE,
		output_file_name.cis = NULL,
		pvOutputThreshold.cis = pvOutputThreshold_cis,
		snpspos = snpspos,
		genepos = genepos,
		cisDist = cisDist,
		pvalue.hist = "qqplot",
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
# cat('Detected local eQTLs:', '\n');
# #show(me$cis$eqtls)
# cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

if(!nullDist) {
	cat('Saving all cis-eQTLs\n')
	#write.table(me$cis$eqtls, output_file_name_all_cis,
	#	row.names = FALSE, quote = FALSE, sep="\t")
  data.table::fwrite(me$cis$eqtls, output_file_name_all_cis,
	 	row.names = FALSE, quote = FALSE, sep="\t")
  saveRDS(me$cis$eqtls, file = gsub(".txt", ".rds", output_file_name_all_cis))
}

cat('Saving best SNP per gene table\n')
cisDT <- data.table(me$cis$eqtls, key="gene")
cisDT$chr <- gsub(":.+$", "", cisDT$snps)
bestCis <- cisDT[,.SD[which.min(pvalue)], by=gene] #more memmory needed
cat("directory exists?")
print(dir.exists(opt$outdir))
write.table(bestCis, output_file_name_cis,
	row.names = FALSE, quote = FALSE, sep="\t")

if(!nullDist) {
  cat('Saving all trans-eQTLs\n')
	transDT <- data.table(me$trans$eqtls, key="gene")
	#
  transDT$chr <- gsub(":.+$", "", transDT$snps)
  bestTrans <- transDT[,.SD[which.min(pvalue)], by=gene]
  # write.table(bestTrans, output_file_name_tra,
	# 	row.names = FALSE, quote = FALSE, sep="\t")
  data.table::fwrite(bestTrans, output_file_name_tra,
	 	row.names = FALSE, quote = FALSE, sep="\t")
  saveRDS(bestTrans, file = gsub(".txt", ".rds", output_file_name_tra))
}

## Plot the qq-plot of all p-values
if(nrow(me$cis$eqtls) > 0 ){
  png(filename = output_file_name_qqplot, width = 650, height = 650)
  plot(me, pch = 16, cex = 0.7)
  dev.off()
}
