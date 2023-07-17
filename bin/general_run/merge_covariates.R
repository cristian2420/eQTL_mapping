#!/usr/bin/R

# ------------------------------------------------------------------------------
# title:  Merged Covariates
# author: Cristian Gonzalez-Colin
# email: cgonzalez@lji.org
# date: mar 11, 2022
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-p", "--pcafile"), action="store", default=NA, type='character',
              help="File with PCA data. From plink output."),
  make_option(c("-f", "--factorfile"), action="store", default=NA, type='character',
              help="Factor file, usually peer factors with donors in rows."),
  make_option(c("-c", "--covariatesFile"), action="store", default=NA, type='character',
              help="File with known covariates with donors in columns."),
  make_option(c("-v", "--covariates"), action="store", default=NA, type='character',
              help="Covariates to use in analysis (included in the covariates file) separate by '-'."),
  make_option(c("-n", "--factors"), action="store", default=20, type = "integer",
              help="Number of factors to merge [default %default]"),
  make_option(c("-s", "--pcs"), action="store", default=6, type = "integer",
              help="Number of PCs to merge [default %default]"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="File to save merge data."),
  make_option(c("-q", "--fastqtl"), action="store_true", default = FALSE, type = 'logical',
              help="Makes the output file compatible with fastqtl pipeline [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

##debug
if(FALSE){
  opt = list(
    pcafile = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/results/matrix_eqtl/covariates/tumor/CD4/cell/genotype.eigenvec",
    factorfile = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/results/matrix_eqtl/covariates/tumor/CD4/cell/peer_factors.csv",
    covariatesFile = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/covariates/LungInfo_covariates.txt",
    factors = 0,
    pcs = 0,
    covariates = 'Gender',
    output = "",
    fastqtl = FALSE
  )
}
##
covariates <- read.table(opt$covariatesFile, header = T, row.names = 1, check.names = F)
if(is.na(opt$covariates)){
  needcov <- FALSE
}else{
  covuse <- strsplit(opt$covariates, '-')[[1]]
  covariates <- covariates[covuse,]
  needcov <- TRUE
}

##
pcaFile <- read.table(opt$pcafile, header = T, check.names = F)
pcaFile$FID <- NULL
pcaname <- pcaFile[,1]
pca <- data.frame(t(pcaFile[,-1]))
colnames(pca) <- pcaname
if(nrow(pca) < opt$pcs) stop("Number of PCs to use greater than PCs in the pca file\n")
pca <- data.frame(pca[1:opt$pcs, ], check.names = F)
############ Try to use the specified peer factos but if not have enough it used the maximum calculated that are less than the specified
if(!is.na(opt$factorfile)){
  peerFile <- read.csv(opt$factorfile, header = T)
  peername <- peerFile[,1]
  peer <- data.frame(t(peerFile[,-1]))
  colnames(peer) <- peername
  #if(nrow(peer) < opt$factors) stop("Number of peer factors to use greater than factors in the peer file\n")
  if(nrow(peer) < opt$factors){
    cat("WARNING: Number of peer factors to use greater than factors in the peer file\n")
    opt$factors <- nrow(peer)
  }
  peer <- peer[1:opt$factors,]
}
##
covariates <- covariates[, colnames(pca)]
#####
if( opt$factors > 0 &  opt$pcs > 0 & needcov){
  dfAll <- Reduce(rbind, list(pca, peer, covariates))
}else if( opt$factors > 0 &  opt$pcs == 0 & needcov){
  dfAll <- Reduce(rbind, list(peer, covariates))
}else if( opt$factors == 0 &  opt$pcs > 0 & needcov){
  dfAll <- Reduce(rbind, list(pca, covariates))
}else if( opt$factors > 0 &  opt$pcs > 0 & !needcov){
  dfAll <- Reduce(rbind, list(pca, peer))
}else if(opt$factors == 0 &  opt$pcs == 0 & needcov){
  dfAll <- covariates
}else if(opt$factors > 0 &  opt$pcs == 0 & !needcov){
  dfAll <- peer
}else if(opt$factors == 0 &  opt$pcs > 0 & !needcov){
  dfAll <- pca
}else{
  dfAll <- covariates[0,]
}
####
###save with fastqtl format
if(opt$fastqtl){
  dfAll$id <- rownames(dfAll)
  dfAll <-cbind(dfAll[,colnames(dfAll) == 'id', drop =F], dfAll[,colnames(dfAll) != 'id'])
  write.table(dfAll, opt$output, quote = F, sep = "\t", row.names = F)
  system(paste0('bgzip ', opt$output))
}else{
  write.table(dfAll, opt$output, quote = F, sep = "\t")
}
