#!/usr/bin/R

# ------------------------------------------------------------------------------
# title:  Peer factors by MOFA2 package
# author: Cristian Gonzalez-Colin
# script modify from https://github.com/bioFAM/PEER
# email: cgonzalez@lji.org
# date:
# Goal: ## Script to train a MOFA model from R - peer factors ##
# if train_model parameter not present, it only subset the expression matrix by donor_list
# ------------------------------------------------------------------------------
#Adapt from peer wrapper with MOFA model: https://github.com/bioFAM/PEER
#based on closed issue (https://github.com/mz2/peer/pull/15) on original peer github

######
Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/bin/python") # line added due to a problem with regular installation please comment or change the path if need it
library(reticulate)
library(edgeR)
# Important: a part from the MOFA2 R package, you also need the python package mofa2 to be installed
library(MOFA2)
library(data.table)

#######
 suppressPackageStartupMessages(require(optparse))

 option_list = list(
   make_option(c("-e", "--expFile"), action="store", default=NA, type='character',
               help="Gene expression file"),
   make_option(c("-d", "--donorFile"), action="store", default=NA, type='character',
               help="List of donors to keep"),
   make_option(c("-t", "--train_model"), action="store_true", default=FALSE,
               help="Should the program get the peer factors? [default %default]"),
   make_option(c("-p", "--peerFactors"), action="store", default=NULL, type='character',
               help="File to store the calculate factors"),
   make_option(c("-n", "--numberPEERFactors"), action="store", default=20, type='numeric',
               help="File to store the calculate factors. [default %default]"),
   make_option(c("-o", "--outputFile"), action="store", default=NULL, type='character',
               help="File to store the gene expression matrix"),
   make_option(c("-q", "--quantileNorm"), action="store_true", default = FALSE, type = 'logical',
               help="Should data need quantile normalization? [default %default]"),
   make_option(c("-f", "--fastqtl"), action="store_true", default = FALSE, type = 'logical',
               help="Makes the output file compatible with fastqtl pipeline [default %default]"),
   make_option(c("-a", "--annotFile"), action="store", default=NA, type='character',
               help="If '--fastqtl' specified, an annotation file should be provided. File must contain at least 3 columns with the following order [geneid, chr, start, end]")
 )
 opt = parse_args(OptionParser(option_list=option_list))
####
#Debug
if( FALSE ){
  opt = list(
    expFile = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/expression/input/tumor/CD4/cell/gene_expression.txt",
    donorFile = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/expression/input/tumor/CD4/cell/donors.txt",
    quantileNorm = TRUE,
    train_model = TRUE,
    peerFactors = "/home/cgonzalez/myscratch/R24/DICE_cancerProject/peerFactors.csv",
    numberPEERFactors = 20,
    outputFile = "/home/cgonzalez/myscratch/R24/DICE_cancerProject/expMatrix.txt",
    fastqtl = FALSE,
    annotFile = '/mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/geneAnnotation/DLCP_geneAnnot_hg38.matrixeQTL.txt')
}
########
cat('***********************\n************Parameters:\n')
cat("Expression File:", opt$expFile, '\n')
cat("Donor File:", opt$donorFile, '\n')
cat("Quantile normalization:", opt$quantileNorm, '\n')
cat("Running peer factor:", opt$train_model, '\n')
cat("Number of peer factors:", opt$numberPEERFactors, '\n')
cat("Peer Factors File:", opt$peerFactors, '\n')
cat("Output File:", opt$outputFile, '\n')
if(opt$fastqtl){
  cat("Output with FastQTL format:", opt$fastqtl, '\n')
  cat("annotFile File:", opt$annotFile, '\n')
}
cat('***********************\n')
#####
donorFile <- readLines(opt$donorFile)
expFile <- read.table(opt$expFile, header = T, check.names = F, row.names = 1)
####filter donors
expFile <- expFile[, c(donorFile)]
####normalization
y <- DGEList(expFile)
#/ calculate TMM normalization factors:
y <- calcNormFactors(y)
#/ get the normalized counts if needed
cpms <- cpm(y, log=FALSE)
if(opt$quantileNorm){
  
  mat = t(apply(cpms, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(cpms)+1));
  cpms = mat
}

###save with fastqtl format
if(opt$fastqtl){
  header_names <- c("TargetID", "#Chr", "start", "end")
  header_order <- c("#Chr", "start", "end", "TargetID")
  annot <- read.table(opt$annotFile, header = T, check.names = F)
  names(annot) <- header_names
  comp_names <- colnames(mat)
  mat <- as.data.frame(mat)
  mat$TargetID <- rownames(mat)
  mat <- merge(annot, mat, by = 'TargetID')
  mat <- mat[, c(header_order, comp_names)]
  ##ordering and removing non-canonical chromosome for compatibilities with genotyping
  mat <- mat[order(mat[,1], mat[,2]),]
###    mat <- mat[mat[,1] %in% c(1:22, 'X', 'Y'),]
  mat <- mat[mat[,1] %in% c(1:22),]

  write.table(mat, opt$outputFile, quote = F, sep = "\t", row.names = F)
  system(paste0('bgzip ', opt$outputFile))
  system(paste0('tabix -f -p bed ', opt$outputFile, '.gz'))
}else{
  write.table(cpms, opt$outputFile, quote = F, sep = "\t")
}


if(opt$train_model){
  flg <- 1
  if(opt$numberPEERFactors < 1 | opt$numberPEERFactors > length(donorFile)){
    cat(paste0('WARNING: Number of Peer Factors should be greater than 1 and/or lesser than number of Donors (',  length(donorFile), ')\n'))
    opt$numberPEERFactors <- 10
    flg <- 0
  }
  if(opt$numberPEERFactors == 1) flg <- 0
  # MOFA is a multi-view factor analysis framework that is a generalisation of PEER
  # The data needs to be input as a list of views. If you have a single data modality, the input data
  # corresponds to a list with a single element.
  data <- list(as.matrix(cpms))

  ########################
  ## Create MOFA object ##
  ########################

  object <- suppressMessages(create_mofa(data))

  #########################
  ## Define data options ##
  #########################

  data_opts <- get_default_data_options(object)

  ##########################
  ## Define model options ##
  ##########################

  model_opts <- get_default_model_options(object)

  # use ARD prior for the factors? (please do not edit this)
  model_opts$ard_factors <- FALSE

  # number of factors
  #model_opts$num_factors <- 10
  model_opts$num_factors <-  opt$numberPEERFactors

  #############################
  ## Define training options ##
  #############################

  train_opts <- get_default_training_options(object)

  # maximum number of iterations
  train_opts$maxiter <- 2000

  # fast, medium, slow
  train_opts$convergence_mode <- "fast"

  # initial iteration to start evaluating convergence using the ELBO (recommended >1)
  train_opts$startELBO <- 1

  # frequency of evaluation of ELBO (recommended >1)
  train_opts$freqELBO <- 1

  # use GPU (needs cupy installed and working)
  train_opts$gpu_mode <- FALSE

  # verbose output?
  train_opts$verbose <- FALSE

  # random seed
  train_opts$seed <- 1

  #############################
  ## Prepare the MOFA object ##
  #############################

  object <- prepare_mofa(
    object = object,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )

  ##############
  ## Run MOFA ##
  ##############
  model_name <- file.path(dirname(opt$peerFactors), gsub('\\..+', '.hdf5', basename(opt$peerFactors)))
  model_plotname <- file.path(dirname(opt$peerFactors), gsub('\\..+', '.pdf', basename(opt$peerFactors)))
  model_var <- file.path(dirname(opt$peerFactors), gsub('\\..+', '_variance.csv', basename(opt$peerFactors)))
  model <- run_mofa(object, outfile=model_name)

  factors <- get_factors(model, factors = "all")
  #weights <- get_weights(model, views = "all", factors = "all")
  ###

  # covariates <- read.table('data/covariates/LungInfo_covariates.txt', sep = '\t', heade = T, row.names = 1)
  # covariates <- covariates[,donorFile]
  # covariates$id <- rownames(covariates)
  # metadata <- data.frame(data.table::dcast(data.table::melt(covariates, id.vars = "id"), variable ~ id) )
  # metadata$group <- 'group1'
  # data.table::setnames(metadata, 'variable', 'sample')
  # metadata$sample <- as.character(metadata$sample)
  # samples_metadata(model) <- metadata
  explvar <- calculate_variance_explained(model, factors = 'all')
  if(flg){
    pdf(model_plotname)
    #############################
    ## Plot variance explained ##
    #############################

    # Plot variance explained using individual factors
    plot_variance_explained(model, factors="all")
    plot_variance_explained(model, factors=c(1,2,3))
    ##elbow plot
    explvar <- calculate_variance_explained(model, factors = 'all')
    plot(1:length(explvar$r2_per_factor$group1), explvar$r2_per_factor$group1, type = "b", pch = 19,
       col = "red", xlab = "x", ylab = "y")
    title(main="Scree plot", xlab="Peer Factors", ylab="Percentage of variance explained %")
    # Plot total variance explained using all factors
    plot_variance_explained(model, plot_total = TRUE)[[2]]

    ########################
    ## Plot factor values ##
    ########################

    plot_factor(model, factor = 1)
    plot_factors(model, factor = c(1,2))

    ###########################
    ## Plot feature loadings ##
    ###########################

    # The weights or loadings provide a score for each gene on each factor.
    # Genes with no association with the factor are expected to have values close to zero
    # Genes with strong association with the factor are expected to have large absolute values.
    # The sign of the loading indicates the direction of the effect: a positive loading indicates that the feature is more active in the cells with positive factor values, and viceversa.

    # Plot the distribution of loadings for Factor 1.
    plot_weights(model,
      view = 1,
      factor = 1,
      nfeatures = 10,     # Top number of features to highlight
      scale = T           # Scale loadings from -1 to 1
    )

    # If you are not interested in the full distribution, but just on the top loadings:
    plot_top_weights(model,
      view = 1,
      factor = 1,
      nfeatures = 10,
      scale = T
    )

    ######################################
    ## Plot correlation between factors ##
    ######################################

    plot_factor_cor(model)
    dev.off()
    ###################
    ## Fetch factors ##
    ###################
  }
  write.csv(explvar$r2_per_factor$group1, model_var)
  write.csv(factors$group1, opt$peerFactors )
  ########
}
##
