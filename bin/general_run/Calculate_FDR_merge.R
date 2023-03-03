###############################################################################
# Libraries ------------------------------------------------------------------#
###############################################################################
library(data.table)
library(ggplot2)
library(qvalue)
library(ggthemes)
library(RColorBrewer)


###############################################################################
# Functions ------------------------------------------------------------------#
###############################################################################

fte_theme <- function() {

	# Generate the colors for the chart procedurally with RColorBrewer
	palette <- brewer.pal("Greys", n=9)
	color.background = "White" #palette[2]
	color.grid.major = palette[5] #palette[3]
	color.axis.text = palette[9] #palette[6]
	color.axis.title = palette[9] #palette[7]
	color.title = palette[9]

	# Begin construction of chart
	theme_bw(base_size=18) +

	# Set the entire chart region to a light gray color
	theme(panel.background=element_rect(fill=color.background, color=color.background)) +
	theme(plot.background=element_rect(fill=color.background, color=color.background)) +
	theme(panel.border=element_rect(color=color.background)) +

	# Format the grid
	theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
	theme(panel.grid.minor=element_blank()) +
	theme(axis.ticks=element_blank()) +

	# Format the legend
	theme(legend.position="right") +
	theme(legend.background = element_rect(fill=color.background)) +
  theme(legend.title = element_text(size=18,color=color.axis.title)) +
	theme(legend.text = element_text(size=14,color=color.axis.text)) +

	# Set title and axis labels, and format these and tick marks
	theme(plot.title=element_text(color=color.title, size=20, vjust=1.25)) +
	theme(axis.text.x=element_text(size=16,color=color.axis.text)) +
	theme(axis.text.y=element_text(size=16,color=color.axis.text)) +
	theme(axis.title.x=element_text(size=18,color=color.axis.title, vjust=0)) +
	theme(axis.title.y=element_text(size=18,color=color.axis.title, vjust=1.25)) +

	# Plot margins
	theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}


ggd.qqplot <- function(pvector, main=NULL, ...) {
  o <- -log10(sort(pvector,decreasing=F))
  e <- -log10(1:length(o)/length(o))
  df <- data.frame(o, e)
  ggplot(data=df, aes(x=e, y=o)) +
    geom_point(size=5) +
    geom_abline(intercept = 0, slope = 1, size = 2, color="red") +
    xlab(expression(Expected~~-log[10](italic(p)))) +
    ylab(expression(Observed~~-log[10](italic(p)))) +
    xlim(c(0,max(e))) + ylim(c(0,max(o))) +
    ggtitle(main) +
		fte_theme()
}


###############################################################################
# Code -----------------------------------------------------------------------#
###############################################################################
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-c", "--cispath"), action="store", default=NA, type='character',
              help="PATH with cis results"),
  make_option(c("-t", "--trapath"), action="store", default=NA, type='character',
              help="PATH with trans results"),
	make_option(c("-i", "--cisname"), action="store", default=NA, type='character',
              help="File with cis results"),
  make_option(c("-r", "--traname"), action="store", default=NA, type='character',
              help="File with trans results"),
  make_option(c("-e", "--eigenMT"), action="store", default=NA, type='character',
              help="File with eigenMT results"),
	make_option(c("-m", "--eigenMTpath"), action="store", default=NA, type='character',
              help="PATH with eigenMT results"),
	make_option(c("-n", "--null"), action="store", default=NA, type='character',
              help="File null distribution data"),
  make_option(c("-o", "--output"), action="store", default="Output", type='character',
              help="File to save output"),
  make_option(c("-u", "--useChrNotation"), action="store_true", default=FALSE,
              help="Add Chr to chrs paths")
)
opt = parse_args(OptionParser(option_list=option_list))
##debug
if(FALSE){
	opt = list(
		cispath = "results/matrix_eqtl/eQTL/healthy/CD4/0/MAF_0.05_covariates_genPCs0_expPEER0_VarGenesFalse_Ngenes_NA/",
		trapath = "results/matrix_eqtl/eQTL/healthy/CD4/0/MAF_0.05_covariates_genPCs0_expPEER0_VarGenesFalse_Ngenes_NA/",
		cisname = "Output_cis.txt",
		traname = "Output_tra.txt",
		null = "results/matrix_eqtl/null_distr/healthy/CD4/0/MAF_0.05_covariates_genPCs0_expPEER0_VarGenesFalse_Ngenes_NA/null_distr.tsv",
		eigenMT = "Output_eigen_cis_sig.tsv",
		eigenMTpath = "results/matrix_eqtl/eQTL/healthy/CD4/0/MAF_0.05_covariates_genPCs0_expPEER0_VarGenesFalse_Ngenes_NA/",
		output = "Output"
	)
}

if (opt$useChrNotation){
	chrs_vec <- paste0("chr", c(1:22, 'X')	)
}else{
	chrs_vec <- c(1:22, 'X')
}
print(paste0("Reading cis files"))
cis <- data.table::rbindlist(lapply(chrs_vec, function(chr){
		cis <- data.table::fread(paste(opt$cispath, chr, opt$cisname, sep = "/"), sep = "\t", header=TRUE)
}))
###
print(paste0("Reading trans files"))
tra <- data.table::rbindlist(lapply(chrs_vec, function(chr){
		tra <- data.table::fread(paste(opt$trapath, chr, opt$traname, sep = "/"), sep = "\t", header=TRUE)
}))
###
print(paste0("Reading null distribution file"))
null <- read.table(opt$null, sep = "\t", header=F, row.names=1)
print(paste0("Reading eigenMT file"))
eigen <- data.table::rbindlist(lapply(chrs_vec, function(chr){
		eigen <- data.table::fread(paste(opt$eigenMTpath, chr, opt$eigenMT, sep = "/"), sep = "\t", header=TRUE)
}))
eigen <- eigen[, c("gene", "BF", "TESTS")]
eigen$BF.FDR <- p.adjust(eigen$BF, method = "fdr")
print("Calculating cis FDR")
#pvaluesCis <- empPvals(stat = abs(as.numeric(cis$statistic)),
#                       stat0 = abs(as.matrix(null)), pool = TRUE)
pvaluesCis <- empPvals(stat = as.numeric(cis$statistic),
                      stat0 = as.matrix(null), pool = TRUE)
qobj <- qvalue(p = pvaluesCis)
#p.adj <- p.adjust(pvaluesCis, method = "BH")
cis$p.value.pool <- pvaluesCis
cis$FDR.pool <- qobj$qvalues
#cis$FDR.pool <- p.adj
cis <- merge(cis, eigen)


print("Saving cis results")
# res <- cis[cis$FDR.pool < 1, ]
write.table(cis, file = paste0(opt$output, "_All_cis_sig.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names=FALSE)

print("Saving trans results")
res <- tra[tra$FDR < 0.1, ]
write.table(tra, file = paste0(opt$output, "_All_tra_sig.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names=FALSE)

print("Plotting QQ-plot")
png(filename = paste0(opt$output, "_All_cis_sig_QQplot.png"), width = 7, height = 7, units = "in", res = 300)
ggd.qqplot(pvaluesCis)
dev.off()
