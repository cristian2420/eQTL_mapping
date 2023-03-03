require(dplyr)
require(tidyr)
args <- commandArgs(trailingOnly = TRUE)

res <- list()

for(i in 1:(length(args)-1)) {
  print(paste0("Reading file ", args[i]))
  dat <- try(read.table(args[i], sep = "\t", header=TRUE))
  if(!inherits(dat, "try-error")) {
      res[[i]] <- dat
  }
}

res <- do.call(rbind, res)
res <- res[, c("gene", "statistic")]
res <- res %>% group_by(gene) %>% mutate(stat_name = paste0("stat", 1:n())) %>% pivot_wider(names_from=stat_name, values_from=statistic)
res <- as.data.frame(res)
rownames(res) <- res$gene
res$gene <- NULL
colnames(res) <- NULL

print("Saving output file")
write.table(res, file = args[length(args)], sep = "\t",
            quote = FALSE, col.names = NA, row.names = TRUE)
