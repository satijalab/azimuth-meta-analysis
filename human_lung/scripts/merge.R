#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)

all.obs <- list()
for (f in args[1:(length(x = args)-1)]) {
  ob <- readRDS(file = f)
  all.obs <- c(all.obs, ob)
}

merged.ob <- merge(x = all.obs[[1]], all.obs[2:length(all.obs)], merge.dr = "proj.umap")
print(merged.ob)
saveRDS(object = merged.ob, file = args[length(x = args)])
