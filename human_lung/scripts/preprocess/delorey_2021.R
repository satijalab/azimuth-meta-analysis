#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(Seurat)
library(dplyr)

data.dir <- file.path("raw_data", args[1])
ob <- readRDS(file = file.path(data.dir, "SCP1052_lung.rds"))
ob <- subset(x = ob, cells = names(x = which(is.na(ob$predicted_celltype))), invert = TRUE)

hybrid.genes <- rownames(x = ob)[grep("-ENSG", x = rownames(x = ob))]
dat <- GetAssayData(object = ob[["RNA"]], slot = "counts")
for (i in hybrid.genes) {
  print(i)
  gene <- strsplit(x = i, split = "-ENSG")[[1]][1]
  if (gene %in% rownames(x = dat)) {
    dat[gene, ] <- dat[gene, ] + dat[i, ]
    dat <- dat[-c(which(rownames(x = dat) == i)), ]
  } else {
    rownames(x = dat)[rownames(x = dat) == i] <- gene
  }
}
ob <- CreateSeuratObject(counts = dat, meta.data = ob[[]])

ob$dataset_origin <- "delorey_2021"
ob$health_status <- "disease"
ob$disease <- "COVID-19"
ob$tissue <-  "lung"
ob$assay <- "10x v3 sequencing"
ob$original_annotation <- ob$predicted_celltype

ob <- NormalizeData(object = ob)
dat <- FetchData(object = ob, vars = c("XIST", "donor"))
dat %>% group_by(donor) %>% summarize(mxist = mean(XIST)) %>%
  mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
sex.labels <- as.data.frame(sex.labels)
rownames(x = sex.labels) <- sex.labels$donor
ob$sex <- sex.labels[ob$donor, "sex"]
ob.list <- SplitObject(object = ob, split.by = "donor")

d12 <- merge(ob.list[["D12_1"]], ob.list[paste0("D12_", 2:5)])
d3 <- merge(ob.list[["D3_1"]], ob.list[paste0("D3_", 2:3)])
d8 <- merge(ob.list[["D8_1"]], ob.list[paste0("D8_", 2:3)])

remove <- c(paste0("D12_", 1:5), paste0("D3_", 1:3), paste0("D8_", 1:3))
for (i in remove) {
  ob.list[[i]] <- NULL
}
ob.list <- c(ob.list, D12 = d12, D3 = d3, D8 = d8)

for(i in 1:length(x = ob.list)) {
    donors <- paste(unique(ob.list[[i]]$donor), collapse = "_")
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))
}

