#!/usr/bin/env Rscript
library(Seurat)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
data.dir <- file.path("raw_data", args[1])

dat <- fread(input = file.path(data.dir, "Raw_exprMatrix.tsv.gz"))
dat <- as(object = as.matrix(x = dat, rownames = 1), Class = "dgCMatrix")
meta <- read.table(file = file.path(data.dir, "meta.tsv"), header = TRUE, sep = "\t", row.names = 1)

ob <- CreateSeuratObject(counts = dat, meta.data = meta, min.cells = 1, min.features = 1)
ob$sex <- ob$Sex
ob$dataset_origin <- "deprez_2020"
ob$health_status <- "healthy"
ob$disease <- "normal"
ob$assay <- "10x v2 sequencing"
ob$tissue <- "respiratory tract epithelium" 
ob$original_annotation <- ob$CellType
ob$donor <- ob$Donor

ob.list <- SplitObject(object = ob, split.by = "donor")
for(i in 1:length(x = ob.list)) {                                                                   
    donors <- paste(unique(ob.list[[i]]$Donor), collapse = "_")
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))               
}  
