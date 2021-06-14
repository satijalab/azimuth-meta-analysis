#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)

data.dir <- file.path("raw_data", args[1])

system2(command = "mv", args = paste0(data.dir, "/2 ", data.dir, "/2.zip"))
system2(command = "unzip", args = paste0(file.path(data.dir, "2.zip"), " -d ", data.dir))

lung.dat <- fread(file = file.path(data.dir, "Counts_lung_cells.csv"), sep = ",")
lung.dat <- as(object = as.matrix(x = lung.dat, rownames = 1), Class = "dgCMatrix")
lung_meta <- read.csv(file = file.path(data.dir, "Metadata_lung_cells_update.csv"), header = TRUE, row.names = 1)
rownames(x = lung_meta) <- make.names(names = rownames(x = lung_meta))

lung <- CreateSeuratObject(counts = lung.dat, meta.data = lung_meta)

lung$sex <- lung$Sex
lung$sex[lung$sex == "F"] <- "female"
lung$sex[lung$sex == "M"] <- "male"
lung$donor <- lung$ID

lung$original_annotation <- lung$Cell.type
lung$health_status <- "healthy"

lung$disease <- "normal"
lung$assay <- "10x v2 sequencing"
lung$tissue <- "lung"
lung$dataset_origin <- "lukassen_2020"
lung.list <- SplitObject(object = lung, split.by = "donor")

hbec.dat <- fread(file = file.path(data.dir, "Counts_HBECs.csv"), sep = ",")
hbec.dat <- as(object = as.matrix(x = hbec.dat, rownames = 1), Class = "dgCMatrix")
hbec_meta <- read.csv(file = file.path(data.dir, "Metadata_HBECs.csv"), header = TRUE, row.names = 1)

hbec <- CreateSeuratObject(counts = hbec.dat, meta.data = hbec_meta)
hbec$sex <- hbec$Sex
hbec$sex[hbec$sex == "f"] <- "female"
hbec$sex[hbec$sex == "m"] <- "male"
hbec$donor <- hbec$ID
hbec$original_annotation <- hbec$Cell.type
hbec$health_status <- "healthy"
hbec$disease <- "normal"
hbec$assay <- "10x v2 sequencing"
hbec$tissue <- "epithelium of bronchus"

hbec$dataset_origin <- "lukassen_2020"
hbec$donor[is.na(x = hbec$donor)] <- sapply(X = names(hbec$donor[is.na(x = hbec$donor)]), FUN = Seurat:::ExtractField )
hbec.list <- SplitObject(object = hbec, split.by = "donor")
ob.list <- c(lung.list, hbec.list)

# some small subset of cells have NA for sex - relabel by donor
ob.list <- lapply(X = ob.list, FUN = function(ob) {
    ob$sex <- names(x = sort(x = table(ob$sex), decreasing = TRUE))[1]
    return(ob)
})

for(i in 1:length(x = ob.list)) {
    donors <- paste(unique(ob.list[[i]]$donor), collapse = "_")
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))
} 

