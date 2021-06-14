#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(Seurat)
library(dplyr)

data.dir <- file.path("raw_data", args[1])
mat <- readMM(file = file.path(data.dir, "GSE136831_RawCounts_Sparse.mtx.gz"))
cells <- readLines(con =  file.path(data.dir, "GSE136831_AllCells.cellBarcodes.txt.gz"))
features <- read.table(file = file.path(data.dir, "GSE136831_AllCells.GeneIDs.txt.gz"), header = TRUE)
# biomaRt was being unreliable
# mart <- useDataset(
#     "hsapiens_gene_ensembl",
#     useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = 'ensembl.org')
# )
#
# genes <- grep(pattern = "ENSG", x = features[, 2], value = TRUE)
# gene_IDs <- rbind(
#   getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = genes[1:5000], mart = mart),
#   getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = genes[5001:10000], mart = mart),
#   getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = genes[10001:length(x = genes)], mart = mart)
# )
gene_IDs <- read.csv(file = file.path(data.dir, "adams_gene_ids.csv"))

rownames(x = mat) <- features[, 2]
colnames(x = mat) <- cells
meta <- read.table(file = file.path(data.dir, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"), row.names = 1, head = TRUE)
mat <- as(mat, "dgCMatrix")

for (i in 1:nrow(x = gene_IDs)) {
  if (gene_IDs[i, "hgnc_symbol"] != "") {
    rn.idx <- which(x = rownames(x = mat) == gene_IDs[i, "ensembl_gene_id"])
    if (gene_IDs[i, "hgnc_symbol"] %in% rownames(x = mat)) {
      mat[gene_IDs[i, "hgnc_symbol"], ] <- mat[gene_IDs[i, "hgnc_symbol"], ] + mat[gene_IDs[i, "ensembl_gene_id"], ]
      mat <- mat[-c(rn.idx), ]
    } else {
      rownames(x = mat)[rn.idx] <- gene_IDs[i, "hgnc_symbol"]
    }
  }
}

gc()

ob <- CreateSeuratObject(
  counts = mat,
  meta.data = meta,
)
ob$dataset_origin <- "adams_2020"
ob$donor <- ob$Subject_Identity

ob$disease <- ob$Disease_Identity
# remap Disease_Identity to common ontology
ob$disease[ob$disease == "Control"] <- "normal"
ob$disease[ob$disease == "COPD"] <- "chronic obstructive pulmonary disease"
ob$disease[ob$disease == "IPF"] <- "idiopathic pulmonary fibrosis"

ob$health_status <- "disease"
ob$health_status[ob$Disease_Identity == "Control"] <- "healthy"


ob$assay <- "10x v2 sequencing"
ob$tissue <- "lung"
ob$original_annotation <- ob$Subclass_Cell_Identity
rm(mat)
gc()

ob <- NormalizeData(object = ob)
dat <- FetchData(object = ob, vars = c("XIST", "donor"))
dat %>% group_by(donor) %>% summarize(mxist = mean(XIST)) %>%
  mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
sex.labels <- as.data.frame(sex.labels)
rownames(x = sex.labels) <- sex.labels$donor
ob$sex <- sex.labels[ob$donor, "sex"]

ob.list <- SplitObject(object = ob, split.by = "Subject_Identity")

small.obs <- ob.list[names(which(table(ob$Subject_Identity) < 1000))]
for (o in names(ob.list)){
  if(o %in% names(x = small.obs)) {
    ob.list[[o]] <- NULL
  }
}
small.obs.merge <- merge(x = small.obs[[1]], y = small.obs[2:length(x = small.obs)])
ob.list[["small"]] <- small.obs.merge

for(i in 1:length(x = ob.list)) {
    donors <- paste(unique(ob.list[[i]]$donor), collapse = "_")
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))
}
