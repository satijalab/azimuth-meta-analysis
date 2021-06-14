#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

ob <- readRDS(file = args[1])
ob[['umap']] <- ob[['proj.umap']]
Key(object = ob[['umap']]) <- "umap_"
DefaultAssay(object = ob[['umap']]) <- "RNA"

DefaultAssay(object = ob) <- "RNA"
ob$annotation.l1 <- ob$predicted.annotation.l1
ob$annotation.l2 <- ob$predicted.annotation.l2

ob$prediction_score_l1_bin <- sapply(X = ob$predicted.annotation.l1.score,  FUN = function(x) {
    if (x < 0.2) {
        "0-0.2"
    } else if (x < 0.4) {
        "0.2-0.4"
    } else if (x < 0.6) {
        "0.4-0.6"
    } else if (x < 0.8) {
        "0.6-0.8"
    } else {
        "0.8-1.0"
    }
})

ob$prediction_score_l2_bin <- sapply(X = ob$predicted.annotation.l2.score,  FUN = function(x) {
    if (x < 0.2) {
        "0-0.2"
    } else if (x < 0.4) {
        "0.2-0.4"
    } else if (x < 0.6) {
        "0.4-0.6"
    } else if (x < 0.8) {
        "0.6-0.8"
    } else {
        "0.8-1.0"
    }
})

ref <- readRDS(file = args[2])
DefaultAssay(ref) <- "RNA"
ref[["SCT"]] <- NULL
DefaultAssay(object = ob[['umap']]) <- "RNA"
ref$dataset_origin <- "travaglini_2020"
ref$disease <- "normal"
ref$health_status <- "healthy"
ref$donor <- ref$patient
ref$assay <- "10x v2 sequencing"
dat <- FetchData(object = ref, vars = c("XIST", "donor"))
dat %>% group_by(donor) %>% summarize(mxist = mean(XIST)) %>%
  mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
sex.labels <- as.data.frame(sex.labels)
rownames(x = sex.labels) <- sex.labels$donor
ref$sex <- sex.labels[ref$donor, "sex"]
ref$original_annotation <- ref$free_annotation
ref$prediction_score_l1_bin <- "0.8-1.0"
ref$prediction_score_l2_bin <- "0.8-1.0"

ob <- merge(ob, ref, merge.dr = "umap")

ob <- NormalizeData(object = ob)
ob <- DietSeurat(
  object = ob,
  dimreducs = "umap",
  assays = "RNA"
)

md.keep <- c("disease", "donor", "sex", "dataset_origin", "health_status", 
             "tissue", "assay", "original_annotation", "annotation.l1", 
             "annotation.l2", "predicted.annotation.l1.score", 
             "predicted.annotation.l2.score", "prediction_score_l1_bin",
             "prediction_score_l2_bin")
for (i in colnames(x = ob[[]])) {
  if (!i %in% md.keep) {
      ob[[i]] <- NULL
  }
}
Misc(object = ob[['umap']], slot = "model") <- NULL


for (i in colnames(x = ob[[]])) {
  if (is.factor(x = ob[[i, drop = TRUE]])) {
    ob[[i]] <- as.character(x = ob[[i, drop = TRUE]])
  }
}


saveRDS(object = ob, file = args[3])
SaveH5Seurat(object = ob, file = args[4], overwrite = TRUE)
Convert(args[4], dest = "h5ad", overwrite = TRUE)


