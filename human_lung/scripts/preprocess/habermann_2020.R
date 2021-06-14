#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(dplyr)

data.dir <- file.path("raw_data", args[1])
dat <- ReadMtx(
  mtx = file.path(data.dir, "GSE135893_matrix.mtx.gz"),
  cells = file.path(data.dir, "GSE135893_barcodes.tsv.gz"),
  features = file.path(data.dir, "GSE135893_genes.tsv.gz"),
  feature.column = 1
)

meta <- read.csv(file = file.path(data.dir, "GSE135893_IPF_metadata.csv.gz"), row.names = 1)
metacols.keep <- c("orig.ident", "Diagnosis", "Sample_Name", "Sample_Source", "Status", "population", "celltype")
for (i in colnames(x = meta)) {
  if (! i %in% metacols.keep) {
    meta[[i]] <- NULL
  }
}

ob <- CreateSeuratObject(
  counts = dat,
  meta.data = meta,
)
ob$dataset_origin <- "habermann_2020"
ob$donor <- ob$Sample_Name
ob <- subset(x = ob, cells = names(which(is.na(ob$celltype))), invert = TRUE)
ob <- NormalizeData(object = ob)

dat <- FetchData(object = ob, vars = c("XIST", "donor"))
dat %>% group_by(donor) %>% summarize(mxist = mean(XIST)) %>%
  mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
sex.labels <- as.data.frame(sex.labels)
rownames(x = sex.labels) <- sex.labels$donor
ob$sex <- sex.labels[ob$donor, "sex"]

ob$disease <- ob$Diagnosis
ob$disease[ob$disease == "cHP"] <- "hypersensitivity pneumonitis"
ob$disease[ob$disease == "Control"] <- "normal"
ob$disease[ob$disease == "ILD"] <- "interstitial lung disease"
ob$disease[ob$disease == "IPF"] <- "idiopathic pulmonary fibrosis"
ob$disease[ob$disease == "NSIP"] <- "non-specific interstitial pneumonia"
ob$disease[ob$disease == "Sarcoidosis"] <- "sarcoidosis"
  

ob$health_status <- "disease"
ob$health_status[ob$disease == "normal"] <- "healthy"

ob$tissue <- "lung"
ob$assay <- "10x sequencing"
ob$original_annotation <- ob$celltype

ob.list <- SplitObject(object = ob, split.by = "Sample_Name")
small.obs <- ob.list[names(which(table(ob$Sample_Name) < 1000))]
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

