#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

data.dir <- file.path("raw_data", args[1])

file <- file.path(data.dir, "munich_cohort_human_dataset_updated.h5ad")

Convert(file, dest = "h5seurat", overwrite = TRUE)
fname <- gsub(pattern = ".h5ad", replacement = ".h5seurat", x = file)
ob <- LoadH5Seurat(
  file = fname,
  assays = list("counts" = c("counts", "data")),
  reductions = FALSE
)
ob[["RNA"]] <- ob[["counts"]]
DefaultAssay(ob) <- "RNA"
ob[["counts"]] <- NULL
ob[["RNA"]] <- SetAssayData(
  object = ob[["RNA"]],
  slot = "counts",
  new.data = GetAssayData(object = ob[["RNA"]], slot = "data")
)

ob$donor <- ob$patient_id
ob <- NormalizeData(object = ob)

dat <- FetchData(object = ob, vars = c("XIST", "donor"))
dat %>% group_by(donor) %>% summarize(mxist = mean(XIST)) %>%
  mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
sex.labels <- as.data.frame(sex.labels)
rownames(x = sex.labels) <- sex.labels$donor
ob$sex <- sex.labels[ob$donor, "sex"]

ob$dataset_origin <- "mayr_2020"
ob$disease <- "interstitial lung disease"
ob$health_status <- as.character(x = ob$health_status)
ob$disease[ob$health_status == "control donor"] <- "normal"
ob$health_status[ob$health_status == "control donor"] <- "healthy"
ob$health_status[ob$health_status != "healthy"] <- "disease"
ob$assay <- "Drop-seq"
ob$tissue <- "lung"
ob$original_annotation <- ob$cell_state_label

ob.list <- SplitObject(object = ob, split.by = "donor")
small.obs <- ob.list[names(which(table(ob$donor) < 1000))]
for (o in names(ob.list)){
  if(o %in% names(x = small.obs)) {
    ob.list[[o]] <- NULL
  }
}
small.obs.merge <- merge(x = small.obs[[1]], y = small.obs[2:length(x = small.obs)])
ob.list[[1]] <- merge(small.obs.merge, ob.list[[1]])

for(i in 1:length(x = ob.list)) {                                                                      
    donors <- paste(unique(ob.list[[i]]$donor), collapse = "_")                                        
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))                  
}  
