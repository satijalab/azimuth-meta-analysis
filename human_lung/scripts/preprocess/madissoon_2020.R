#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

data.dir <- file.path("raw_data", args[1])                                                          
dat <- readRDS(file = file.path(data.dir, "lung_ts.rds"))
ob <- CreateSeuratObject(counts = GetAssayData(object = dat[["RNA"]], slot = "counts"), meta.data = dat[[]])
ob <- NormalizeData(object = ob)

dat <- FetchData(object = ob, vars = c("XIST", "Donor"))
dat %>% group_by(Donor) %>% summarize(mxist = mean(XIST)) %>%
  mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
sex.labels <- as.data.frame(sex.labels)
rownames(x = sex.labels) <- sex.labels$Donor
ob$sex <- sex.labels[ob$Donor, "sex"]
ob$donor <- ob$Donor

ob$dataset_origin <- "madissoon_2020"
ob$health_status <- "healthy"
ob$disease <- "normal"

ob$tissue <- "lung"
ob$assay <- "10x v2 sequencing"
ob$original_annotation <- ob$Celltypes_updated_July_2020

ob.list <- SplitObject(object = ob, split.by = "donor")
for(i in 1:length(x = ob.list)) {                                                                      
    donors <- paste(unique(ob.list[[i]]$Donor), collapse = "_")                                        
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))
}  
