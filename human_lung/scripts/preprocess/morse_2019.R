#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)

data.dir <- file.path("raw_data", args[1])
system2(command = "tar", args = paste0("-xvf ", data.dir, "/'index.html?acc=GSE128033&format=file' -C ", data.dir))

donors <- list.files(path = data.dir, pattern = "matrix.mtx.gz")
files <- gsub(pattern = "matrix.mtx.gz", replacement = "", x = donors)
ob.list <- list()

for (f in files) {
  fpath <- file.path(data.dir, f)
  dat <- ReadMtx(mtx = paste0(fpath, "matrix.mtx.gz"), 
                 cells = paste0(fpath, "barcodes.tsv.gz"),
                 features = paste0(fpath, "genes.tsv.gz")
  )
  print(f)
  # match filtering done according to GEO entry
  ob <- CreateSeuratObject(counts = dat, min.cells = 1, min.features = 200)
  ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^MT-")
  ob <- subset(x = ob, subset = percent.mt < 35 & nCount_RNA < 50000)
  ob <- NormalizeData(object = ob)
  
  if ("XIST" %in% rownames(x = ob)) {
    dat <- FetchData(object = ob, vars = "XIST")
    dat %>% summarize(mxist = mean(XIST)) %>%
      mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
    ob$sex <- sex.labels$sex
  } else {
    ob$sex <- "male"
  }
  donor <- Seurat:::ExtractField(string = f, field = 2)
  ob$donor <- donor
  ob$health_status <- ifelse(test = grepl(pattern = "NOR", x = donor), yes = "healthy", no = "disease")
  ob$disease <- ob$health_status
  ob$disease[ob$health_status == "disease"] <- "idiopathic pulmonary fibrosis"
  ob$disease[ob$health_status == "healthy"] <- "normal"
  ob$assay <- "10x v2 sequencing"
  ob$tissue <- ifelse(test = grepl(pattern = "bal", x = donor), yes = "bronchial alveolar lavage", "lung")
  ob$original_annotation <- NA
  ob$dataset_origin <- "morse_2019"
  ob.list <- c(ob.list, ob)
}

for(i in 1:length(x = ob.list)) {                                                                      
    donors <- paste(unique(ob.list[[i]]$donor), collapse = "_")                                        
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(donors, ".rds")))                  
}  
