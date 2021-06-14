#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(dplyr)

data.dir <- file.path("raw_data", args[1])

system2(command = "tar", args = paste0("-xvf ", data.dir, "/'index.html?acc=GSE122960&format=file' -C ", data.dir))

in.files <- list.files(data.dir, pattern = "filtered_gene")
ob.list <- list()

for (i in 1:length(x = in.files)) {
    fname <- in.files[i]
    meta <- unlist(strsplit(x = fname, split = "_"))
    dat <- Read10X_h5(filename = file.path(data.dir, in.files[i]))
    ob.list[[i]] <- CreateSeuratObject(counts = dat)
    ob.list[[i]]$sample <- paste(strsplit(x = in.files[i], split = "_")[[1]][2:3], collapse = "_")
    ob.list[[i]] <- NormalizeData(object = ob.list[[i]])
    ob.list[[i]]$dataset_origin <- "reyfman_2019"
    ob.list[[i]]$donor <- paste(meta[2], meta[3], sep = "_")
    
    disease.dict <- c("Cryobiopsy" = "interstitial lung disease", 
        "Donor" = "normal", "HP" = "hypersensitivity pneumonitis", 
        "IPF" = "idiopathic pulmonary fibrosis", 
        "Myositis-ILD" = "interstitial lung disease", 
        "SSc-ILD" = "interstitial lung disease"
    )
    ob.list[[i]]$disease <- unname(obj = disease.dict[meta[2]])
    ob.list[[i]]$health_status <- "disease"
    ob.list[[i]]$health_status[ob.list[[i]]$disease == "normal"] <- "healthy"

    dat <- FetchData(object = ob.list[[i]], vars = c("XIST", "donor"))
    dat %>% group_by(donor) %>% summarize(mxist = mean(XIST)) %>%
        mutate(sex = ifelse(test = mxist > 0.1, yes = "female", no = "male")) -> sex.labels
    sex.labels <- as.data.frame(sex.labels)
    rownames(x = sex.labels) <- sex.labels$donor
    ob.list[[i]]$sex <- sex.labels[ob.list[[i]]$donor, "sex"]
    ob.list[[i]]$assay <- '10x v2 sequencing'
    ob.list[[i]]$tissue <- "lung"
    saveRDS(object = ob.list[[i]], file = file.path(args[2], paste0(ob.list[[i]]$donor[1], ".rds")))
}


