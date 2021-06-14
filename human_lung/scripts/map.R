#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(Azimuth)

reference <- Azimuth:::LoadReference(path = args[2])
query.all <- readRDS(file = args[1])
ref.xfer <- args[3:(length(args)-1)]

if (inherits(query.all, "Seurat")) {
    query.all <- list(query.all)
}

query.mapped <- lapply(X = 1:length(x = query.all), FUN = function(x) {
    message("Query: ", x )
    query <- query.all[[x]]
    # Preprocess with SCTransform
    query <- SCTransform(
      object = query,
      assay = "RNA",
      new.assay.name = "refAssay",
      residual.features = rownames(x = reference$map),
      reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
      method = 'glmGamPoi',
      ncells = 2000,
      n_genes = 2000,
      do.correct.umi = FALSE,
      do.scale = FALSE,
      do.center = TRUE
    )

    # Find anchors between query and reference
    anchors <- FindTransferAnchors(
      reference = reference$map,
      query = query,
      k.filter = NA,
      reference.neighbors = "refdr.annoy.neighbors",
      reference.assay = "refAssay",
      query.assay = "refAssay",
      reference.reduction = "refDR",
      normalization.method = "SCT",
      features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
      dims = 1:50,
      n.trees = 20,
      mapping.score.k = 100
    )

    # Transfer cell type labels 
    refdata <- lapply(X = ref.xfer, function(x) {
      reference$map[[x, drop = TRUE]]
    })
    names(x = refdata) <- ref.xfer
    query <- TransferData(
      reference = reference$map,
      query = query,
      dims = 1:50,
      anchorset = anchors,
      refdata = refdata,
      n.trees = 20,
      store.weights = TRUE
    )
    query <- IntegrateEmbeddings(
      anchorset = anchors,
      reference = reference$map,
      query = query,
      reductions = "pcaproject",
      reuse.weights.matrix = TRUE
    )
    query[["query_ref.nn"]] <- FindNeighbors(
      object = Embeddings(reference$map[["refDR"]]),
      query = Embeddings(query[["integrated_dr"]]),
      return.neighbor = TRUE,
      l2.norm = TRUE
    )
    query <- Azimuth:::NNTransform(
      object = query,
      meta.data = reference$map[[]]
    )
    # Project the query to the reference UMAP.
    query[["proj.umap"]] <- RunUMAP(
      object = query[["query_ref.nn"]],
      reduction.model = reference$map[["refUMAP"]],
      reduction.key = 'UMAP_'
    )
    # Calculate mapping score and add to metadata
    query <- AddMetaData(
      object = query,
      metadata = MappingScore(anchors = anchors),
      col.name = "mapping.score"
    )
    DefaultAssay(object = query) <- "RNA"
    query[["refAssay"]] <- NULL
    gc()
    return(query)
})

saveRDS(object = query.mapped, file = args[length(x = args)])

