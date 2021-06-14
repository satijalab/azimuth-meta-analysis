# Meta-analysis of community datasets

This repository contains code to map publicly available single-cell datasets to Azimuth references. 

# Overview

The mapping workflow is setup as a snakemake workflow. This can be run locally with `snakemake --use-singularity --cores X RULE`. 

## Workflow

The main steps of the workflow are as follows:

1. Download Azimuth references and query data (urls provided in the `links` directory)
2. Preprocess newly downloaded query data and save into single Seurat objects per donor (if > 1000 cells).
3. Map each donor to the reference individually
4. Merge results from each donor back into a single Seurat object per dataset

The mapped results can be found as Seurat objects in the `seurat_objects/mapped` directory.

## Adding new query datasets to be mapped

In order to add new datasets into the workflow, you will need to:

1. Decide on a name for the new query dataset. For existing datasets, this is done by the author and year of publication (e.g. `reyfman_2019`).
2. Add the download links to a plain text file in the `links/query` directory with your chosen name.
3. Add a preprocessing script to `scripts/preprocess`. This should be written as an R script that accepts command line arguments where the first argument is the name of the directory in `raw_data` containing the downloaded data and the second argument is the output directory. This script is expected to create individual Seurat objects for each donor in the dataset and save them as `.rds` files. It should also annotate the following metadata fields when possible: 

| Field | Description |
| ----- | ----------- |
| dataset_origin | Unique identifier for the dataset |
| donor | Identifier for the Donor |
| health_status | Specify either "disease" or "healthy" |
| disease | Specific type of lung disease (unless "normal") |
| assay | Single-cell technology used in data collection |
| tissue | Source of the tissue sample |
| sex | Sex of the donor, can be determined by XIST expression if not provided |
| original_annotation | Original cell type annotation |

