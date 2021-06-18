# WG2 pipeline: cell-type classification

This repo contains the singularity image, recipe and scripts to perform automatic 
cell-type classification using the `Seurat` label transfer based on a 
joint RNA + Protein space implemented as part of `azimuth`.

Cell type classification is supervised using a [CITE-seq](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3) reference 
containing:

- 161,764 human PBMCs
- 20,953 genes
- 224 protein markers (CITE-seq)
- 31 cell types

| **Lymphoid populations**                                        | **Myeloid populations** |
|---------------------------------------------------------------	|-----------------------	|
| CD4 Naive, CD4 TCM, CD4 TEM, CD4 CTL, Treg, CD4 Proliferating 	| CD14 Mono, CD16 Mono  	|
| CD8 Naive, CD8 TCM, CD8 TEM, CD8 Proliferating                	| cDC1, cDC2, pDC, ASDC 	|
| gdT, MAIT, ILC, dnT                                           	| HSPC                  	|
| NK, NK_CD56bright, NK Proliferating                           	| Platelet, Eryth       	|
| B naive, B intermediate, B memory, Plasmablast                	| Doublet               	|

> This PBMC reference dataset was generated as part of the Hao and Hao et al, Cell 2021 paper. It is comprised of data from eight volunteers enrolled in an HIV vaccine trial from which three time point samples were taken at day 0, 3, and 7 following vaccination. All 24 samples were processed with a CITE-seq panel of 228 TotalSeq A antibodies to generate single-cell RNA and ADT data. The data were then integrated using metholody described in the pre-print linked above to generate a weighted nearest neighbor (WNN) representation of the RNA and protein data jointly. This WNN representation is used in the Azimuth app to assign celltypes, embed in the reference UMAP, and impute protein levels for the query dataset.


# Input

Query RDS **Seurat object** containing:

- **Gene expression matrix (counts)**. [Mandatory]
- Batch information (e.g. pool / well).  [Optional]

To improve alignment to reference and cell annotation, we advise to integrate each batch to the reference separately.
Cell type classification can be processed in parallel via the `future` package.


# Output

Query RDS **Seurat object** containing:

- Cell type classification (`celltype.l2`) + prediction scores (`celltype.l2.scores`)
- Reference-based reductions:
  + sPCA
  + sUMAP
  
All results are contained in a new Seurat object with

- Two extra columns in the metadata slot (`predicted.celltype.l2` and `predicted.celltype.l2.scores`)
- Two dimensional reduction objects corresponding to the query supervised PCA and 
UMAP embeddings computed based on the sPCA. 
  
Additionally, plots for the `sPCA` and `sUMAP` reductions are included as outputs 
for exploratory evaluation.


# Running cell-type classification method

The `azimuth.sif` container has all components to perform cell-type classification.
The singularity container can be downloaded from this [dropbox](https://www.dropbox.com/s/x3ru1h02qn6oqoo/azimuth.sif?dl=0)

The container can be run using the following parameters:

- `-f` or `--file`: Query RDS Seurat object containing the gene count matrix for
the cells to classify
- `-o` or `--out`: Name of output directory containing results
- `-b` or `--batch`: Column name in Seurat object metadata. If this value is 
provided, cell-type classification is performed independently by batch. This 
approach reduces computation time and allows parallel computation. 
- `-p` or `--plan`: Parallelization approach via the `future` package:
  + **Multicore**: forked processes (recommended)
  + **Multisession**: background sessions
  + **Cluster**: external sessions on current, local, and/or remote machines
  + **Remote**: Simple access to remote sessions
- `-w` or `--workers`: Number of cores to be used for parallelization
- `-m` or `--mem`: memory size limit to be exported to each worker in _GB_. By default, 
this parameter is undefined (limitless memory export) and refers to the `future.globals.maxSize` 
setting.

Example:

```sh
singularity run -B $PWD azimuth.sif             \
                        --file query.RDS        \
                        --batch pool            \
                        --plan multicore        \
                        --workers 10
```


# Authors

- Jose Alquicira Hernandez (j.alquicira@garvan.org.au)
- Lieke Michelsen (L.C.M.Michielsen@tudelft.nl)

