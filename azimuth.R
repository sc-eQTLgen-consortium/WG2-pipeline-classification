#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Cell type classification
# author: Jose Alquicira Hernandez, Lieke Michelsen
# date: 2021-06-04
# description: Classifies cells from scRNA-seq data following Azimuth classifi-
# cation approach. 

#   ____________________________________________________________________________
#   Import libraries                                                        ####

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("SeuratDisk"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("optparse"))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option(c("-f", "--file"), 
              type = "character", 
              default = NULL, 
              help = crayon::green("RDS object file name"), 
              metavar = "character"),
  make_option(c("-o", "--out"), 
              type = "character", 
              default = "azimuth_classification", 
              help = crayon::green("output file name [default= %default]"), 
              metavar = "character"),
  make_option(c("-b", "--batch"), 
              type = "character", 
              default = NULL, 
              help = crayon::yellow("Batch column. If provided, each group in from the batch columns is mapped to reference independently"), 
              metavar = "character"),
  make_option(c("-p", "--plan"), 
              type = "character", 
              default = "sequential", 
              help = crayon::yellow("Strategy to resolve future [default= %default]:
                multicore
                multisession
                cluster
                remote
                transparent"), 
              metavar = "character"),
  make_option(c("-w", "--workers"), 
              type = "integer", 
              default = 1, 
              help = crayon::yellow("Number of workers used for parallelization
                [default= %default]"), 
              metavar = "numeric"),
  make_option(c("-m", "--mem"), 
              type = "numeric", 
              default = Inf, 
              help = crayon::yellow("Maximum allowed total size (in GB) of global variables identified
                [default= %default]"), 
              metavar = "numeric")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop(crayon::red("Query file name is missing"), call. = FALSE)
}


#   ____________________________________________________________________________
#   Define future                                                           ####


options(future.globals.maxSize = opt$mem * 1024^3)
if(opt$plan != "sequential"){
  plan(opt$plan, workers = opt$workers)
}

#   ____________________________________________________________________________
#   Helper functions                                                        ####

echo <- function(text, color = c("green", "red", "yellow", "blue")){
  
  text <- paste0(text, "\n")
  
  color <- match.arg(color)
  
  if(color == "green")
    cat(crayon::green(text))
  else if(color == "red")
    cat(crayon::red(text))
  else if(color == "yellow")
    cat(crayon::yellow(text))
  else if(color == "blue")
    cat(crayon::blue(text))
  
}

#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading query data......................................................", 
     "blue")

data <- readRDS(opt$file)
data <- UpdateSeuratObject(data)

echo("DONE....................................................................", 
     "blue")


#   ____________________________________________________________________________
#   Import CITE-seq reference                                               ####

echo("Loading CITE-seq reference..............................................", 
     "yellow")

reference <- LoadH5Seurat("/pbmc_multimodal.h5seurat")

echo("DONE....................................................................", 
     "yellow")


#   ____________________________________________________________________________
#   Split data by batch                                                     ####

if (is.null(opt$bath)){
  batches <- list(data)
}else{
  echo("Splitting data............................................................", 
       "green")
  
  batches <- SplitObject(data, opt$batch)
  
  echo("DONE......................................................................", 
       "green")
}

#   ____________________________________________________________________________
#   Normalize data                                                          ####

echo("Applying SCTransform to each batch......................................", 
     "green")

batches <- future_lapply(batches, SCTransform)

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Find anchors                                                            ####

echo("Finding anchors for each batch..........................................",
     "green")

anchors <- future_lapply(batches, function(x)
  FindTransferAnchors(reference = reference,
                      query = x,
                      normalization.method = "SCT",
                      reference.reduction = "spca",
                      dims = 1:50), 
  future.seed = TRUE)

echo("DONE....................................................................", 
     "green")

#   ____________________________________________________________________________
#   Map cells                                                               ####

echo("Finding anchors for each batch..........................................", 
     "green")

batches <- future_mapply(function(anchor, x){
  x <- MapQuery(
    anchorset = anchor,
    query = x,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap")
  
}, anchors, batches, SIMPLIFY = FALSE, future.seed = TRUE)

echo("DONE....................................................................", 
     "green")


#   ____________________________________________________________________________
#   Gather cell type classification and store in main object                ####

echo("Gathering results........................................................",
     "green")

celltype_l2 <- lapply(batches, 
                      \(x) x[[]][, 
                                 c("predicted.celltype.l2", 
                                   "predicted.celltype.l2.score"), drop = FALSE
                      ])

celltype_l2 <- lapply(celltype_l2, \(x){
  x$barcode <- row.names(x)
  x
})

celltype_l2 <- do.call(rbind, celltype_l2)
rownames(celltype_l2) <- celltype_l2$barcode
celltype_l2$barcode <- NULL
data <- AddMetaData(data, celltype_l2)

echo("DONE....................................................................", 
     "green")

#   ____________________________________________________________________________
#   Gather dimensionality reductions and store in main object               ####

echo("Storing reference-based dimensionality reductions in query object........",
     "green")

spca <- lapply(batches, \(x) x[["ref.spca"]])
umap <- lapply(batches, \(x) x[["ref.umap"]])

spca <- merge(merge(spca[[1]], spca[-1]))
umap <- merge(merge(umap[[1]], umap[-1]))

spca@assay.used <- "RNA"
umap@assay.used <- "RNA"

data[["azimuth_spca"]] <- spca
data[["azimuth_umap"]] <- umap

echo("DONE....................................................................", 
     "green")

#   ____________________________________________________________________________
#   Plot reductions                                                         ####

dir.create(opt$out)

echo("Plotting data............................................................",
     "green")

p_umap <- DimPlot(data, 
                  group.by = "predicted.celltype.l2", 
                  reduction = "azimuth_umap", 
                  label = TRUE, 
                  repel = TRUE, 
                  raster = TRUE)



p_spca <- DimPlot(data, 
                  group.by = "predicted.celltype.l2", 
                  reduction = "azimuth_spca", 
                  label = TRUE, 
                  repel = TRUE, 
                  raster = TRUE)


ggsave(file.path(opt$out, "umap.png"), p_umap, width = 10)
ggsave(file.path(opt$out, "spca.png"), p_spca, width = 10)

echo("DONE....................................................................", 
     "green")
#   ____________________________________________________________________________
#   Export data                                                             ####

echo("Saving query data........................................................",
     "yellow")

saveRDS(data, file.path(opt$out, "query.RDS"))

echo("DONE....................................................................", 
     "yellow")