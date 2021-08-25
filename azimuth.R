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
suppressPackageStartupMessages(library("progressr"))
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
                multisession
                multicore
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
#   Input information                                                       ####

echo("Input information.......................................................", 
     "yellow")

echo(paste0(crayon::bold("Input file:\n"), opt$file), "yellow")
echo(paste0(crayon::bold("Output directory for results:\n"), opt$out), "yellow")
echo(paste0(crayon::bold("Batch variable:\n"), opt$batch), "yellow")
echo(paste0(crayon::bold("Parallelization plan:\n"), opt$plan), "yellow")
echo(paste0(crayon::bold("Number of workers: \n"), opt$workers), "yellow")
echo(paste0(crayon::bold("maxSize future global: \n"), opt$mem), "yellow")


echo("DONE....................................................................", 
     "yellow")

#   ____________________________________________________________________________
#   Define future                                                           ####

echo("Future settings.........................................................", 
     "blue")

options(future.globals.maxSize = opt$mem * 1024^3)
handlers(global = TRUE)
handlers("progress")


if(opt$plan != "sequential"){
  plan(opt$plan, workers = opt$workers)
}

echo("DONE....................................................................", 
     "blue")

#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading query data......................................................", 
     "blue")

data <- readRDS(opt$file)
if(!inherits(data, "Seurat")) stop("Input query data is not a Seurat object")
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

if (is.null(opt$batch)){
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

apply_sctransform <- function(xs){
  p <- progressor(along = xs)
  future_lapply(xs, function(x){
    x <- SCTransform(x, verbose = FALSE)
    p()
    x
  })
}

batches <- apply_sctransform(batches)

echo("DONE....................................................................",
     "green")

###

echo("Applying SCTransform to each batch......................................", 
     "green")

apply_sctransform <- function(xs){
  p <- progressor(along = xs)
  future_lapply(xs, function(x, i, n){
    x <- SCTransform(x, verbose = FALSE)
    p(message = sprintf("| Batch %g", x))
    x
  }, length(xs))
}

batches <- apply_sctransform(batches)

echo("DONE....................................................................",
     "green")


###

echo("Applying SCTransform to each batch......................................", 
     "green")

apply_sctransform <- function(xs){
  p <- progressor(along = xs)
  future_imap(xs, function(x, i, n){
    x <- SCTransform(x, verbose = FALSE)
    p(message = sprintf("| Batch %d/%d", i, n))
    x
  }, length(xs))
}

batches <- apply_sctransform(batches)

echo("DONE....................................................................",
     "green")

####

echo("Applying SCTransform to each batch......................................", 
     "green")

handlers(global = TRUE)
handlers("progress")
batches <- readRDS("batches.RDS")

apply_sctransform <- function(xs){
  p <- progressor(along = xs)
   mapply(function(x, i, n){
    x <- SCTransform(x, verbose = FALSE)
    p(message = sprintf("| Batch %d/%d", i, n))
    x
  }, xs, seq_along(xs), MoreArgs = list(n = length(xs)), SIMPLIFY = FALSE)
}

batches <- apply_sctransform(batches)

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Find anchors                                                            ####

echo("Finding anchors for each batch..........................................",
     "green")

find_anchors <- function(xs){
  p <- progressor(along = xs)
  
  future_lapply(xs, function(x){
    
    x <- FindTransferAnchors(reference = reference,
                             query = x,
                             normalization.method = "SCT",
                             reference.reduction = "spca",
                             dims = 1:50) 
    p()
    x
    
  }, future.seed = TRUE)
  
}

anchors <- find_anchors(batches)

echo("DONE....................................................................", 
     "green")

#   ____________________________________________________________________________
#   Map cells                                                               ####

echo("Mapping query cells for each batch......................................", 
     "green")

map_cells <- function(anchors, batches){
  p <- progressor(along = anchors)
  
  future_mapply(function(anchor, x){
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
  
}

batches <- map_cells(anchors, batches)

echo("DONE....................................................................", 
     "green")

#   ____________________________________________________________________________
#   Gather cell type classification and store in main object                ####

echo("Gathering results........................................................",
     "green")

celltype_l2 <- lapply(batches, 
                      function(x) x[[]][, 
                                        c("predicted.celltype.l2", 
                                          "predicted.celltype.l2.score"), drop = FALSE
                      ])

celltype_l2 <- lapply(celltype_l2, function(x){
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

echo("Storing reference-based dimensionality reductions in query object.......",
     "green")

spca <- lapply(batches, function(x) x[["ref.spca"]])
umap <- lapply(batches, function(x) x[["ref.umap"]])

spca <- merge(merge(spca[[1]], spca[-1]))
umap <- merge(merge(umap[[1]], umap[-1]))

spca@assay.used <- "RNA"
umap@assay.used <- "RNA"

data[["azimuth_spca"]] <- spca
data[["azimuth_umap"]] <- umap

echo("DONE....................................................................", 
     "green")


#   ____________________________________________________________________________
#   Gather imputed ADT                                                      ####

echo("Storing imputed ADT assay in query object...............................",
     "green")

predicted_ADT <- lapply(batches, function(x) x[["predicted_ADT"]])

if(length(predicted_ADT) == 1){
  predicted_ADT <- predicted_ADT[[1]]
}else{
  predicted_ADT <- merge(predicted_ADT[[1]], predicted_ADT[-1])
}

data[["predicted_ADT"]] <- predicted_ADT

echo("DONE....................................................................", 
     "green")

#   ____________________________________________________________________________
#   Plot reductions                                                         ####

dir.create(opt$out)

echo("Plotting data...........................................................",
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

echo("Saving query data.......................................................",
     "yellow")

saveRDS(data, file.path(opt$out, "query.RDS"))

echo("DONE....................................................................", 
     "yellow")