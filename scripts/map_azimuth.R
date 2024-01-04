#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Cell type classification
# author: Jose Alquicira Hernandez, Lieke Michelsen, Martijn Vochteloo
# date: 2021-06-04
# description: Classifies cells from scRNA-seq data following Azimuth classifi-
# cation approach.

#   ____________________________________________________________________________
#   Import libraries                                                        ####

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SeuratDisk)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(progressr)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(tidyverse)))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <- list(
  make_option("--file",
              type = "character",
              default = NULL,
              help = crayon::green("RDS object file name"),
              metavar = "character"),
  make_option("--batch",
              type = "character",
              default = NULL,
              help = crayon::yellow("Batch column. If provided, each group in from the batch columns is mapped to reference independently"),
              metavar = "character"),
  make_option("--reference",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--refdata",
              type = "character",
              default = "celltype.l2=celltype.l2",
              help = crayon::green("Output file name [default= %default]"),
              metavar = "character"),
  make_option("--plan",
              type = "character",
              default = "sequential",
              help = crayon::yellow("Strategy to resolve future [default= %default]:
                multisession
                multicore
                cluster
                remote
                transparent"),
              metavar = "character"),
  make_option("--workers",
              type = "integer",
              default = 1,
              help = crayon::yellow("Number of workers used for parallelization
                [default= %default]"),
              metavar = "numeric"),
  make_option("--mem",
              type = "numeric",
              default = Inf,
              help = crayon::yellow("Maximum allowed total size (in GB) of global variables identified
                [default= %default]"),
              metavar = "numeric"),
  make_option(c("--palette"),
              type = "character",
              default = NULL,
              help = crayon::green("The palette to use"),
              metavar = "character"),
  make_option("--out",
              type = "character",
              default = "azimuth",
              help = crayon::green("Output file name [default= %default]"),
              metavar = "character"),
  make_option("--path",
              type = "character",
              default = ".",
              help = crayon::green("Output path to store results [default= %default]"),
              metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop(crayon::red("Query file name is missing"), call. = FALSE)
}

#   ____________________________________________________________________________
#   Helper functions                                                        ####

echo <- function(text, color = c("green", "red", "yellow", "blue")) {

  text <- paste0(text, "\n")

  color <- match.arg(color)

  if (color == "green")
    cat(crayon::green(text))
  else if (color == "red")
    cat(crayon::red(text))
  else if (color == "yellow")
    cat(crayon::yellow(text))
  else if (color == "blue")
    cat(crayon::blue(text))

}

#   ____________________________________________________________________________
#   Input information                                                       ####

echo("Input information.......................................................",
     "yellow")

echo(paste0(crayon::bold("Input file:\n"), opt$file), "yellow")
echo(paste0(crayon::bold("Batch variable:\n"), opt$batch), "yellow")
echo(paste0(crayon::bold("Reference file:\n"), opt$reference), "yellow")
echo(paste0(crayon::bold("Reference data:\n"), opt$refdata), "yellow")
echo(paste0(crayon::bold("Parallelization plan:\n"), opt$plan), "yellow")
echo(paste0(crayon::bold("Number of workers: \n"), opt$workers), "yellow")
echo(paste0(crayon::bold("maxSize future global: \n"), opt$mem), "yellow")
echo(paste0(crayon::bold("Palette: "), opt$palette), "yellow")
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory: "), opt$path), "yellow")


echo("DONE....................................................................",
     "yellow")

#   ____________________________________________________________________________
#   Define future                                                           ####

echo("Future settings.........................................................",
     "blue")

handlers(global = TRUE)
handlers("progress")


if (opt$plan != "sequential") {
  options(future.globals.maxSize=(args$mem * 1000 * 1024^2))
  plan(opt$plan, workers = opt$workers)
}

echo("DONE....................................................................",
     "blue")

#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading query data......................................................",
     "blue")

data <- readRDS(opt$file)
if (!inherits(data, "Seurat")) stop("Input query data is not a Seurat object")
data <- UpdateSeuratObject(data)

echo("DONE....................................................................",
     "blue")


#   ____________________________________________________________________________
#   Import seq reference                                               ####

echo("Loading reference..............................................",
     "yellow")

reference <- tryCatch({
	print("Loading reference using readRDS()")
	reference <- readRDS(opt$reference)
},error = function(e){
	print("Failed, trying to load reference using LoadH5Seurat()")
	reference <- LoadH5Seurat(opt$reference)
})

echo("DONE....................................................................",
     "yellow")


#   ____________________________________________________________________________
#   Split data by batch                                                     ####

if (is.null(opt$batch)) {
  batches <- list(data)
}else {
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

apply_sctransform <- function(xs) {
  p <- progressor(along = xs)
  mapply(function(x, i, n) {
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

find_anchors <- function(xs) {
  p <- progressor(along = xs)

  lapply(xs, function(x) {

    x <- FindTransferAnchors(reference = reference,
                             query = x,
                             normalization.method = "SCT",
                             reference.reduction = "pca",
                             dims = 1:50)
    p()
    x

  })

}

anchors <- find_anchors(batches)

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Map cells                                                               ####

echo("Mapping query cells for each batch......................................",
     "green")

refdata <- list()
for (refcolumn in unlist(strsplit(opt$refdata, ";"))) {
  keyvalue <- unlist(strsplit(refcolumn, "="))
  refdata[[keyvalue[1]]] <- keyvalue[[2]]
}

map_cells <- function(anchors, batches) {
  p <- progressor(along = anchors)

  mapply(function(anchor, x) {
    x <- MapQuery(
      anchorset = anchor,
      query = x,
      reference = reference,
      refdata = refdata,
      reference.reduction = "pca",
      reduction.model = "umap")

  }, anchors, batches, SIMPLIFY = FALSE)

}

batches <- map_cells(anchors, batches)

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Gather cell type classification and store in main object                ####

echo("Gathering results........................................................",
     "green")

metadata.columns <- c(paste0("predicted.", names(refdata)), paste0("predicted.", names(refdata), ".score"))
if (!is.null(opt$batch)) {
  metadata.columns <- c(opt$batch, metadata.columns)
}

metadata <- lapply(batches, function(x) x[[]][, metadata.columns, drop = FALSE])
metadata <- lapply(metadata, function(x) {
  x$barcode <- row.names(x)
  x
})

metadata <- do.call(rbind, metadata)
metadata <- metadata[, c("Pool", "barcode", metadata.columns)]
write_delim(metadata, file=gzfile(paste0(opt$path, opt$out, ".metadata.tsv.gz")), delim="\t")

rownames(metadata) <- metadata$barcode
metadata$barcode <- NULL
data <- AddMetaData(data, metadata)

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Gather dimensionality reductions and store in main object               ####

echo("Storing reference-based dimensionality reductions in query object.......",
     "green")

pca <- lapply(batches, function(x) x[["ref.pca"]])
umap <- lapply(batches, function(x) x[["ref.umap"]])

pca <- merge(merge(pca[[1]], pca[-1]))
umap <- merge(merge(umap[[1]], umap[-1]))

pca@assay.used <- "RNA"
umap@assay.used <- "RNA"

data[["azimuth_pca"]] <- pca
data[["azimuth_umap"]] <- umap

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Plot reductions                                                          ####


echo("Plotting data...........................................................",
     "green")

palette <- NULL
if(!is.null(opt$palette)){
  palette <- c()
  palette_names <- c()
  for (ct in unlist(strsplit(opt$palette, ";"))) {
    keyvalue <- unlist(strsplit(ct, "="))
    palette <- c(palette, keyvalue[[2]])
    palette_names <- c(palette_names, keyvalue[[1]])
  }
  names(palette) <- palette_names
}

for (dim_reduction in c("pca", "umap")) {
  plot <- NULL
  for (group_by in c("Pool", paste0("predicted.", names(refdata)))) {
    cols <- palette
    for (value in unique(data[[group_by]])) {
      if (!value %in% names(palette)) {
        cols <- NULL
      }
    }

    p <- DimPlot(data,
                 group.by = group_by,
                 cols = cols,
                 reduction = paste0("azimuth_", dim_reduction),
                 label = TRUE,
                 repel = TRUE,
                 raster = TRUE) + NoLegend()
    if (is.null(plot)) {
      plot <- p
    } else {
      plot <- plot + p
    }
  }
  ncols <- ceiling(length(plot) / 2)
  nrows <- ceiling(length(plot) / ncols)
  ggsave(paste0(opt$path, opt$out, "_ref_", dim_reduction, ".png"),
         plot,
         width = min(10 * ncols, 49),
         height = min(10 * nrows, 49))
}

echo("DONE....................................................................",
     "green")
#   ____________________________________________________________________________
#   Export data                                                             ####

# echo("Saving query data.......................................................",
#      "yellow")
#
# saveRDS(data, paste0(opt$path, opt$out, ".RDS"))
#
# echo("DONE....................................................................",
#      "yellow")