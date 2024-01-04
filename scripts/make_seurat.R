#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Make seurat object
# author: Martijn Vochteloo
# date: 2023-12-14
# description:

#   ____________________________________________________________________________
#   Import libraries                                                        ####

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(progressr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(scCustomize)))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option("--poolsheet",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--batch",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--out",
              type = "character",
              default = "split_object",
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

if (is.null(opt$poolsheet)) {
  print_help(opt_parser)
  stop(crayon::red("Poolsheet file name is missing"), call. = FALSE)
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

echo(paste0(crayon::bold("Poolsheet file:\n"), opt$poolsheet), "yellow")
echo(paste0(crayon::bold("Batch variable:\n"), opt$batch), "yellow")
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory: "), opt$path), "yellow")


echo("DONE....................................................................",
     "yellow")

#   ____________________________________________________________________________
#   Create Seurat object

echo("Creating seurat object..............................................",
     "yellow")

pools <- read_delim(opt$poolsheet, delim = "\t")
if (!is.null(opt$batch)) {
  pools <- pools %>% filter(Pool == opt$batch)
}

load_seurat_object <- function(row){
  print(paste0("  Creating Seurat object for pool ", row[["Pool"]]))
  counts <- tryCatch({
      Read10X_h5(row[["Counts"]])
  }, error = function(e) {
      Read_CellBender_h5_Mat(row[["Counts"]])
  })

  meta.data <- as.data.frame(read_delim(row[["Barcodes"]], delim="\t", col_names=c("Barcode")))
  print(paste0("    Loaded assignments with shape (", length(rownames(meta.data)), ", ", length(colnames(meta.data)), ")"))

  if (!identical(colnames(counts), meta.data$Barcode)) {
    print("Error, barcodes do not match count matrix.")
    quit()
  }

  meta.data$Pool <- row[["Pool"]]
  meta.data$Barcode <- gsub("-1", "", meta.data$Barcode)
  meta.data$Barcode <- paste0(meta.data$Barcode, "_", row[["Pool"]])
  rownames(meta.data) <- meta.data$Barcode

  # change to a barcode unique across lanes
  colnames(counts) <- meta.data$Barcode

  seurat <- CreateSeuratObject(counts, min.cells = 0, min.features = 0, meta.data = meta.data)
  print("  Done")
  return(seurat)
}
seurat_objects <- apply(pools, 1, load_seurat_object)

if (length(seurat_objects) == 1) {
  seurat <- seurat_objects[[1]]
} else {
  seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])
}

echo("DONE....................................................................",
     "yellow")

#   ____________________________________________________________________________
#   Export data

echo("Saving query data.......................................................",
     "yellow")

saveRDS(seurat, paste0(opt$path, opt$out, ".RDS"))

echo("DONE....................................................................",
     "yellow")