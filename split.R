#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Split Seurat object
# author: Jose Alquicira Hernandez, Lieke Michelsen
# date: 2021-11-17
# description: Split data by batches based on a column in the Seurat object 
# metadata

#   ____________________________________________________________________________
#   Import libraries                                                        ####

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("progressr"))
suppressPackageStartupMessages(library("optparse"))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option("--file",
              type = "character",
              default = NULL,
              help = crayon::green("RDS object file name"),
              metavar = "character"),
  make_option("--batch",
              type = "character",
              default = NULL,
              help = crayon::yellow("Batch column to split the data"),
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

if (is.null(opt$file)){
  print_help(opt_parser)
  stop(crayon::red("Query file name is missing"), call. = FALSE)
}

if (is.null(opt$batch)){
  print_help(opt_parser)
  stop(crayon::red("Batch is missing"), call. = FALSE)
}

#   ____________________________________________________________________________
#   Validate output                                                         ####

if(!dir.exists(opt$path)) stop("Output path cannot be found: ", opt$path)

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

echo(paste0(crayon::bold("Input file: "), opt$file), "yellow")
echo(paste0(crayon::bold("Batch variable: "), opt$batch), "yellow")
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory: "), opt$path), "yellow")


echo("DONE....................................................................\n",
     "yellow")


#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading query data......................................................",
     "blue")

data <- readRDS(opt$file)
if(!inherits(data, "Seurat")) stop("Input query data is not a Seurat object")

echo("DONE....................................................................\n",
     "blue")


#   ____________________________________________________________________________
#   Validate batch                                                          ####

# Verify that columns exist in metadata
if(!opt$batch %in% names(data[[]]))
  stop("Batch column '", xaxis, "' is not present in metadata")

#   ____________________________________________________________________________
#   Plot data                                                               ####

echo("Splitting data..........................................................", 
     "green")

batches <- SplitObject(data, split.by = opt$batch)

echo("DONE....................................................................\n", 
     "green")


#   ____________________________________________________________________________
#   Export data                                                             ####

echo("Saving results..........................................................",
     "yellow")

handlers(global = TRUE)
handlers("progress")

save_batches<- function(xs){
  p <- progressor(along = xs)
  mapply(function(x, name, i, n){
    saveRDS(x, file = file.path(opt$path, paste0(opt$out, "-", name, ".RDS")))
    p(message = sprintf("| Batch %d/%d", i, n))
  }, xs, names(xs), seq_along(xs), MoreArgs = list(n = length(xs)), SIMPLIFY = FALSE)
}

dump <- save_batches(batches)

echo("DONE....................................................................\n",
     "yellow")
