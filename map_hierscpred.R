#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Cell type classification
# author: Jose Alquicira Hernandez, Lieke Michelsen
# date: 2021-11-17
# description: Classifies cells from scRNA-seq data following hierarchical scPred
# classification approach.

#   ____________________________________________________________________________
#   Import libraries                                                        ####

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("HierscPred"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("progressr"))

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
              help = crayon::yellow("Batch column. If provided, each group in from the batch columns is mapped to reference independently"), 
              metavar = "character"),
  make_option("--thr", 
              type = "numeric", 
              default = 0, 
              help = crayon::yellow("Threshold for rejection. By default no rejection is implemented"), 
              metavar = "character"),
  make_option("--iter", 
              type = "integer", 
              default = 10, 
              help = crayon::yellow("Maximum number or Harmony iterations"), 
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
  make_option("--out", 
              type = "character", 
              default = "hier_scpred", 
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
echo(paste0(crayon::bold("Batch variable:\n"), opt$batch), "yellow")
echo(paste0(crayon::bold("Probability threshold:\n"), opt$thr), "yellow")
echo(paste0(crayon::bold("Harmony iterations:\n"), opt$iter), "yellow")
echo(paste0(crayon::bold("Parallelization plan:\n"), opt$plan), "yellow")
echo(paste0(crayon::bold("Number of workers: \n"), opt$workers), "yellow")
echo(paste0(crayon::bold("maxSize future global: \n"), opt$mem), "yellow")
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


if(opt$plan != "sequential"){
  options(future.globals.maxSize = opt$mem * 1024^3)
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

echo("Loading scPred reference................................................",
     "yellow")

reference <- readRDS("/hier_scpred.RDS")

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

echo("Classifiying cells .....................................................",
     "green")

classify <- function(xs){
  p <- progressor(along = xs)
  future_mapply(function(x, i, n){
    x <- predictTree(reference, newData = x, threshold = opt$thr, max.iter.harmony = opt$iter)
    p(message = sprintf("| Batch %d/%d", i, n))
    x
  }, xs, seq_along(xs), MoreArgs = list(n = length(xs)), SIMPLIFY = FALSE, future.seed = TRUE)
}

batches <- classify(batches)

echo("DONE....................................................................",
     "green")


#   ____________________________________________________________________________
#   Gather cell type classification and store in main object                ####

echo("Gathering results........................................................",
     "green")

scpred_prediction <- lapply(batches,
                            function(x) x[[]][,
                                              "scpred_prediction", drop = FALSE
                            ])



scpred_prediction <- do.call(rbind, scpred_prediction)
data <- AddMetaData(data, scpred_prediction)

echo("DONE....................................................................",
     "green")

#   ____________________________________________________________________________
#   Export data                                                             ####

echo("Saving query data.......................................................",
     "yellow")

saveRDS(data, file.path(opt$path, paste0(opt$out,".RDS")))

echo("DONE....................................................................",
     "yellow")
