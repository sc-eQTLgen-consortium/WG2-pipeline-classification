#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Merge multiple Seurat objects
# author: Jose Alquicira Hernandez, Lieke Michelsen
# date: 2021-11-17
# description: Aggregates input Seurat objects into a single object

#   ____________________________________________________________________________
#   Import libraries                                                        ####

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("optparse"))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option(c("--file"),
              type = "character",
              default = NULL,
              help = crayon::green("Directory filename containing Seurat RDS objects"),
              metavar = "character"),
  make_option(c("--out"), 
              type = "character", 
              default = "reduced_data", 
              help = crayon::green("Output file name [default= %default]"), 
              metavar = "character"),
  make_option(c("--path"), 
              type = "character", 
              default = ".", 
              help = crayon::green("Output path to store results [default= %default]"), 
              metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop(crayon::red("Directory name is missing"), call. = FALSE)
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
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory: "), opt$path), "yellow")


echo("DONE....................................................................",
     "yellow")

#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading query data......................................................",
     "blue")

list_dirs <- list.files(opt$file, full.names = TRUE)
list_dirs <- grep(".RDS$", list_dirs, value = TRUE)

data <- lapply(list_dirs, function(x){
  message("Reading ", x)
  x <- readRDS(x)
  if(!inherits(x, "Seurat")) stop("Input query data is not a Seurat object")
  x
})

echo("DONE....................................................................",
     "blue")

#   ____________________________________________________________________________
#   Merge data                                                              ####

echo("Merging objects.........................................................",
     "blue")

reductions <- lapply(data, Reductions)
reductions <- unique(unlist(reductions))

data <- merge(data[[1]], data[-1], merge.dr = reductions)

echo("DONE....................................................................",
     "blue")
#   ____________________________________________________________________________
#   Export data                                                             ####

echo("Saving results..........................................................",
     "yellow")

saveRDS(data, file = file.path(opt$path, paste0(opt$out, ".RDS")))

echo("DONE....................................................................",
     "yellow")
