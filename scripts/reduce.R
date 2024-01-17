#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Merge multiple Seurat objects
# author: Jose Alquicira Hernandez, Lieke Michelsen, Martijn Vochteloo
# date: 2021-11-17
# description: Aggregates input Seurat objects into a single object

#   ____________________________________________________________________________
#   Import libraries                                                        ####

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(tidyverse)))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option("--poolsheet",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--indir",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--out",
              type = "character", 
              default = NULL,
              help = crayon::green("Output file name"),
              metavar = "character"),
  make_option("--path",
              type = "character", 
              default = ".",
              help = crayon::green("Output path to store results [default= %default]"), 
              metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$poolsheet)){
  print_help(opt_parser)
  stop(crayon::red("Poolsheet file is missing"), call. = FALSE)
}
if (is.null(opt$indir)){
  print_help(opt_parser)
  stop(crayon::red("Input directory name is missing"), call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop(crayon::red("Output file name is missing"), call. = FALSE)
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

echo(paste0(crayon::bold("Poolsheet file:\n"), opt$poolsheet), "yellow")
echo(paste0(crayon::bold("Input directory:\n"), opt$indir), "yellow")
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory: "), opt$path), "yellow")


echo("DONE....................................................................",
     "yellow")

#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading metadata......................................................",
     "blue")

pools <- read_delim(opt$poolsheet, delim = "\t")

load_metadata <- function(row){
  message("Reading ", row[["Pool"]])
  metadata <- read_delim(paste0(opt$indir, opt$out, "_", row[["Pool"]], ".metadata.tsv.gz"), delim="\t")
  metadata <- cbind(Pool = row[["Pool"]], metadata)
  return(metadata)
}
metadata <- apply(pools, 1, load_metadata)

echo("DONE....................................................................",
     "blue")

#   ____________________________________________________________________________
#   Merge data                                                              ####

echo("Merging objects.........................................................",
     "blue")

metadata <- as.data.frame(do.call(rbind,metadata))

echo("DONE....................................................................",
     "blue")
#   ____________________________________________________________________________
#   Export data                                                             ####

echo("Saving results..........................................................",
     "yellow")

write_delim(metadata, gzfile(paste0(opt$path, opt$out, ".metadata.tsv.gz")), "\t")

echo("DONE....................................................................",
     "yellow")
