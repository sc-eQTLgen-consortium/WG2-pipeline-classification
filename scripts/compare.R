#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Plot comparison results
# author: Jose Alquicira Hernandez, Lieke Michelsen, Martijn Vochteloo
# date: 2021-11-17
# description: Plots a heatmap based on a contingency table from two columns from
# the Seurat metadata object slot

#   ____________________________________________________________________________
#   Import libraries                                                        ####

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(viridis)))
suppressMessages(suppressWarnings(library(tidyverse)))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option("--metadata1",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--metadata2",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--xaxis",
              type = "character",
              default = "predicted.celltype.l2",
              help = crayon::green("Column in metadata"),
              metavar = "character"),
  make_option("--yaxis",
              type = "character",
              default = "scpred_prediction",
              help = crayon::green("Column in metadata"),
              metavar = "character"),
  make_option("--sort",
              type = "logical",
              default = TRUE,
              help = crayon::green("Sort labels in both axes to match cell type hierarchy?"),
              metavar = "character"),
  make_option("--out",
              type = "character",
              default = "comp",
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

if (is.null(opt$metadata1)){
  print_help(opt_parser)
  stop(crayon::red("Metadata file 1 is missing"), call. = FALSE)
}
if (is.null(opt$metadata2)){
  print_help(opt_parser)
  stop(crayon::red("Metadata file 2 is missing"), call. = FALSE)
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


createheatmap <- function(x, order = TRUE){
  
  xaxis <- x[, 1]
  yaxis <- x[, 2]
  
  cont_table <- table(xaxis, yaxis)
  prop_table <- cont_table / rowSums(cont_table)
  cont_table <-  data.frame(cont_table)
  prop_table <-  data.frame(prop_table)
  cont_table$Prop <- prop_table$Freq
  
  axes_names <- names(x)
  grey_col <-  grey.colors(5000, start = 1, end = 0) # Color of the text in the heatmap
  grey_col[1:2500] <- grey_col[1]
  grey_col[2501:5000] <- grey_col[5000]
  
  
  if(order){
    label_order <- c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "Treg", "CD4 Proliferating", 
                     "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "gdT", "MAIT", "ILC", 
                     "dnT", "NK", "NK_CD56bright", "NK Proliferating", "B naive", "B intermediate", 
                     "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC", 
                     "ASDC", "HSPC", "Platelet", "Eryth", "Doublet")
    
    cont_table[,1] <- factor(cont_table[,1], levels = label_order)
    cont_table[,2] <- factor(cont_table[,2], levels = label_order)
  }
  
  
  # Colormap
  hmcol <- viridis(500)
  grey_col[1] <- hmcol[1] #Plot 0 as invisible
  
  p <- ggplot(cont_table, aes(x = xaxis, y = yaxis, fill = Prop)) +
    geom_tile(colour = "white", size = 0.25) +
    geom_text(aes(label = Freq, colour = Prop), size = 4) +
    scale_colour_gradientn(colours = grey_col, limits = c(0, 1.1), guide = "none") +
    scale_fill_gradientn(colours = hmcol, 
                         limits = c(0, 1), 
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         labels = c(0, 0.25, 0.5, 0.75, 1),
                         guide = "none") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.background = element_blank(),
          panel.border = element_blank(),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    xlab(axes_names[1]) +
    ylab(axes_names[2])
  
  p2 <- ggplot(cont_table, aes(x = xaxis, y = yaxis, fill = Prop)) +
    geom_tile(colour = "white", size = 0.25) +
    geom_text(aes(label = round(Prop, 2), colour = Prop), size = 4) +
    scale_colour_gradientn(colours = grey_col, limits = c(0, 1.1), guide = "none") +
    scale_fill_gradientn(colours = hmcol, 
                         limits = c(0, 1), 
                         breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         labels = c(0, 0.25, 0.5, 0.75, 1),
                         guide = "none") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.background = element_blank(),
          panel.border = element_blank(),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    xlab(axes_names[1]) +
    ylab(axes_names[2])
  
  list(count = p, prop = p2, cont_table = cont_table)
  
}

#   ____________________________________________________________________________
#   Input information                                                       ####

echo("Input information.......................................................",
     "yellow")

echo(paste0(crayon::bold("Metadata file 1:\n"), opt$metadata1), "yellow")
echo(paste0(crayon::bold("Metadata file 2:\n"), opt$metadata2), "yellow")
echo(paste0(crayon::bold("X axis:\n"), opt$xaxis), "yellow")
echo(paste0(crayon::bold("Y axis:\n"), opt$yaxis), "yellow")
echo(paste0(crayon::bold("Sorted labels:\n"), opt$sort), "yellow")
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory for results:\n"), opt$path), "yellow")


echo("DONE....................................................................",
     "yellow")


#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading data......................................................",
     "blue")

metadata1 <- read_delim(opt$metadata1, delim="\t")
metadata2 <- read_delim(opt$metadata2, delim="\t")

echo("DONE....................................................................",
     "blue")

#   ____________________________________________________________________________
#   Merge data                                                              ####

echo("Merging objects.........................................................",
     "blue")

metadata <- merge(metadata1, metadata2, by = 'barcode')

echo("DONE....................................................................",
     "blue")


#   ____________________________________________________________________________
#   Validate columns                                                        ####


# Verify that columns exist in metadata
if(!opt$xaxis %in% colnames(metadata))
  stop("Column '", opt$xaxis, "' is not present in metadata")

if(!opt$yaxis %in% colnames(metadata))
  stop("Column '", opt$yaxis, "' is not present in metadata")

# Verify data type for variables of interest

if(!(is.character(metadata[, opt$xaxis]) | is.factor(metadata[, opt$xaxis])))
  stop("Metadata column '",  opt$xaxis, "' has to be of type character or factor")
if(!(is.character(metadata[, opt$yaxis]) | is.factor(metadata[, opt$yaxis])))
  stop("Metadata column '",  opt$yaxis, "' has to be of type character or factor")

#   ____________________________________________________________________________
#   Plot data                                                               ####

md <- metadata[, c(opt$xaxis, opt$yaxis)]
res <- createheatmap(md, order = opt$sort)

#   ____________________________________________________________________________
#   Export data                                                             ####

echo("Saving results..........................................................",
     "yellow")

ggsave(paste0(opt$path, opt$out, "_heatmap_counts.pdf"), plot = res$count, width = 12, height = 9, dpi = "print")
ggsave(paste0(opt$path, opt$out, "_heatmap_prop.pdf"), plot = res$prop, width = 12, height = 9, dpi = "print")
write_delim(res$cont_table, gzfile(paste0(opt$path, opt$out, "_contingency_table.tsv.gz")), "\t")


echo("DONE....................................................................",
     "yellow")
