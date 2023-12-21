#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Visualise
# author: Martijn Vochteloo
# date: 2023-12-21
# description: Aggregates input Seurat objects into a single object

#   ____________________________________________________________________________
#   Import libraries                                                        ####

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(dplyr)))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option("--metadata",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option(c("--columns"),
              type = "character",
              default = NULL,
              help = crayon::green("The cell type columns"),
              metavar = "character"),
  make_option(c("--palette"),
              type = "character",
              default = NULL,
              help = crayon::green("The palette to use"),
              metavar = "character"),
  make_option(c("--out"),
              type = "character",
              default = "visualise",
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

if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop(crayon::red("Metadata file is missing"), call. = FALSE)
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

echo(paste0(crayon::bold("Metadata file:\n"), opt$metadata), "yellow")
echo(paste0(crayon::bold("Cell type columns: "), opt$columns), "yellow")
echo(paste0(crayon::bold("Palette: "), opt$palette), "yellow")
echo(paste0(crayon::bold("Output base filename: "), opt$out), "yellow")
echo(paste0(crayon::bold("Output directory: "), opt$path), "yellow")


echo("DONE....................................................................",
     "yellow")

#   ____________________________________________________________________________
#   Import query data                                                       ####

echo("Loading data...........................................................",
     "blue")

metadata <- read_delim(opt$metadata, delim="\t")

if(opt$out == "azimuth"){
  ct_columns <- c()
  for (col in unlist(strsplit(opt$columns, ";"))) {
    keyvalue <- unlist(strsplit(col, "="))
    ct_columns <- c(ct_columns, paste0("predicted.", keyvalue[[1]]))
  }
} else {
  ct_columns <- c(opt$columns)
}

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



echo("DONE....................................................................",
     "blue")

#   ____________________________________________________________________________
#   Cell proportions                                                        ####

echo("Plotting cell proportions..............................................",
     "blue")

n_total <- dim(metadata)[1]
n_fun <- function(x) {
  return(data.frame(y = x, label = paste0(x, "\n(", round(x / n_total,2) * 100, "%)")))
}
for (ct_col in ct_columns) {
  occurrences <- as.data.frame(table(metadata[[ct_col]]))
  colnames(occurrences) <- c(ct_col, "N")
  p <- ggplot(data = occurrences,
              aes(
                x = .data[[ct_col]],
                y = .data[["N"]],
                fill = .data[[ct_col]],
              )
  ) +
    geom_bar(stat="identity") +
    ylab("#cells") +
    xlab("cell type") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    stat_summary(fun.data = n_fun, geom = "text", vjust = -0.2)
  if (!is.null(palette)) {
    use_palette <- TRUE
    for (value in unique(occurrences[[ct_col]])) {
      if (!value %in% names(palette)) {
        use_palette <- FALSE
      }
    }
    if (use_palette) {
      p <- p + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
    }
  }
  ggsave(paste0(opt$path, opt$out, "_", ct_col, "_proportion.png"),
         p,
         width = 20,
         height = 10,
  )
}

echo("DONE....................................................................",
     "blue")

#   ____________________________________________________________________________
#   Prediction scores

echo("Plotting prediction scores..............................................",
     "blue")

n_fun <- function(x) {
  return(data.frame(y = 1.02, label = paste0("n = ", length(x))))
}

for (ct_col in ct_columns) {
  if (!paste0(ct_col, ".score") %in% colnames(metadata)) {
    next
  }
  p <- ggplot(data = metadata,
              aes(
                x = .data[[ct_col]],
                y = .data[[paste0(ct_col, ".score")]],
                fill = .data[[ct_col]],
              )
  ) +
    geom_boxplot(alpha = 0.5) +
    geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2)) +
    ylab("prediction score") +
    xlab("cell type") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    stat_summary(fun.data = n_fun, geom = "text")
  if (!is.null(palette)) {
    use_palette <- TRUE
    for (value in unique(metadata[[ct_col]])) {
      if (!value %in% names(palette)) {
        use_palette <- FALSE
      }
    }
    if (use_palette) {
      p <- p + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
    }
  }
  ggsave(paste0(opt$path, opt$out, "_", ct_col, "_score.png"),
         p,
         width = 20,
         height = 10,
  )
}

echo("DONE....................................................................",
     "blue")