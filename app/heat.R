library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
library(readr, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)

args <- commandArgs(trailingOnly = T)
path <- args[1:2]
dims <- as.numeric(args[3:4])
unit <- args[5]
dpi <- as.integer(args[6])
exts <- args[7:length(args)]

p <- (
  read_tsv(path[1], col_types = cols(psim = "d", .default = "c")) %>%
    mutate(
      call = factor(call, levels = c("TP", "FN", "FP", "TN", "XX")),
      com = factor(com, levels = rev(c("F", "P", "R", "F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3")))
    ) %>%
    ggplot(aes(acc, com, fill = psim)) +
    geom_tile() +
    facet_grid(assay ~ call, space = "free", scales = "free") +
    scale_x_discrete(labels = NULL) +
    scale_fill_gradientn(
      colors = c("#440154", "#3B528B", "#21918C", "#5ec962", "#fde725"),
      labels = c("0", "75", " |\n80", "|\n99", "   100"),
      values = c(0.0, 0.75, 0.8, 0.99, 1.0),
      breaks = c(0.0, 0.75, 0.8, 0.99, 1.0),
      limits = c(0, 1)
    ) +
    labs(x = "accession", y = "component") +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      panel.background = element_rect(fill = "#F3F3F3"),
      panel.border = element_rect(fill = NA, linetype = "dotted"),
      strip.background = element_rect(color = "black"),
      strip.text.x = element_text(hjust = 0.5, vjust = 0),
      strip.text.y = element_text(hjust = 0, vjust = 0.5, angle = 0)
    )

)

cat(
  paste(
    lapply(
      exts, 
      function(ext) {
        ggsave(path <- paste(path[2], ext, sep="."), plot = p, width = dims[1], height = dims[2], units = unit, dpi=dpi)
        path
      }
    ),
    collapse="\n"
  )
)
