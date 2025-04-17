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
  read_tsv(path[1], col_types = cols(pos = "i", count = "i", .default = "c")) %>%
    mutate(
      call = factor(call, levels = c("TP", "FN", "FP", "TN", "XX")),
      com = factor(com, levels = c("F", "P", "R", "F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3")),
      mut = factor(mut, levels = c("TRS", "TRV", "INS", "DEL", "DIS", "GAP", "UNK"))
    ) %>%
    ggplot(aes(x=pos, y=count)) +
    geom_bar(aes(fill = mut), stat = "identity") +
    facet_grid(assay ~ com) +
    guides(fill = guide_legend(nrow=1)) +
    theme(
        panel.border = element_rect(fill = NA, linetype = "dotted"),
        strip.background = element_rect(color = "black"),
        strip.text.x = element_text(hjust = 0.5, vjust = 0),
        strip.text.y = element_text(hjust = 0, vjust = 0.5, angle = 0),
        legend.position = "bottom"
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
