#!/usr/bin/env Rscript

# categorical_heatmap.R
# This script generates a categorical SNP genotype heatmap and saves it as PDF and PNG.

# (0) Install and load required packages if not already installed:
required <- c("ComplexHeatmap", "circlize", "grid")
installed <- rownames(installed.packages())
for(pkg in required) {
    if (!pkg %in% installed) {
        install.packages(pkg, repos="https://cran.rstudio.com/")
    }
}
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})

# (1) Define the genotype matrix:
geno_values <- matrix(
    c(0,   NA, 2, 1, 0, 0, 0, 2, 1, 2,
      0,   NA, 2, 0, 0, 0, 0, 1, 1, 0,
      1,   0,  0, 0, 1, 2, 0, 1, 1, 0,
      1,   0,  0, NA, NA,2, 0, 0, 1, 1,
      0,   1,  0, NA, 0,  1, 0, 1, 1, 0,
      NA,  2, NA,  2, 1,  NA,0, 0, 0, 0),
    nrow = 6, byrow = TRUE
)
rownames(geno_values) <- paste0("Sample_", 1:6)
colnames(geno_values) <- paste0("SNP_", 1:10)

# (2) Define color mapping: 0->blue, 1->orange, 2->red, NA->lightgray
col_mapping <- colorRamp2(c(0, 1, 2), c("#1f77b4", "#ff7f0e", "#d62728"))

# (3) Create the heatmap:
ht <- Heatmap(
    geno_values,
    name = "Genotype",
    col = col_mapping,
    na_col = "lightgray",
    rect_gp = gpar(col = "black"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        val <- geno_values[i, j]
        if (!is.na(val)) {
            txt_col <- ifelse(val == 2, "white", "black")
            grid.text(as.character(val), x = x, y = y, gp = gpar(fontsize = 10, col = txt_col))
        }
    }
)

# (4) Draw and save as PDF with footnote:
pdf("categorical_heatmap_with_footnote.pdf", width = 10, height = 6)
draw(ht, heatmap_legend_side = "right")
pushViewport(viewport(y = unit(0, "npc"), just = "bottom"))
grid.text("Heatmap of SNP Genotypes (0=Ref, 1=Het, 2=Alt, NA=Gray)",
          y = unit(0.02, "npc"), x = unit(0.5, "npc"),
          gp = gpar(fontsize = 10, fontface = "italic"))
popViewport()
dev.off()

# (5) Draw and save as PNG with footnote:
png("categorical_heatmap_with_footnote.png", width = 1000, height = 600)
draw(ht, heatmap_legend_side = "right")
pushViewport(viewport(y = unit(0, "npc"), just = "bottom"))
grid.text("Heatmap of SNP Genotypes (0=Ref, 1=Het, 2=Alt, NA=Gray)",
          y = unit(0.02, "npc"), x = unit(0.5, "npc"),
          gp = gpar(fontsize = 10, fontface = "italic"))
popViewport()
dev.off()
