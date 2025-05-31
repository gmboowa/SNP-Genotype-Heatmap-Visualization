#!/usr/bin/env Rscript

# categorical_heatmap.R
# This script generates a categorical SNP genotype heatmap (no intermediate values)
# and saves it as PDF and PNG, with a footnote.

# (0) Install and load required packages if not already installed:
required <- c("ComplexHeatmap", "circlize", "grid")
installed <- rownames(installed.packages())
for (pkg in required) {
    if (!pkg %in% installed) {
        install.packages(pkg, repos = "https://cran.rstudio.com/")
    }
}
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})

# (1) Define the genotype matrix (only 0, 1, 2 or NA — no 0.5 or 1.5):
geno_values <- matrix(
    c(0,   NA, 2, 1, 0, 0, 0, 2, 1, 2,
      0,   NA, 2, 0, 0, 0, 0, 1, 1, 0,
      1,   0,  0, 0, 1, 2, 0, 1, 1, 0,
      1,   0,  0, NA, NA, 2, 0, 0, 1, 1,
      0,   1,  0, NA, 0,  1, 0, 1, 1, 0,
      NA,  2, NA, 2,  1, NA, 0, 0, 0, 0),
    nrow = 6, byrow = TRUE
)
rownames(geno_values) <- paste0("Sample_", 1:6)
colnames(geno_values) <- paste0("SNP_", 1:10)

# (2) Define discrete color mapping:
#     0 -> blue (#1f77b4)
#     1 -> orange (#ff7f0e)
#     2 -> red (#d62728)
#     NA -> lightgray
discrete_cols <- c("0" = "#1f77b4", "1" = "#ff7f0e", "2" = "#d62728")

# We wrap this in a function so ComplexHeatmap will treat values as exact categories
get_col <- function(x) {
    # x is the numeric value (0,1,2) or NA; return named color
    if (is.na(x)) {
        return("lightgray")
    }
    # Convert to character to match names(discrete_cols)
    return(discrete_cols[as.character(x)])
}

# (3) Create a matrix of colors of the same shape as geno_values:
color_matrix <- matrix(
    nrow = nrow(geno_values), 
    ncol = ncol(geno_values),
    dimnames = dimnames(geno_values)
)
for (i in seq_len(nrow(geno_values))) {
    for (j in seq_len(ncol(geno_values))) {
        color_matrix[i, j] <- get_col(geno_values[i, j])
    }
}

# (4) Build a Heatmap object using the pre‐computed color matrix:
ht <- Heatmap(
    geno_values,
    name = "Genotype",
    col = discrete_cols,
    na_col = "lightgray",
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        val <- geno_values[i, j]
        if (!is.na(val)) {
            txt_col <- ifelse(val == 2, "white", "black")
            grid.text(
                as.character(val), 
                x = x, 
                y = y, 
                gp = gpar(fontsize = 10, col = txt_col)
            )
        }
    }
)

# (5) Draw and save as PDF with footnote:
pdf("categorical_heatmap_with_footnote.pdf", width = 10, height = 6)
draw(ht, heatmap_legend_side = "right")
pushViewport(viewport(y = unit(0, "npc"), just = "bottom"))
grid.text(
    "Heatmap of SNP Genotypes (0=Ref, 1=Het, 2=Alt, NA=Gray)",
    y = unit(0.02, "npc"),
    x = unit(0.5, "npc"),
    gp = gpar(fontsize = 10, fontface = "italic")
)
popViewport()
dev.off()

# (6) Draw and save as PNG with footnote:
png("categorical_heatmap_with_footnote.png", width = 1000, height = 600)
draw(ht, heatmap_legend_side = "right")
pushViewport(viewport(y = unit(0, "npc"), just = "bottom"))
grid.text(
    "Heatmap of SNP Genotypes (0=Ref, 1=Het, 2=Alt, NA=Gray)",
    y = unit(0.02, "npc"),
    x = unit(0.5, "npc"),
    gp = gpar(fontsize = 10, fontface = "italic")
)
popViewport()
dev.off()
