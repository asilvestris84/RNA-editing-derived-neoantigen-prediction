# ───────────────────────────────────────────────
# Librerie
# ───────────────────────────────────────────────
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(grid)

options(bitmapType = "cairo")

# ───────────────────────────────────────────────
# Lettura file e selezione colonne SRR/ICGC
# ───────────────────────────────────────────────
df <- read.table("recoding_heatmap_matrix.txt", header = TRUE, sep = "\t", check.names = FALSE)
sample_cols <- grep("^(SRR|ICGC)", colnames(df), value = TRUE)
mat <- as.matrix(df[, sample_cols]) * 100  # livelli % di editing

# ───────────────────────────────────────────────
# Rinomina righe nel formato Gene_p.XYZ
# ───────────────────────────────────────────────
# Funzione per estrarre la prima variante "p." trovata
extract_main_p_change <- function(aa_change) {
    p_changes <- unlist(str_extract_all(aa_change, "p\\.[A-Za-z0-9]+"))
    if (length(p_changes) > 0) {
        return(names(sort(table(p_changes), decreasing = TRUE))[1])
    } else {
        return("p.?")
    }
}

# Applica su colonna AAChange
gene_names <- df$Gene
aa_changes <- sapply(df$AAChange, extract_main_p_change)
rownames(mat) <- make.unique(paste(gene_names, aa_changes, sep = "_"))
head(mat)
condition <- ifelse(grepl("^SRR", sample_cols), "Ctrl", "pNET")
condition <- factor(condition, levels = c("Ctrl", "pNET"))
group_colors <- c("Ctrl" = "#377EB8", "pNET" = "#E41A1C")
ha_col <- HeatmapAnnotation(
  Condition = condition,
  col = list(Condition = group_colors),
  annotation_name_side = "left"
)
head(condition)

ctrl_mat <- mat[, grepl("^SRR", colnames(mat))]
tumor_mat <- mat[, grepl("^ICGC", colnames(mat))]
# Trasforma NA in 0 o lascia NA (dipende dalla scelta)
ctrl_mat[is.na(ctrl_mat)] <- NA
tumor_mat[is.na(tumor_mat)] <- NA
# Crea liste di data.frame per boxplot
ctrl_data_list <- lapply(1:nrow(ctrl_mat), function(i) data.frame(value = as.numeric(ctrl_mat[i, ])))
tumor_data_list <- lapply(1:nrow(tumor_mat), function(i) data.frame(value = as.numeric(tumor_mat[i, ])))

ctrl_data_list <- lapply(1:nrow(ctrl_mat), function(i) as.numeric(ctrl_mat[i, ]))
tumor_data_list <- lapply(1:nrow(tumor_mat), function(i) as.numeric(tumor_mat[i, ]))

max_val <- max(mat, na.rm = TRUE)
ylim_range <- c(0, max_val)
row_anno <- rowAnnotation(
  Ctrl = anno_boxplot(
    ctrl_data_list,
    which = "row",
    axis_param = list(direction = "reverse"),
    gp = gpar(fill = "#377EB8"),
    width = unit(1.2, "cm"),
    border = TRUE,
    ylim = ylim_range
  ),
  pNET = anno_boxplot(
    tumor_data_list,
    which = "row",
    axis_param = list(direction = "reverse"),
    gp = gpar(fill = "#E41A1C"),
    width = unit(1.2, "cm"),
    border = TRUE,
    ylim = ylim_range
  )
)
png("heatmap_editing_SRR_ICGC_gene_AApchange_boxplot.png", width = 3000, height = 3200, res = 300)
Heatmap(
  mat,
  name = "Editing Level (%)",
  col = colorRamp2(c(0, 50, 100), c("white", "skyblue", "darkblue")),
  na_col = "grey90",
  top_annotation = ha_col,
  left_annotation = row_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = "Editing Level (%)",
    title_position = "leftcenter-rot"
  )
)
dev.off()

