#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

options(bitmapType = "cairo")

## =========================================================
## PARAMETRO DA MODIFICARE
## =========================================================
## File di testo con UNO per riga i full path delle tabelle
## (ognuna con la colonna PriorityScore già calcolata)
file_list <- "ICGC_PS_MHCII_path"  # <--- cambia solo questo

## Nome del PNG in output
output_png <- "neoAg_priority_heatmap.png"

## =========================================================
## LETTURA LISTA FILE
## =========================================================

if (!file.exists(file_list)) {
  stop("File lista non trovato: ", file_list)
}

paths <- readLines(file_list, warn = FALSE)
paths <- paths[nzchar(paths)]

if (!length(paths)) stop("Lista percorsi vuota.")

message("Tabelle da processare: ", length(paths))

## =========================================================
## Funzione per estrarre ID campione (ICGC_XXXX) dal nome file
## =========================================================

get_sample_id_from_path <- function(path) {
  fname <- basename(path)
  m <- regexpr("ICGC_[0-9A-Za-z]+", fname)
  if (m > 0) {
    return(regmatches(fname, m))
  } else {
    # fallback: nome file senza estensione
    return(sub("\\..*$", "", fname))
  }
}

## =========================================================
## LETTURA DI TUTTE LE TABELLE E COSTRUZIONE TABELLA LONG
## =========================================================

all_list <- list()

for (f in paths) {
  message("Leggo: ", f)
  if (!file.exists(f)) {
    warning("File non trovato, salto: ", f)
    next
  }
  
  dt <- fread(f, sep = "\t", header = TRUE, check.names = FALSE)
  
  needed <- c("Gene", "Evaluated_Mutant_Peptide", "PriorityScore")
  miss <- setdiff(needed, colnames(dt))
  if (length(miss)) {
    warning("Nel file ", f, " mancano colonne: ",
            paste(miss, collapse = ", "),
            "  --> file ignorato.")
    next
  }
  
  sample_id <- get_sample_id_from_path(f)
  
  ## ---- CHIAVE UNIVoca: GENE PULITO + PEPTIDE ----
  ## Gene contiene roba tipo "29_CCNI": tengo solo "CCNI"
  dt[, Gene_clean := sub("^[0-9]+_", "", Gene)]
  
  ## Peptide ID = Gene_clean + "_" + sequenza mutante
  dt[, peptide_id := paste(Gene_clean, Evaluated_Mutant_Peptide, sep = "_")]
  dt[, sample_id  := sample_id]
  
  all_list[[length(all_list) + 1]] <-
    dt[, .(peptide_id, sample_id, PriorityScore)]
}

if (!length(all_list)) stop("Nessun file valido letto.")

all_dt <- rbindlist(all_list, use.names = TRUE)

## Se nello stesso campione ho più righe con stesso peptide_id
## (es. isoforme), prendo il PriorityScore massimo
agg_dt <- all_dt[
  ,
  .(PriorityScore = max(PriorityScore, na.rm = TRUE)),
  by = .(peptide_id, sample_id)
]

## =========================================================
## COSTRUZIONE MATRICE PEPTIDE × CAMPIONE (LEFT JOIN globale)
## =========================================================

peptides <- sort(unique(agg_dt$peptide_id))   # 63 attesi
samples  <- sort(unique(agg_dt$sample_id))

mat <- matrix(NA_real_,
              nrow = length(peptides),
              ncol = length(samples),
              dimnames = list(peptides, samples))

for (i in seq_len(nrow(agg_dt))) {
  r <- agg_dt$peptide_id[i]
  c <- agg_dt$sample_id[i]
  s <- agg_dt$PriorityScore[i]
  if (is.finite(s)) {
    mat[r, c] <- s
  }
}

message("Peptidi unici (righe):   ", nrow(mat))
message("Campioni (colonne):      ", ncol(mat))

## =========================================================
## HEATMAP: colore = PriorityScore, NA = grigio
## =========================================================

col_fun <- colorRamp2(c(0, 0.5, 1),
                      c("white", "orange", "red"))

png(output_png, width = 3200, height = 2400, res = 300)

ht <- Heatmap(
  mat,
  name = "PriorityScore",
  col = col_fun,
  na_col = "grey82",           # peptide non presente nel campione
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Priority score",
    title_position = "leftcenter-rot"
  )
)

draw(
  ht,
  padding = unit(c(3, 10, 3, 3), "mm")  # top, right, bottom, left
)

dev.off()

message("\nHeatmap salvata in: ", output_png, "\n")

