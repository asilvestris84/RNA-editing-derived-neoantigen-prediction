#!/usr/bin/env Rscript

## -------------------------------------------------------
## Calcolo PriorityScore gerarchico per neoantigeni editati
## RNA_norm × Freq_norm + 0.25 × ΔRank_norm
## -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0)
  stop("Uso: Rscript priority_score_neoAg.R file1.txt [file2.txt ...]\n")

## Normalizzatori
normalize_pos <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r) == 0)
    return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
}

normalize_inv <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r) == 0)
    return(rep(0.5, length(x)))
  1 - ((x - r[1]) / diff(r))
}

process_file <- function(infile) {
  message("\n>>> Processo file: ", infile)

  df <- read.table(infile, header = TRUE, sep = "\t",
                   quote = "", check.names = FALSE,
                   stringsAsFactors = FALSE)

  needed <- c("Mut_Rank", "Wt_Rank", "Frequency", "Total_RNA")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0)
    stop("Nel file mancano le colonne: ", paste(missing, collapse=", "))

  ## DeltaRank = quanto migliora l'affinità il mutante rispetto al WT
  df$DeltaRank <- df$Mut_Rank / df$Wt_Rank

  ## Normalizzazioni
  df$RNA_norm       <- normalize_pos(df$Total_RNA)
  df$Freq_norm      <- normalize_pos(df$Frequency)
  df$DeltaRank_norm <- normalize_inv(df$DeltaRank)

  ## PriorityScore gerarchico
  df$PriorityScore <-
      (df$RNA_norm * df$Freq_norm) +
      0.25 * df$DeltaRank_norm

  ## Output
  outfile <- sub("\\.txt$|\\.tsv$|\\.tab$", "", infile)
  outfile <- paste0(outfile, "_priorityScore.txt")

  write.table(df, file = outfile, sep = "\t",
              quote = FALSE, row.names = FALSE)

  message("   -> OUTPUT salvato in: ", outfile)
}

## Loop files
for (f in args)
  process_file(f)

message("\n✔️ Finito. Score generati con α = 0.25\n")

q()

library(ggplot2)

df <- df[order(-df$PriorityScore), ]
top20 <- head(df, 20)

png("neoAg_priority_top20.png", width=1200, height=900, res=150)
ggplot(top20, aes(x=reorder(Gene, PriorityScore), y=PriorityScore)) +
  geom_col(fill="#2255AA") +
  coord_flip() +
  theme_bw(base_size=14) +
  labs(title="Top prioritized RNA-editing-derived neoantigens",
       x="Gene / Peptide", y="Priority score")
dev.off()













