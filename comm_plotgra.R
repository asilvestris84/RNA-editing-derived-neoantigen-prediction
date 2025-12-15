#!/bin/env Rscript

# Carica le librerie necessarie
library(dplyr)
library(reshape2)
library(ggplot2)

options(bitmapType='cairo')
# Definisci la cartella contenente i file di output di Neoantimon
data_folder <- "/lustre/biomed/asilvestris/dati_Cives/b_mandriani_pancreas/neoAg_editing_def/result.ID.SNV1/"  # Sostituisci con il percorso corretto

# Elenca tutti i file .ALL.txt nella cartella
all_files <- list.files(path = data_folder, pattern = "^filtered_.*\\.ALL.txt$", full.names = TRUE)

# Inizializza una lista per memorizzare i neoantigeni per campione
neoantigen_list <- list()

# Funzione per filtrare i "strong binders" in base solo a Mut_Rank
filter_strong_binders <- function(data) {
  data <- data %>%
    mutate(Mut_Rank = as.numeric(Mut_Rank)) %>%
    filter(Mut_Rank <= 1)  # Soglia per strong binders basata solo su Mut_Rank
  return(data)
}

# Leggi ciascun file e filtra i "strong binders"
for (file in all_files) {
  sample_name <- gsub(".ALL.txt$", "", basename(file))  # Estrai il nome del campione
  data <- read.delim(file, sep = "\t", header = TRUE)
  
  # Filtra i "strong binders"
  strong_binders <- filter_strong_binders(data)
  
  # Memorizza i peptidi mutanti valutati per il campione
  neoantigen_list[[sample_name]] <- unique(strong_binders$Evaluated_Mutant_Peptide)
}

# Crea una lista di tutti i neoantigeni unici
all_neoantigens <- unique(unlist(neoantigen_list))

# Crea una matrice binaria di presenza/assenza
presence_matrix <- sapply(neoantigen_list, function(x) all_neoantigens %in% x)
rownames(presence_matrix) <- all_neoantigens

# Calcola la matrice di co-occorrenza
cooccurrence_matrix <- t(presence_matrix) %*% presence_matrix

# Imposta la diagonale a NA per escluderla dalla visualizzazione
diag(cooccurrence_matrix) <- NA
#cooccurrence_matrix[lower.tri(cooccurrence_matrix)] <- NA #crea una maschera booleana mantenendo soltanto la parte superiore della matrice per i valori che si ripetono
cooccurrence_matrix[upper.tri(cooccurrence_matrix)] <- NA

# Dopo aver creato cooccurrence_matrix
rownames(cooccurrence_matrix) <- gsub(".*_(ICGC_\\d+)\\..*", "\\1", rownames(cooccurrence_matrix))
colnames(cooccurrence_matrix) <- gsub(".*_(ICGC_\\d+)\\..*", "\\1", colnames(cooccurrence_matrix))

# Converte la matrice in un data frame lungo per ggplot2
cooccurrence_df <- melt(cooccurrence_matrix, na.rm = TRUE)

head(cooccurrence_df)

# Crea l'emi-matrice di correlazione

png("MHC_I_cooccurrences.png", w=24, h=20, pointsize=33, res = 300, units = 'in')

ggplot(cooccurrence_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 22),  # Aumenta la dimensione delle etichette
        axis.text.y = element_text(size = 22),                         # Aumenta la dimensione delle etichette
        panel.grid = element_blank(),
        legend.title = element_text(size=20),
        plot.subtitle=element_text(size=30),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25)) + #dimensione assi x e y
  labs(fill = "Common\nNeoantigens", subtitle = "Derived from recoding RNA editing") +
  xlab(" ") +
  ylab("PanNETs") +
  coord_fixed() +
  ggtitle("Class I Neoantigens") +
  theme(plot.title = element_text(hjust = 0.5, size = 35)) +  # Aumenta la dimensione del titolo
  geom_text(aes(label =value), size = 10, color = "black",data=cooccurrence_df[cooccurrence_df$value != 0,])

dev.off()

###############################################################################################################

# Inizializza un dataframe per accumulare le informazioni
all_strong_binders <- data.frame()

for (file in all_files) {
  sample_name <- gsub(".ALL.txt$", "", basename(file))
  data <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Filtra i strong binders
  strong_binders <- filter_strong_binders(data)
  
  if(nrow(strong_binders) > 0){
    strong_binders$Sample <- sample_name
    all_strong_binders <- rbind(all_strong_binders, strong_binders)
  }
}

# Raggruppa per Evaluated_Mutant_Peptide e HLA_Class per ottenere una lista non ridondante.
# - Escludiamo Mut_Rank.
# - Nel nome del gene, rimuoviamo il prefisso numerico e l'underscore (ad es. "123_GENE" diventa "GENE").
# - Inseriamo NM_ID e Change (presi dal primo record del gruppo).
nonredundant_strong_binders <- all_strong_binders %>%
  group_by(Evaluated_Mutant_Peptide) %>%
  summarise(
    Gene = sub("^[0-9]+_", "", first(Gene)),
    NM_ID = first(NM_ID),
    Change = first(Change),
    Mutation_Position = first(Mutation_Position),
    Peptide_Sequence = first(Evaluated_Mutant_Peptide),
    Samples = paste(unique(Sample), collapse = ","),
    N_Samples = n_distinct(Sample)
  ) %>%
  ungroup()

# Separa in due tabelle: una per HLA Classe I e una per HLA Classe II
#strong_binders_classI <- nonredundant_strong_binders %>% filter(HLA_Class == "I")
#strong_binders_classII <- nonredundant_strong_binders %>% filter(HLA_Class == "II")

# Salva le tabelle in file TSV
#write.table(strong_binders_classI, file = "strong_binders_classI.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(nonredundant_strong_binders, file = "strong_binders_classI.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#cat("Numero di strong binders unici per HLA Classe I:", nrow(strong_binders_classI), "\n")
#cat("Numero di strong binders unici per HLA Classe II:", nrow(strong_binders_classII), "\n")





























