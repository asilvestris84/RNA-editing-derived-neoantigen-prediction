#!/usr/bin/env Rscript 
library(ggplot2)
library(dplyr)
library(tidyr)

# Carica i dati
data <- read.delim("/lustrehome/pietrolucamazzacuva/share/Pancreas_Alessandro_Miriam/common_strong_binders filtered_classII_with_GTEx_data.tsv", stringsAsFactors = FALSE)
data <- data[, -1]
data_unica <- data %>%
  distinct(Mutation_Position_Hg38, .keep_all = TRUE)
df_long <- data_unica %>%
  separate_rows(Healthy_Tissues, Median_Editing_Frequency_., SRR_Samples_., sep = ",") %>%
  mutate(
    Median_Editing_Frequency_. = as.numeric(Median_Editing_Frequency_.),
    SRR_Samples_. = as.numeric(SRR_Samples_.)
  ) %>%
  filter(Healthy_Tissues != "", !is.na(Median_Editing_Frequency_.))
df_long$Healthy_Tissues <- factor(df_long$Healthy_Tissues, levels = sort(unique(df_long$Healthy_Tissues)))
df_long[[1]] <- paste0(df_long[[1]], "_", df_long[[4]])

# Crea il plot
png("bubble_plot_high_res_MHC_II.png", w=12, h=10, pointsize=33, res = 300, units = 'in')
ggplot(df_long, aes(
  x = Gene,
  y = Healthy_Tissues,
  )) +
  geom_point(alpha = 0.7,aes(color=Median_Editing_Frequency_., size=SRR_Samples_.)) +
  scale_size(range = c(1, 7)) +
  theme_minimal() +
  ylab("GTEx tissues") +
  xlab("RNA editing sites") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), plot.title = element_text(hjust = 0, face="bold", size=12),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position="right") +
  scale_color_gradient(low = "lightblue", high = "blue")+
  guides(color =guide_colourbar(title="Median editing frequency", order=1), size = guide_legend(title = "Percentage of healthy samples", override.aes = list(alpha = 0.3))) +
  labs(
    fill = "Median editing frequency",
    size = "Percentage of healthy samples",
    subtitle = "from recoding RNA editing") +
  ggtitle("Class I neoantigens") 

dev.off()

#ggsave(
#  filename = "bubble_plot_high_res.png",   # nome file
#  plot = bubble_plot,                      # oggetto grafico
#  width = 12, height = 8,                  # dimensioni in pollici
#  dpi = 300                                # risoluzione
#)
  

# Visualizza
bubble_plot
