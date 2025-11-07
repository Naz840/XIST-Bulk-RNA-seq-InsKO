# ====================================================== 
# Beta Cell Pathway Analysis
# Identifying pathways related to beta cell function, metabolism, and health
# ======================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Set working directory
setwd("C:/Users/mqadir/Box/Fahd shared to MBRLab/Figures/raw/BetaCell_function/ORA")

# Load data
down_go <- read.csv("DOWN_GO_enrichment_full.csv", stringsAsFactors = FALSE)
up_go <- read.csv("UP_GO_enrichment_full.csv", stringsAsFactors = FALSE)

# Add regulation status
down_go$regulation <- "DOWN"
up_go$regulation <- "UP"

# Combine datasets
all_go <- rbind(down_go, up_go)

# Define beta cell-related pathway keywords
beta_cell_keywords <- list(
  `Beta Cell Function` = c("insulin", "beta.cell", "β.cell", "beta cell", "islet", "pancrea", 
                           "glucagon", "somatostatin", "endocrine", "incretin", "glp1"),
  
  `Glucose Metabolism` = c("glucose", "glycol", "glucokinase", "hexokinase", "pyruvate", 
                           "glycogen", "pentose", "gluconeogenesis", "g6p", "carbohydrate"),
  
  `Secretion & Exocytosis` = c("secretion", "exocytosis", "vesicle", "granule", "hormone secretion",
                               "insulin secretion", "regulated secretion", "peptide hormone"),
  
  `Mitochondrial & Energy` = c("mitochondr", "oxidative phosphorylation", "oxphos", "respiratory chain",
                               "atp synthesis", "electron transport", "energy", "respiration"),
  
  `Ion Channels & Transport` = c("calcium", "potassium", "sodium", "ion channel", "ion transport",
                                 "membrane potential", "depolarization", "channel"),
  
  `ER Stress & Protein Processing` = c("endoplasmic reticulum", "er stress", "unfolded protein",
                                       "protein folding", "upr", "chaperone"),
  
  `Cell Survival & Apoptosis` = c("apoptosis", "cell death", "survival", "necrosis", "autophagy",
                                  "programmed cell death", "caspase"),
  
  `Stress Response` = c("response to stress", "oxidative stress", "hypoxia", "inflammation",
                        "response to oxygen"),
  
  `Signaling` = c("signal transduction", "signaling", "cascade", "receptor", "kinase"),
  
  `Transcription & Gene Expression` = c("transcription", "gene expression", "chromatin", "rna",
                                        "mrna", "translation")
)

# Function to categorize pathways
categorize_pathway <- function(term_name) {
  term_lower <- tolower(term_name)
  
  for (category in names(beta_cell_keywords)) {
    keywords <- beta_cell_keywords[[category]]
    if (any(sapply(keywords, function(kw) grepl(kw, term_lower, fixed = FALSE)))) {
      return(category)
    }
  }
  return(NA)
}

# Categorize all pathways
all_go$category <- sapply(all_go$term_name, categorize_pathway)

# Filter for beta cell-related pathways
beta_pathways <- all_go %>%
  filter(!is.na(category)) %>%
  filter(significant == TRUE) %>%
  select(term_name, p_value, regulation, category, term_id) %>%
  mutate(neg_log10_pval = -log10(p_value)) %>%
  arrange(regulation, desc(neg_log10_pval))

# Print summary
cat("\n=== BETA CELL PATHWAY SUMMARY ===\n")
cat("Total beta cell-related pathways found:", nrow(beta_pathways), "\n")
cat("\nBreakdown by regulation:\n")
print(table(beta_pathways$regulation))
cat("\nBreakdown by category:\n")
print(table(beta_pathways$category, beta_pathways$regulation))

# Select top pathways for visualization (top 25 from each direction)
top_down <- beta_pathways %>%
  filter(regulation == "DOWN") %>%
  slice_max(order_by = neg_log10_pval, n = 25)

top_up <- beta_pathways %>%
  filter(regulation == "UP") %>%
  slice_max(order_by = neg_log10_pval, n = 25)

plot_data <- rbind(top_down, top_up)

# Prepare data for plotting - FIX for duplicate names
plot_data <- plot_data %>%
  mutate(
    # Make DOWN values negative for diverging plot
    plot_value = ifelse(regulation == "DOWN", -neg_log10_pval, neg_log10_pval),
    # Shorten long pathway names
    term_short = ifelse(nchar(term_name) > 50, 
                        paste0(substr(term_name, 1, 47), "..."), 
                        term_name)
  ) %>%
  arrange(plot_value)

# Make shortened names unique by adding numbers to duplicates
plot_data$term_short <- make.unique(plot_data$term_short, sep = " ")

# Now convert to factor with unique levels
plot_data <- plot_data %>%
  mutate(term_short = factor(term_short, levels = term_short))

# Define colors for categories
category_colors <- c(
  `Beta Cell Function` = "#D62728",        # Strong Red - CORE FUNCTION
  `Glucose Metabolism` = "#1F77B4",        # Deep Blue - METABOLISM
  `Secretion & Exocytosis` = "#FF7F0E",    # Orange - SECRETION
  `Mitochondrial & Energy` = "#2CA02C",    # Green - ENERGY/MITO
  `Ion Channels & Transport` = "#9467BD",  # Purple - CHANNELS
  `ER Stress & Protein Processing` = "#8C564B",  # Brown - ER STRESS
  `Cell Survival & Apoptosis` = "#E377C2", # Pink - CELL DEATH
  `Stress Response` = "#BCBD22",           # Olive - STRESS
  `Signaling` = "#17BECF",                 # Cyan - SIGNALING
  `Transcription & Gene Expression` = "#7F7F7F"  # Gray - TRANSCRIPTION
)

# Create the plot
p1 <- ggplot(plot_data, aes(x = term_short, y = plot_value, fill = category)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = category_colors,
                    name = "Pathway Category") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "Beta Cell-Related Pathways in Dynorphin KO",
    subtitle = "Top Pathways Related to Beta Cell Function, Metabolism, and Health",
    x = NULL,
    y = "← Downregulated  -log10(p-value)  Upregulated →"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save the plot
ggsave("/mnt/user-data/outputs/beta_cell_pathways_diverging.png", 
       plot = p1, width = 14, height = 12, dpi = 300)

cat("\n✓ Diverging plot saved!\n")

# Create separate UP and DOWN plots for clarity
# DOWN pathways
top_down_plot <- top_down %>%
  mutate(term_short = ifelse(nchar(term_name) > 50, 
                             paste0(substr(term_name, 1, 47), "..."), 
                             term_name))

# Make unique
top_down_plot$term_short <- make.unique(top_down_plot$term_short, sep = " ")

p_down <- ggplot(top_down_plot, aes(x = reorder(term_short, neg_log10_pval), 
                                    y = neg_log10_pval, fill = category)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = category_colors, name = "Pathway Category") +
  labs(
    title = "Downregulated Beta Cell Pathways in Dynorphin KO",
    x = NULL,
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave("/mnt/user-data/outputs/beta_cell_pathways_DOWN.png", 
       plot = p_down, width = 12, height = 10, dpi = 300)

# UP pathways
top_up_plot <- top_up %>%
  mutate(term_short = ifelse(nchar(term_name) > 50, 
                             paste0(substr(term_name, 1, 47), "..."), 
                             term_name))

# Make unique
top_up_plot$term_short <- make.unique(top_up_plot$term_short, sep = " ")

p_up <- ggplot(top_up_plot, aes(x = reorder(term_short, neg_log10_pval), 
                                y = neg_log10_pval, fill = category)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = category_colors, name = "Pathway Category") +
  labs(
    title = "Upregulated Beta Cell Pathways in Dynorphin KO",
    x = NULL,
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave("/mnt/user-data/outputs/beta_cell_pathways_UP.png", 
       plot = p_up, width = 12, height = 10, dpi = 300)

cat("✓ Separate UP/DOWN plots saved!\n")

# Save the filtered pathways to CSV
write.csv(beta_pathways, 
          "/mnt/user-data/outputs/beta_cell_related_pathways_all.csv", 
          row.names = FALSE)

write.csv(plot_data, 
          "/mnt/user-data/outputs/beta_cell_pathways_plotted.csv", 
          row.names = FALSE)

cat("\n✓ CSV files saved!\n")

# Print detailed list of pathways
cat("\n=== TOP DOWNREGULATED BETA CELL PATHWAYS ===\n")
top_down_display <- top_down %>%
  select(term_name, category, p_value, neg_log10_pval) %>%
  mutate(p_value = format(p_value, scientific = TRUE, digits = 3))
print(top_down_display, n = 25)

cat("\n=== TOP UPREGULATED BETA CELL PATHWAYS ===\n")
top_up_display <- top_up %>%
  select(term_name, category, p_value, neg_log10_pval) %>%
  mutate(p_value = format(p_value, scientific = TRUE, digits = 3))
print(top_up_display, n = 25)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Plots saved to /mnt/user-data/outputs/\n")
cat("- beta_cell_pathways_diverging.png\n")
cat("- beta_cell_pathways_DOWN.png\n")
cat("- beta_cell_pathways_UP.png\n")
cat("- beta_cell_related_pathways_all.csv\n")
cat("- beta_cell_pathways_plotted.csv\n")
