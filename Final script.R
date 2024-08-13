# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(readxl)

# Set file paths
series_matrix_file <- "C:/Users/User/Desktop/Nineesha/Plan/Research/INFG/Parkinsons/Seriesmatricx.xlsx"
data_dir <- "C:/Users/User/Desktop/Nineesha/Plan/Research/INFG/Parkinsons/plots"

# Ensure the directory exists
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Load the Series Matrix File (Excel format)
exprs_data <- read_excel(series_matrix_file)

# Convert to a base data.frame
exprs_data <- as.data.frame(exprs_data)

# Set row names as Probe IDs
rownames(exprs_data) <- exprs_data$ID_REF

# Remove the ID_REF column from the data frame
exprs_data_numeric <- exprs_data[, -1]

# Convert all columns to numeric
exprs_data_numeric[] <- lapply(exprs_data_numeric, as.numeric)

# Log2 normalization
exprs_data_log2 <- log2(exprs_data_numeric + 1)

# List of genes we want
ifng_genes <- c("Stat1", "Irf7", "Irf1", "Cxcl10", "Isg15", "Il6", "Il1b", "Irf9", "Stat2", "Gbp3", "Psmb8", "Psmb9", "Gbp2", "Dhx58", "Tapbp")

# Mapping probe IDs to gene names
probe_to_gene <- c("1007_s_at" = "Stat1", "1053_at" = "Irf7", "117_at" = "Irf1",
                   "121_at" = "Cxcl10", "1255_g_at" = "Isg15", "1294_at" = "Il6",
                   "1316_at" = "Il1b", "1320_at" = "Irf9", "1405_i_at" = "Stat2",
                   "1431_at" = "Gbp3", "1438_at" = "Psmb8", "1487_at" = "Psmb9")

# Filter expression data for the specific genes
matching_probes <- names(probe_to_gene)[probe_to_gene %in% ifng_genes]
filtered_data <- exprs_data_log2[rownames(exprs_data_log2) %in% matching_probes, ]

# Create and save box plots for each gene
for (gene in ifng_genes) {
  gene_probes <- names(probe_to_gene)[probe_to_gene == gene]
  gene_data <- filtered_data[rownames(filtered_data) %in% gene_probes, ]
  melted_data <- melt(as.matrix(gene_data))
  
  # Create the box plot
  p <- ggplot(data = melted_data, aes(x = Var2, y = value)) +
    geom_boxplot(aes(fill = Var2)) +
    ggtitle(paste("Box Plot for", gene)) +
    xlab("Samples") + ylab("Expression Level (Log2)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot to a file
  ggsave(file = file.path(data_dir, paste0(gene, "_box_plot.png")), plot = p)
}

# Create and save comprehensive heatmap
png(filename = file.path(data_dir, "heatmap.png"))
pheatmap(filtered_data, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Comprehensive Heatmap of IFN-?? Signaling Genes")
dev.off()

     