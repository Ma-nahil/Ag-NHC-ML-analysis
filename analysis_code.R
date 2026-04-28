##Artificial intelligence -based structural analysis of Ag-NHC complexes 
#Required libraries 
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(factoextra)
library(pheatmap)
#Import Dataset
#read Excel file
data <- read_excel("ligand_complex_dataset.xlsx")
#Data Reshaping 
long_data <- data %>%
  pivot_longer(cols = c(CN, CC, C_N, C_C),
               names_to = "BondType",
               values_to = "BondLength")
#Faceted Boxplots (Bond-wise comparison)
ggplot(long_data, aes(x = Type, y = BondLength, fill = Type)) +
  geom_boxplot() +
  facet_nested(~ BondType + Type, scales = "free_y") +
  scale_fill_manual(values = c("ligand" = "skyblue", "complex" = "darkorange")) +
  theme_bw() +
  labs(title = "Bond Length Comparison by Type and Feature",
       x = "Molecule Type",
       y = "Bond Length (Å)")
#PCA: Global Bond Structure Separation
scaled_data <- data %>%
  select(CN, CC, C_N, C_C, Type)   # Keep only features + Type
pca <- prcomp(scaled_data[,1:4], center = TRUE, scale. = TRUE)
fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = scaled_data$Type,
             palette = c("skyblue", "darkorange"),
             addEllipses = TRUE,
             legend.title = "Type",
             title = "PCA: Ligand vs Complex (Bond Structure)")
#Heatmap: Global Bond Pattern
heat_matrix <- as.matrix(scaled_data[,1:4])
rownames(heat_matrix) <- data$Complex_ID 
annotation <- data.frame(Type = data$Type)
rownames(annotation) <- data$Complex_ID
pheatmap(heat_matrix,
         annotation_row = annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of Scaled Bond Lengths")
#Individual Bond Comparisons
#C=N Bond
ggplot(data, aes(x = Type, y = CN, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = c("ligand" = "skyblue", "complex" = "darkorange")) +
  theme_bw() +
  labs(title = "Bond Length Comparison: C=N (Ligand vs Complex)",
       x = "Molecule Type",
       y = "C=N Bond Length (Å)")
#C=C Bond
ggplot(data, aes(x = Type, y = CC, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = c("ligand" = "skyblue", "complex" = "darkorange")) +
  theme_bw() +
  labs(title = "Bond Length Comparison: C=C (Ligand vs Complex)",
       x = "Molecule Type",
       y = "C=C Bond Length (Å)")
#C–C Bond
ggplot(data, aes(x = Type, y = C_C, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = c("ligand" = "skyblue", "complex" = "darkorange")) +
  theme_bw() +
  labs(title = "Bond Length Comparison: C–C (Ligand vs Complex)",
       x = "Molecule Type",
       y = "C–C Bond Length (Å)")
#C–N Bond
ggplot(data, aes(x = Type, y = C_N, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = c("ligand" = "skyblue", "complex" = "darkorange")) +
  theme_bw() +
  labs(title = "Bond Length Comparison: C–N (Ligand vs Complex)",
       x = "Molecule Type",
       y = "C–N Bond Length (Å)")




