# This script performs correlation and network analysis on publicly available 
# CF respiratory microbiome data from the following studies:
# https://doi.org/10.1128/jcm.00432-15
# https://doi.org/10.1128/jcm.00432-15
# https://doi.org/10.1371/journal.ppat.1006798
# https://doi.org/10.1007/s12223-017-0562-3
# https://doi.org/10.1371/journal.pone.0172811
# https://doi.org/10.1016/j.heliyon.2018.e00795
# https://doi.org/10.1186/s40168-019-0636-3
# https://doi.org/10.3390/microorganisms9030492
# https://doi.org/10.1186/s40168-020-00810-3
# https://doi.org/10.1016/j.jcf.2023.04.017
# https://doi.org/10.1165/rcmb.2021-0359OC
# https://doi.org/10.3389/fmicb.2020.01463
# https://doi.org/10.1093/infdis/jiaa655
# https://www.proquest.com/dissertations-theses/prevotella-phylogeny-genomic-molecular-insights/docview/2694431131/se-2
# https://doi.org/10.1016/j.jcf.2021.11.003
# https://doi.org/10.1128/spectrum.00787-24
# https://doi.org/10.1172/JCI167957
# https://doi.org/10.1002/alr.23288
# This analysis computes Spearman correlation matrices, generates heatmaps,
# and builds correlation networks with community detection.

# Load packages
library(igraph)
library(Hmisc)
library(corrplot)
library(readxl)
library(RColorBrewer)

# Define list of input files
input_files <- c(
  "PostETI_Adults_Raw_RespOnly.xlsx",
  "PreETI_Peds_Raw_RespOnly.xlsx",
  "PreETI_Adults_Raw_RespOnly.xlsx",
  "PostETI_Sputum_Raw_RespOnly.xlsx",
  "PostETI_Sinus_Raw_RespOnly.xlsx",
  "PostETI_Oro_Raw_RespOnly.xlsx",
  "PreETI_Oro_Raw_RespOnly.xlsx",
  "PreETI_Sputum_Raw_RespOnly.xlsx",
  "PreETI_BAL_Raw_RespOnly.xlsx",
  "PostETI_Overall_Raw_RespOnly.xlsx",
  "PreETI_Overall_Raw_RespOnly.xlsx")

# Create output directory
dir.create("NetAnalysis_Results")

# Loop over each file
for (file in input_files) {
  cat("Processing:", file, "\n")
  
  # Extract file base name without extension
  base_name <- tools::file_path_sans_ext(basename(file))
  output_dir <- file.path("NetAnalysis_Results", base_name)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Load and prepare data
  counts <- read_excel(file)
  counts <- as.data.frame(counts)
  rownames(counts) <- make.unique(as.character(counts[,1]))
  counts <- counts[,-1]
  # Remove rows with all zeros
  counts <- counts[rowSums(counts) > 0, ]
  # Convert to relative abundance by column sums
  col_sums <- colSums(counts)
  rel_abund <- sweep(counts, 2, col_sums + 1e-6, FUN = "/")
  # Sum relative abundance per genus and rank
  genus_totals <- rowSums(rel_abund)
  top15_genera <- names(sort(genus_totals, decreasing = TRUE))[1:min(15, length(genus_totals))]
  # Subset top 15 genera only
  rel_abund_top15 <- rel_abund[top15_genera, ]
  # Transpose data
  rel_abund_top15 <- t(rel_abund_top15)
  
  # Correlation analysis
  cor_mat <- rcorr(as.matrix(rel_abund_top15), type = "spearman")
  cor_mat_sig <- ifelse(cor_mat$P < 0.05, cor_mat$r, 0)
  cor_mat_sig[is.na(cor_mat_sig)] <- 1
  cor_mat_sig[!(cor_mat_sig > 0.5 | cor_mat_sig < -0.5)] <- 0
  
  # Save correlation matrix
  write.csv(cor_mat_sig, file.path(output_dir, "CorrelationMatrix_Filtered.csv"))
  
  # Correlation heatmap
  png(file.path(output_dir, "Correlation_Heatmap.png"),
      width = 5, height = 5, units = "in", res = 600)
  col3 <- colorRampPalette(c("purple", "white", "orange"))
  corrplot(as.matrix(cor_mat_sig),
           col = col3(20),
           tl.cex = 0.8,
           method = "circle",
           tl.col = "black",
           is.corr = TRUE)
  dev.off()
  
  # Build correlation network
  adj_mat_filtered <- cor_mat_sig
  diag(adj_mat_filtered) <- 0
  net <- graph_from_adjacency_matrix(adj_mat_filtered,
                                     mode = "undirected",
                                     diag = FALSE,
                                     weighted = TRUE)
  net <- delete_edges(net, E(net)[weight == 0])
  
  # Edge and node attributes
  E(net)$abs_weight <- abs(E(net)$weight)
  E(net)$color <- ifelse(E(net)$weight > 0, "orange", "purple")   # orange = positive, purple = negative
  E(net)$width <- E(net)$abs_weight * 5
  V(net)$size <- degree(net) * 2 + 5
  
  # Layout using absolute weights
  layout_net <- layout_with_fr(net, weights = E(net)$abs_weight)
  layout_net <- layout_net + matrix(runif(length(layout_net), -0.1, 0.1), ncol = 2)
  
  # Plot correlation network
  png(file.path(output_dir, "Correlation_Network.png"),
      width = 5, height = 5, units = "in", res = 600)
  
  plot(net,
       layout = layout_net,
       vertex.color = "white",
       vertex.label.color = "black",
       vertex.label.cex = 0.6,
       edge.width = E(net)$width,
       edge.color = E(net)$color,
       vertex.size = V(net)$size,
       main = "Correlation Network")
  
  legend("topright", legend = c("Positive", "Negative"),
         col = c("orange", "purple"), lwd = 2, bty = "n")
  
  dev.off()
  
  # Network metrics
  net_density <- edge_density(net, loops = FALSE)
  write(paste("Network density:", net_density),
        file = file.path(output_dir, "Network_Density.txt"))
  
  top_degree <- sort(degree(net), decreasing = TRUE)[1:min(20, vcount(net))]
  
  png(file.path(output_dir, "top15_genera_by_Degree.png"),
      width = 5, height = 5, units = "in", res = 600)
  barplot(top_degree,
          ylab = "Degree",
          cex.axis = 0.8,
          cex.names = 0.8,
          col = "orange",
          las = 2,
          main = "Top Genera by Degree")
  dev.off()
  
  # Community detection with Louvain
  ceb <- cluster_louvain(net, weights = E(net)$abs_weight)
  
  # Assign node colors by community
  n_clusters <- length(unique(membership(ceb)))
  V(net)$color <- brewer.pal(min(n_clusters, 12), "Set3")[membership(ceb)]
  
  # Plot community detection
  png(file.path(output_dir, "Community_Detection.png"),
      width = 5, height = 5, units = "in", res = 600)
  
  plot(net,
       layout = layout_net,
       vertex.label.color = "black",
       vertex.label.cex = 0.9,
       vertex.label.font = 3,
       edge.width = E(net)$width,
       edge.color = E(net)$color,
       vertex.size = V(net)$size)
  
  dev.off()
}
