# This script generates prevalence–abundance plots of publicly available CF respiratory microbiome data
# from the following studies:
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

# Load required packages
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Generate prevalence–abundance plot
generate_abundance_plot <- function(
    file_path,
    title,
    label_genera = c(
      "Pseudomonas", "Prevotella", "Staphylococcus", "Streptococcus",
      "Fusobacterium", "Rothia", "Haemophilus", "Veillonella", "Porphyromonas"
    )
) {
  # Read and clean
  df <- as.data.frame(read_excel(file_path, sheet = 1))
  rownames(df) <- df[[1]]
  df <- df[, -1]
  df <- df[rowSums(df, na.rm = TRUE) > 0, , drop = FALSE]
  df <- df[, colSums(df, na.rm = TRUE) > 0, drop = FALSE]
  
  # Normalize
  rel_abund <- apply(df, 2, function(s) s / sum(s))
  
  # Compute prevalence and mean abundance
  prevalence <- apply(rel_abund, 1, function(x) sum(x > 0) / length(x))
  mean_abund <- apply(rel_abund, 1, mean)
  
  # Filter
  filtered <- data.frame(
    Genus = names(prevalence[prevalence > 0.2]),
    Prev = prevalence[prevalence > 0.2],
    Ab = log10(pmax(mean_abund[prevalence > 0.2]))
  )
  
  csv_name <- paste0(tools::file_path_sans_ext(basename(file_path)), "_PrevAbnd_data.csv")
  write.csv(filtered, csv_name, row.names = FALSE)
  
  # Select top 15 genera
  top15 <- filtered %>%
    arrange(desc(Ab)) %>%
    head(15)
  
  top15$Genus <- factor(top15$Genus, levels = top15$Genus)
  label_data <- top15 %>% filter(Genus %in% label_genera)
  
  # Plot prevalence vs abundance
  p <- ggplot(top15, aes(x = Prev * 100, y = Ab)) +
    geom_point(aes(color = Genus), size = 7) +
    geom_text_repel(
      data = label_data,
      aes(label = Genus),
      size = 8,
      max.overlaps = 15,
      box.padding = 1,
      point.padding = 1,
      fontface = "italic"
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x)
    ) +
    scale_y_continuous(
      limits = c(-3, -0.2),
      breaks = seq(-3, -0.2, 0.25),
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 25),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25, face = "italic")
    ) +
    labs(
      title = title,
      x = "Prevalence (%)",
      y = "Relative Abundance (log10)",
      color = "Top Genera"
    )
  return(p)
}

  files <- list(
  "PreETI_Overall_Raw_RespOnly.xlsx" = "Pre-ETI: All samples",
  "PreETI_Sputum_Raw_RespOnly.xlsx"  = "Pre-ETI: Sputum samples",
  "PreETI_BAL_Raw_RespOnly.xlsx"     = "Pre-ETI: BAL samples",
  "PreETI_Oro_Raw_RespOnly.xlsx"     = "Pre-ETI: Oropharynx samples",
  "PreETI_Adults_Raw_RespOnly.xlsx"  = "Pre-ETI: Adult patients",
  "PreETI_Peds_Raw_RespOnly.xlsx"    = "Pre-ETI: Pediatric patients",
  "PostETI_Overall_Raw_RespOnly.xlsx" = "Post-ETI: All samples",
  "PostETI_Adults_Raw_RespOnly.xlsx" = "Post-ETI: Adult patients",
  "PostETI_Oro_Raw_RespOnly.xlsx"    = "Post-ETI: Oropharynx samples",
  "PostETI_Sinus_Raw_RespOnly.xlsx"  = "Post-ETI: Sinonasal samples",
  "PostETI_Sputum_Raw_RespOnly.xlsx" = "Post-ETI: Sputum samples")


  plots <- lapply(names(files), function(f) generate_abundance_plot(f, files[[f]]))
  names(plots) <- names(files)

  for (i in seq_along(plots)) {
  print(plots[[i]])
  }
