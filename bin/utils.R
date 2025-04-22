# ---
# author: "Laura Jiménez Gracia"
# date: 2021-03-20
# ---
# This R script is used to define functions used in multiple scripts.

################################################################################
############################ QC Cellranger mapping #############################
################################################################################

clean_cellranger_metrics_df <- function(metrics_df) {
  colnames(metrics_df) <- str_remove(colnames(metrics_df), ":")
  colnames(metrics_df) <- str_remove(colnames(metrics_df), "\\(")
  colnames(metrics_df) <- str_remove(colnames(metrics_df), "\\)")
  colnames(metrics_df) <- str_replace_all(colnames(metrics_df), "-", "_")
  metrics_df <- as.data.frame(metrics_df)
  for (col in colnames(metrics_df)) {
    metrics_df[, col] <- str_remove(metrics_df[, col], regex("_\\(\\d+.\\d+\\%\\)"))
    if (any(str_detect(na.omit(metrics_df[, col]), "%"))) {
      metrics_df[, col] <- as.double(str_remove(metrics_df[, col], "%"))
    }
  }
  metrics_df
}

group_hashing_cellranger_metrics_df <- function(
  metrics_df, not_hashed_gemid_list, hashed_gemid_list) {
  
  not_hashed_df <- subset.data.frame(metrics_df,
                                     subset = gem_id %in% not_hashed_gemid_list)
  hashed_df <- subset.data.frame(metrics_df,
                                 subset = gem_id %in% hashed_gemid_list)
  metrics_df <- bind_rows(list("not_hashed" = not_hashed_df,
                               "CellPlex" = hashed_df),
                          .id = "hashing")
  metrics_df <- metrics_df[order(metrics_df$gem_id), ]
  metrics_df
}


table_cellranger_metrics_gex_lib_hash <- function(metrics_gex_df, cellranger_version) {
  table_metrics_gex <- metrics_gex_df[, c("hashing",
                                          "library_name",
                                          "Number_of_reads_in_the_library", 
                                          "Estimated_number_of_cells", 
                                          "Fraction_reads_in_cells",
                                          "Mean_reads_per_cell",
                                          "Confidently_mapped_to_exonic_regions"
                                          )]
  
  table_metrics_gex %>%
    gt::gt() %>%
    fmt_percent(columns = c("Fraction_reads_in_cells", "Confidently_mapped_to_exonic_regions"), 
                scale_values = FALSE, decimals = 1) %>%
    fmt_number(columns = "Number_of_reads_in_the_library", scale_by = 1 / 1E6, pattern = "{x}M") %>% 
    tab_header(
      title = md("**Library GEX QC metrics**"),
      subtitle = (cellranger_version)
    ) %>%
    cols_label(
      hashing = md("**Hashing**"),
      library_name = md("**Library**"),
      Number_of_reads_in_the_library = md("**Number of Reads**"),
      Estimated_number_of_cells = md("**Estimated Number of Recovered Cells**"),
      Fraction_reads_in_cells = md("**Fraction of Reads in Cells**"),
      Mean_reads_per_cell = md("**Mean Reads per Cell**"),
      Confidently_mapped_to_exonic_regions = md("**Fraction of Reads Mapped to Exonic Reads**")
      )  
}


table_cellranger_metrics_hashing <- function(metrics_gex_df, cellranger_version) {
  table_metrics_hashing <- metrics_gex_df[, c("library_name",
                                              "Number_of_reads_CMO",
                                              "Mean_reads_per_cell_associated_barcode_CMO",
                                              "Fraction_CMO_reads_usable_CMO",
                                              "Median_CMO_UMIs_per_cell_associated_barcode_CMO")]
  
  table_metrics_hashing %>%
    gt::gt() %>%
    fmt_percent(columns = c("Fraction_CMO_reads_usable_CMO"), 
                scale_values = FALSE, decimals = 1) %>%
    fmt_number(columns = "Number_of_reads_CMO", scale_by = 1 / 1E6, pattern = "{x}M") %>% 
    tab_header(
      title = md("**CMO QC metrics**"),
      subtitle = (cellranger_version)
    ) %>%
    cols_label(
      library_name = md("**Library**"),
      Number_of_reads_CMO = md("**Number of CMO Reads**"),
      Mean_reads_per_cell_associated_barcode_CMO = md("**Mean CMO Reads per Cell**"),
      Fraction_CMO_reads_usable_CMO = md("**Fraction of CMO Reads (usable)**"),
      Median_CMO_UMIs_per_cell_associated_barcode_CMO = md("**Median CMO UMIs per Cell**")
    )
}


table_demultiplexing_summary <- function(metrics_gex_df_hashed, cellranger_version) {
  
  demux_vars <- c("library_name",
                  "Estimated_number_of_cells",
                  "Cells_assigned_to_a_sample",
                  "Cell_associated_barcodes_identified_as_multiplets",
                  "Cell_associated_barcodes_not_assigned_any_CMOs")
  
  demux_global_cells <- metrics_gex_df_hashed[, demux_vars]
  
  demux_global_df <- demux_global_cells %>%
    mutate(Singlet = Cells_assigned_to_a_sample) %>% 
    mutate(Singlet_pct = Cells_assigned_to_a_sample/Estimated_number_of_cells) %>%
    mutate(Multiplet = Cell_associated_barcodes_identified_as_multiplets) %>% 
    mutate(Multiplet_pct = Cell_associated_barcodes_identified_as_multiplets/Estimated_number_of_cells) %>%    
    mutate(Unassigned = Cell_associated_barcodes_not_assigned_any_CMOs) %>% 
    mutate(Unassigned_pct = Cell_associated_barcodes_not_assigned_any_CMOs/Estimated_number_of_cells) %>%     
    mutate(Multiplet_rate = Estimated_number_of_cells*0.08/10000)
  demux_global_df <- demux_global_df[c("library_name", "Estimated_number_of_cells", 
                                       "Singlet", "Singlet_pct",
                                       "Multiplet", "Multiplet_pct", "Multiplet_rate",
                                       "Unassigned", "Unassigned_pct")]
  rownames(demux_global_df) <- 1:nrow(demux_global_df)
  
  demux_global_df %>% 
    gt() %>%
    fmt_percent(
      columns = vars(Multiplet_pct, Multiplet_rate, Singlet_pct, Unassigned_pct),
      decimals = 1
    ) %>% 
    tab_header(
      title = md("**CMO Multiplexing QC**"),
      subtitle = (cellranger_version)
    ) %>%
    cols_label(
      library_name = md("**Library**"),
      Estimated_number_of_cells = md("**# Total cells**"),
      Multiplet = md("**# cells**"),
      Multiplet_pct = md("**Fraction**"),
      Multiplet_rate = md("*Expected*"),
      Singlet = md("**# cells**"),
      Singlet_pct = md("**Fraction**"),
      Unassigned = md("**# cells**"),
      Unassigned_pct = md("**Fraction**"),
    ) %>% 
    tab_spanner(
      label = md("**Singlets**"),
      columns = vars(
        Singlet,
        Singlet_pct)) %>% 
    tab_spanner(
      label = md("**Multiplets**"),
      columns = vars(
        Multiplet,
        Multiplet_pct,
        Multiplet_rate))  %>% 
    tab_spanner(
      label = md("**Unassigned**"),
      columns = vars(
        Unassigned,
        Unassigned_pct)
    )
}

table_demultiplexing_CMO <- function(metrics_demux_df, cellranger_version) {
  metrics_demux_df_sub <- metrics_demux_df[c("library_name", "sample_name", 
                                             "CMO_name", "CMO_signal_to_noise_ratio",
                                             "Cells_assigned_to_CMO", "Cells_assigned_to_this_sample"
  )]
  metrics_demux_df_sub %>% 
    gt(rowname_col = "row", groupname_col = "library_name") %>%
    tab_header(
      title = md("**CMO signal-to-noise QC**"),
      subtitle = (cellranger_version=cellranger_version)
    ) %>%
    fmt_percent(columns = c("Cells_assigned_to_CMO"), 
                scale_values = FALSE, decimals = 1) %>%
    cols_label(
      #library_name = md("**Library**"),
      sample_name = md("**Sample**"),
      CMO_name = md("**CMO**"),
      CMO_signal_to_noise_ratio = md("**CMO Signal-to-Noise Ratio**"),
      Cells_assigned_to_this_sample = md("**# cells**"),
      Cells_assigned_to_CMO = md("**Fraction**")
    ) %>% 
    tab_spanner(
      label = md("**Singlets assigned to sample**"),
      columns = vars(
        Cells_assigned_to_this_sample,
        Cells_assigned_to_CMO))
}



table_cellranger_metrics_gex_sample_hash <- function(metrics_gex_df, cellranger_version) {
  table_metrics_demux <- metrics_demux_df[, c("library_name",
                                              "sample_name",
                                              "Number_of_reads_assigned_to_the_sample", 
                                              "Cells_assigned_to_this_sample", 
                                              "Fraction_reads_in_cell_associated_barcodes",
                                              "Median_reads_per_cell",
                                              "Confidently_mapped_to_exonic_regions",
                                              "Median_genes_per_cell"
  )]
  
  table_metrics_demux %>%
    gt::gt(rowname_col = "row", groupname_col = "library_name") %>%
    fmt_percent(columns = c("Fraction_reads_in_cell_associated_barcodes", "Confidently_mapped_to_exonic_regions"), 
                scale_values = FALSE, decimals = 1) %>%
    fmt_number(columns = "Number_of_reads_assigned_to_the_sample", scale_by = 1 / 1E6, pattern = "{x}M") %>% 
    tab_header(
      title = md("**Sample GEX QC metrics**"),
      subtitle = (cellranger_version=cellranger_version)
    ) %>%
    cols_label(
      library_name = md("**Library**"),
      sample_name = md("**Sample**"),
      Number_of_reads_assigned_to_the_sample = md("**Number of Reads**"),
      Cells_assigned_to_this_sample = md("**Cells assigned to sample**"),
      Fraction_reads_in_cell_associated_barcodes = md("**Fraction of Reads in Cells**"),
      Median_reads_per_cell = md("**Median Reads per Cell**"),
      Confidently_mapped_to_exonic_regions = md("**Fraction of Reads Mapped to Exonic Reads**"),
      Median_genes_per_cell = md("**Median Genes per Cell**")
    )  
}



gg_demultiplex_percent_glob <- function(dataframe) {
  ###
  # Plots the percentage of singlets, multiplets and unassigned cell-barcodes
  ###
  dataframe[, c("library_name", "Estimated_number_of_cells", "Cells_assigned_to_a_sample", "Cell_associated_barcodes_identified_as_multiplets", "Cell_associated_barcodes_not_assigned_any_CMOs")] %>% 
    mutate(Singlets = Cells_assigned_to_a_sample/Estimated_number_of_cells*100,
           Multiplets = Cell_associated_barcodes_identified_as_multiplets/Estimated_number_of_cells*100,   
           Unassigned_cells = Cell_associated_barcodes_not_assigned_any_CMOs/Estimated_number_of_cells*100) %>%  
    pivot_longer(c(Singlets, Multiplets, Unassigned_cells), names_to = "demux_classification", values_to = "pct") %>% 
    ggplot(aes(x = "", y = pct, fill = demux_classification)) +
    geom_bar(width = 1, stat = "identity") +
    facet_grid(~library_name, scales = "free_x", space = "free") +
    geom_text(
      aes(label = paste0(pct * Estimated_number_of_cells/100, " (", percent(pct / 100, accuracy = 0.1), ")")),
      position = position_stack(vjust = 0.5),
      size = 5
    ) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = "Demultiplexing: Singlets, Multiplets and Unassigned cells",
         x = "Libraries",
         y = "Fraction (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 10),
      strip.text.x = element_text(size = 16),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
    )
}

gg_demultiplex_sample <- function(dataframe) {
  ###
  # Plots the number / percentage of sample-specific cell-barcodes within library
  ###
  metrics_demux_df %>%
    ggplot(aes(x = "", y = Cells_assigned_to_CMO, fill = sample_name)) +
    geom_bar(width = 1, stat = "identity") +
    facet_grid(~library_name, scales = "free_x", space = "free") +
    geom_text(
      aes(label = paste0(Cells_assigned_to_this_sample, " (", round(Cells_assigned_to_CMO, 1), "%)")),
      position = position_stack(vjust = 0.5),
      size = 5
    ) +
    scale_fill_brewer(palette = "Paired") +
    labs(title = "Demultiplexing: Sample specific cells",
         x = "Libraries",
         y = "Fraction (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 10),
      strip.text.x = element_text(size = 16),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
    )
}



################################################################################
################ QC Hashtag Demultiplexing & Doublet prediction ################
################################################################################


gg_doublet_scores_scrublet <- function(seurat_obj_metadata, scrublet_threshold) {
  ###
  # Plot the scrublet doublet score of singlets, doublets and negative 
  #
  # Parameters
  # seurat_obj_metadata: seurat metadata containing the HTO_classification.global and scrublet_doublet_scores
  # scrublet_threshold: Approx. scrublet doublet predicted threshold (by default)
  #
  # Return
  # It returns a ggplot boxplot object with scrublet doublet scores for each HTO global classificator
  ###
  seurat_obj@meta.data %>%
    ggplot(aes(x = HTO_classification.global,
               y = scrublet_doublet_scores,
               fill = HTO_classification.global)) +
    geom_violin() +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.3) +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_brewer(palette = "Set2") +
    geom_hline(yintercept = scrublet_threshold, linetype = "dashed", colour = "red") +
    geom_text(aes(0, scrublet_threshold, label = scrublet_threshold,
                  vjust = -1, size = 3, hjust = -0.4, colour = "red")) +
    labs(x = "",
         y = "Scrublet Doublet Score") +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none"
    )
}

gg_demultiplexing_qc_snr_glob <- function(seurat_obj_metadata, hashed_gemids) {
  seurat_obj_metadata %>%
    ggplot(aes(x = gem_id,
               y = hashing_snr,
               fill = HTO_classification.global)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = "Hashing QC",
         x = "Libraries (GEM id)",
         y = "Signal-to-Noise Ratio",
         fill = "Barcode classifiaction") +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(hjust = 1)
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    coord_flip() +
    scale_x_discrete(limits = rev(hashed_gemids))  
}


gg_demultiplexing_qc_snr_sing <- function(seurat_obj_metadata, hashed_gemids) {
  seurat_obj_metadata %>%
    filter(HTO_classification.global == "Singlet") %>% 
    ggplot(aes(x = gem_id,
               y = hashing_snr,
               fill = HTO_classification.global)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = "Hashing QC  [Singlets]",
         x = "Libraries (GEM id)",
         y = "Signal-to-Noise Ratio") +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(hjust = 1),
          legend.position = "none"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    coord_flip() +
    scale_x_discrete(limits = rev(hashed_gemids))
}




gg_demultiplex_svn_libsize <- function(seurat_metadata, hashed_gemids) {
  seurat_metadata %>%
    dplyr::filter(HTO_classification.global != "Doublet") %>%
    dplyr::mutate(HTO_classification.global = factor(
      HTO_classification.global,
      levels = c("Singlet", "Negative")
    )) %>% 
    ggplot(aes_string(x = "gem_id", y ="nCount_RNA", 
                      fill = "HTO_classification.global")) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = "Library Size",
         x = "Libraries (GEM ID)", 
         y = "Total UMIs",
         fill = "Barcode classification") +
    theme_bw() +
    scale_fill_manual(values = c("#66C2A5", "#8DA0CB")) +
    scale_y_log10() +
    theme(axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          legend.position = "bottom") + 
    coord_flip() +
    scale_x_discrete(limits = rev(hashed_gemids))
}


gg_demultiplex_svn_libcomplex <- function(seurat_metadata, hashed_gemids) {
  seurat_metadata %>%
    dplyr::filter(HTO_classification.global != "Doublet") %>%
    dplyr::mutate(HTO_classification.global = factor(
      HTO_classification.global,
      levels = c("Singlet", "Negative")
    )) %>% 
    ggplot(aes_string(x = "gem_id", y ="nFeature_RNA", 
                      fill = "HTO_classification.global")) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = "Library Complexity",
         x = "", 
         y = "Number of Detected Genes",
         fill = "Barcode classification") +
    theme_bw() +
    scale_fill_manual(values = c("#66C2A5", "#8DA0CB")) +
    scale_y_log10() +
    theme(axis.text.x = element_text(size = 11),
          axis.title.y = element_text(size = 14),
          legend.position = "bottom") + 
    coord_flip() +
    scale_x_discrete(limits = rev(hashed_gemids))
}
  
gg_scrublet_qc_glob <-

gg_scrublet_qc_singlets <- 


gg_demultiplex_scrublet_scores_glob <- function(seurat_metadata) {
  seurat_metadata %>%
    ggplot(aes(x = HTO_classification.global,
               y = scrublet_doublet_scores,
               fill = HTO_classification.global)) +
    geom_violin() +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.3) +
    theme_bw() +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, 1) +
    labs(title = "Scrublet doublet scores",
         x = "Barcode classification",
         y = "Doublet Score") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "none"
    )
}
  
  
gg_demultiplex_scrublet_pred_glob <- function(seurat_metadata) {
  seurat_metadata %>%
    ggplot(aes(x = HTO_classification.global,
               y = scrublet_doublet_scores,
               fill = scrublet_predicted_doublet)) +
    geom_violin() +
    theme_bw() +
    ylim(0, 1) +
    labs(title = "Scrublet predicted doublets",
      x = "Barcode classification",
         y = "Doublet Score",
         fill = "Predicted Doublet") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "bottom"
    )
}


gg_demultiplex_scrublet_scores_sing <- function(seurat_metadata) {
  seurat_metadata %>%
  filter(HTO_classification.global == "Singlet") %>% 
  ggplot(aes(x = scrublet_doublet_scores)) + 
  geom_density() + 
  xlim(0, 1) +
  labs(title = "[Singlets] Scrublet doublet scores",
       x = "Doublet Score",
       y = "Density") +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none"
  )
}


gg_demultiplex_scrublet_pred_sing <- function(seurat_metadata) {
  seurat_metadata %>%
    filter(HTO_classification.global == "Singlet") %>% 
    ggplot(aes(x = scrublet_doublet_scores, fill = scrublet_predicted_doublet)) + 
    geom_density(alpha = 0.4) + 
    xlim(0, 1) +
    labs(title = "[Singlets] Scrublet predicted doublets",
         x = "Doublet Score",
         y = "Density",
         fill = "Predicted Doublet") +
    theme_pubr() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "bottom"
    )
  
}


gg_nothashed_scrublet_pred <- function(seurat_metadata, gem_id, scrublet_threshold) {
  seurat_metadata %>%
    ggplot(aes(x = gem_id,
               y = scrublet_doublet_scores,
               fill = scrublet_predicted_doublet,
               colour = scrublet_predicted_doublet)) +
    geom_violin(alpha = 0.1) +
    geom_jitter(height = 0, width = 0.1) + 
    theme_bw() +
    ylim(0, 1) +    
    geom_hline(yintercept = scrublet_threshold, linetype = "dashed", colour = "red") +
    geom_text(aes(0, scrublet_threshold, label = scrublet_threshold,
                  vjust = -1, hjust = -0.4)) +
    labs(title = gem_id,
         x = "",
         y = "Doublet Score",
         fill = "Predicted Doublet",
         colour = "Predicted Doublet") +
    theme(axis.title = element_text(size = 14),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom"
    )
}

gg_nothashed_umis <- function(seurat_obj) {
  ###
  # Plots the RNA UMIs distribution of singlets, doublets and negative 10X barcodes
  #
  # Parameters
  # seurat_obj
  #
  # Return
  # It returns a ggplot violin plot object with UMI distribution for HTO-global classification
  ###
  Idents(seurat_obj) <- "scrublet_predicted_doublet"
  VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1, log = TRUE) +
    labs(x = "Predicted Doublet", y = "Number of UMIs") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          plot.title = element_blank(),
          legend.position = "none"
    )
}


  
################################################################################
############################## QC Gene Expression ############################## 
################################################################################

gg_gex_horizontal_boxplot <- function(df,
                               categorical_var,
                               continuous_var,
                               fill,
                               title,
                               ylab,
                               decreasing = FALSE) {
  unordered_lev <- unique(df[[categorical_var]])
  means_cont <- purrr::map_dbl(unordered_lev, function(x) {
    mean(df[[continuous_var]][df[[categorical_var]] == x])
  })
  names(means_cont) <- unordered_lev
  ordered_lev <- names(means_cont)[order(means_cont, decreasing = decreasing)]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = ordered_lev)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var, fill = fill)) +
    geom_boxplot() +
    labs(title = title, x = "", y = ylab) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = "bottom") +
    coord_flip() +
    scale_y_log10()
}


table_qc_gex <- function(qc_metadata_df, subtitle) {
  qc_metadata_df %>%
    group_by(sample_id) %>%
    dplyr::summarise(num_cells = n(),
                     mean_library_size = round(mean(nCount_RNA)),
                     median_num_detected_genes = round(median(nFeature_RNA)),
                     average_mitochondrial_fraction = mean(pct_mt) / 100
    ) %>%
    gt() %>%
    tab_header(
      title = md("**GEX QC**"),
      subtitle = subtitle
    ) %>%
    cols_label(
      sample_id = md("**Sample**"),
      num_cells = md("**Number of Cells**"),
      mean_library_size = md("**Mean Reads per Cell**"),
      median_num_detected_genes = md("**Median Genes per Cell**"),
      average_mitochondrial_fraction = md("**Mean Mitochondrial Fraction**")
    ) %>%
    fmt_percent(
      columns = vars(average_mitochondrial_fraction),
      decimals = 1
    )
}


################################################################################
############################# QC Doublet exclusion #############################
########################### by Ramon Massoni-Badosa ############################
################################################################################

plot_histogram_doublets <- function(df, x, x_lab, bins) {
  df %>%
    ggplot(aes_string(x)) +
    geom_histogram(bins = bins) +
    xlab(x_lab) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

plot_density_doublets <- function(df, x, x_lab, color, color_lab) {
  df %>%
    ggplot(aes_string(x = x, color = color)) +
    geom_density() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = x_lab, color = color_lab) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    )
}


plot_boxplot_doublets <- function(df, x, y, fill, y_lab) {
  df %>%
    ggplot(aes_string(x, y, fill = fill)) +
    geom_boxplot() +
    labs(x = "", y = y_lab) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 13),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 11)
    )
}


plot_scatter_doublets <- function(df, x, y, x_lab, y_lab) {
  df %>%
    ggplot(aes_string(x, y)) +
    geom_point(size = 0.1, alpha = 0.25) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = x_lab, y = y_lab) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}


feature_plot_doublets <- function(seurat_obj, feature) {
  p <- FeaturePlot(
    seurat_obj,
    features = feature,
    cols = viridis::inferno(10),
    pt.size = 0.3
  ) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.text = element_text(size = 11)
    )
  p
}


################################################################################
################################# GEX analysis #################################
################################################################################

find_DEgenes <- function(seurat, ident_1, ident_2, test_de, threshold_pvaladj) {
  # Find DE genes
  de_genes <- Seurat::FindMarkers(
    seurat,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = test_de, # default wilcox
    logfc.threshold = 0, # default 0.25
    min.pct = 0.1 # default 0.1
  )
  
  # Tag and sort significant DE genes
  de_genes <- de_genes %>% 
    rownames_to_column(var = "gene") %>%
    mutate(is_significant = ifelse(p_val_adj < threshold_pvaladj, TRUE, FALSE)) %>% 
    dplyr::arrange(desc(is_significant), desc(abs(avg_log2FC)))
  
  # Return DE genes
  de_genes
}


find_DEgenes_celltypes <- function(seurat, ident_1, ident_2, test_de, threshold_pvaladj, min_cells=3, latent_vars=NULL) {
  # Find DE genes
  de_genes <- Seurat::FindMarkers(
    seurat,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = test_de, # default wilcox
    logfc.threshold = 0, # default 0.25
    min.pct = 0.1, # default 0.1
    min.cells.group = min_cells, # default 3
    latent.vars = latent_vars #default NULL
  )
  
  # Tag and sort significant DE genes
  de_genes <- de_genes %>% 
    rownames_to_column(var = "gene") %>%
    mutate(is_significant = ifelse(p_val_adj < threshold_pvaladj, TRUE, FALSE)) %>% 
    dplyr::arrange(desc(is_significant), desc(abs(avg_log2FC)))
  
  # Return DE genes
  de_genes
}


gg_volcano_DEgenes <- function(DEgenes_df, title, ident_1, ident_2, threshold_log2FC, threshold_adjpval) {
  # Pre-process data
  DEgenes_df <- DEgenes_df %>% 
    mutate(is_sig_direction = case_when(
      is_significant & avg_log2FC > threshold_log2FC ~ "UP",
      is_significant & avg_log2FC < -threshold_log2FC ~ "DOWN",
      TRUE ~ "FALSE"))
  
  # Subset sDE gene names
  subset_sig <- DEgenes_df %>% 
    dplyr::filter(is_significant & abs(avg_log2FC) > threshold_log2FC)
  top_up <- as.character(subset_sig$gene[subset_sig$avg_log2FC > 0][1:15])
  top_down <- subset_sig %>% 
    dplyr::arrange(avg_log2FC)
  top_down <- as.character(top_down$gene[1:15])
  subset_sig <- subset_sig %>% 
    dplyr::filter(gene %in% c(top_up, top_down))

  # Volcano plot
  volcano_plot <- DEgenes_df %>%
    ggplot(aes(avg_log2FC, -1 * log10(p_val_adj), color = is_sig_direction)) +
    geom_point(size = 3, alpha=0.7) +
    scale_color_manual(
      values = c("UP" = "red", "FALSE" = "grey", "DOWN" = "blue"),
      labels = c("UP" = "Up-regulated genes",
                 "FALSE" ="non-sDE genes",
                 "DOWN" = "Down-regulated genes"
      )) +      
    geom_hline(yintercept = threshold_adjpval, color = "grey", linetype = "dashed") +
    geom_vline(xintercept = threshold_log2FC, color = "grey", linetype = "dashed") +
    geom_vline(xintercept = -threshold_log2FC, color = "grey", linetype = "dashed") +
    geom_text_repel(data = subset_sig, aes(label = gene), color = "black", size=5, max.overlaps = Inf) +
    labs(title = title,
         x = "Log2FC",
         y = "-Log10(Adj p-value)",
         color = "") +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  

  volcano_plot
}


run_enrichR_forGO <- function(database, DE_genes, adj_pval, number_enriched_genes, num_GO_size_min, num_GO_size_max, OR) {
  
  # Run enrichR
  enrichR_results <- enrichr(DE_genes, databases = database)[[1]]
  
  # Term name
  enrichR_results$term_name <- lapply(enrichR_results[["Term"]], function(x) {
    x <- str_split(x, pattern = "\\(", n = Inf)[[1]][1]
    x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    x
  }
  )
  
  # Term GO code
  enrichR_results$code_GO <- lapply(enrichR_results[["Term"]], function(x) {
    str_extract_all(x, "(?<=\\().*(?=\\))")[[1]]}
  )  
  
  # Number of DE genes enriched for a specific GO
  enrichR_results$number_genes <- lapply(enrichR_results[["Overlap"]], function(x) {
    str_split(x, pattern = "/", n = Inf)[[1]][1]})
  
  # Number of genes within a specific GO
  enrichR_results$GO_size <- lapply(enrichR_results[["Overlap"]], function(x) {
    str_split(x, pattern = "/", n = Inf)[[1]][2]})
  
  # Percentage overlap
  enrichR_results$pct_overlap <- sapply(enrichR_results$Overlap, function(x) 
    eval(parse(text=x)))
  
  # Convert dataframe results
  enrichR_results <- as.tibble(enrichR_results) %>% 
    transform(term_name = as.character(term_name),
              code_GO = as.character(code_GO),
              number_genes = as.numeric(number_genes),
              GO_size = as.numeric(GO_size))
  
  # Tag GO as reliable TRUE/FALSE
  enrichR_results <- enrichR_results %>%
    mutate(is_reliable = ifelse(
      (Adjusted.P.value < adj_pval & 
         number_genes >= number_enriched_genes & 
         GO_size >= num_GO_size_min &
         GO_size <= num_GO_size_max & 
         Odds.Ratio  >= OR), TRUE, FALSE)) %>%
    dplyr::arrange(desc(is_reliable), Adjusted.P.value)
  
  enrichR_results
}

gg_enrichR_GO <- function(comparison, database, enrichR_results_filtered, n_top) {
  
  enrichR_results_filtered <- enrichR_results_filtered[1:n_top,][, c("Odds.Ratio", "term_name", "Adjusted.P.value", "number_genes")]
  
  gg_plot <- enrichR_results_filtered %>%
    ggplot(., aes(
      x = Odds.Ratio,
      y = fct_reorder(term_name, -Adjusted.P.value),
      size = number_genes,
      color = Adjusted.P.value)) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "red") +
    scale_color_gradient2(
      name = "Adjusted p-value",
      low = "#bf212f",
      mid = "#006f3c",
      high = "lightgrey") +
    theme_minimal() +
    labs(title = comparison,
         x = "Odds Ratio",
         y = paste0(database, "Term"),
         size = "# Genes") +
    theme(
      axis.text.y = element_text(size=12)
    )
  
  gg_plot
}






#   # Functions
#   
 
# DEanalysis_MAplot_general <- function(DEgenes_df, ident, ident_1, ident_2, threshold_log2FC, threshold_avexpr) {
#   # Subset sDE gene names
#   deg_subset_sig <- dplyr::filter(DEgenes_df, is_significant & 
#                                     log(average_expression + 1) > threshold_avexpr & 
#                                     abs(log2_fc) > threshold_log2FC)
#   
#   deg_subset_sig_up <- deg_subset_sig %>% 
#     dplyr::arrange(desc(log2_fc)) %>% 
#     filter(log2_fc > 0, preserve = TRUE)
#   deg_subset_sig_up_names <- as.character(
#     slice_max(deg_subset_sig_up, order_by = log2_fc, n = 10)$gene)
#   
#   deg_subset_sig_down <- deg_subset_sig %>% 
#     dplyr::arrange(log2_fc) %>% 
#     filter(log2_fc < 0, preserve = TRUE)
#   deg_subset_sig_down_names <- as.character(
#     slice_min(deg_subset_sig_down, order_by = log2_fc, n = 10)$gene)
#   
#   deg_subset_sig <- deg_subset_sig %>% 
#     dplyr::filter(gene %in% c(deg_subset_sig_up_names, deg_subset_sig_down_names))
#   
#   # MA plot
#   maplot_general <- DEgenes_df %>% 
#     ggplot(aes(log2(average_expression + 1), log2_fc, color = is_significant)) +
#     geom_point(size = 1) +
#     geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
#     geom_smooth(method = "loess", color = "darkblue") +
#     labs(title = paste0(ident_1, " vs ", ident_2),
#          x = "Log2(Average expression)",
#          y = "Log2FC",
#          color = "") +
#     theme_classic() +
#     scale_color_manual(
#       values = c("darkgrey", "lightgreen"), 
#       labels = c("non-sDE genes",
#                  "sDE genes")) +
#     theme(axis.title = element_text(size = 11),
#           legend.text = element_text(size = 11),
#           plot.title = element_text(hjust = 0.5, face = "bold")) +
#     geom_text_repel(data = deg_subset_sig, aes(label = gene), color = "black", max.overlaps = 20)  
#   
#   # Save ggplot
#   ggsave(filename = paste0("results/DEanalysis/cd45pos_maplot_DEAgeneral_", ident, ".png"),
#          plot = maplot_general)
#   
#   print(maplot_general)
# }


################################################################################
############################### General Functions ############################## 
################################################################################


convert_seurat_to_h5ad <- function(seurat_obj, intermediate_h5_file){
  # Create new Seurat object
  seurat_obj_new <- Seurat::CreateSeuratObject(
    counts = seurat_obj@assays$RNA@counts,
    meta.data = seurat_obj@meta.data
  )
  seurat_obj_new@assays$RNA@data <- seurat_obj@assays$RNA@data
  seurat_obj_new@reductions <- seurat_obj@reductions
  seurat_obj_new@assays$RNA@meta.features <- seurat_obj@assays$RNA@meta.features
  
  # Convert Seurat to Anndata object
  # First, save it as a Seurat h5 object
  SeuratDisk::SaveH5Seurat(seurat_obj_new, filename = intermediate_h5_file, overwrite = TRUE)
  
  # Then, convert the h5 file into an AnnData object
  SeuratDisk::Convert(intermediate_h5_file, dest = "h5ad", overwrite = TRUE)
  
  # Remove Seurat h5 intermediate object
  base::system2('rm', args = intermediate_h5_file)
}


convert_seuratdiet_to_h5ad <- function(seurat_obj, intermediate_h5_file){

  # Convert Seurat to Anndata object
  # First, save it as a Seurat h5 object
  SeuratDisk::SaveH5Seurat(seurat_obj, filename = intermediate_h5_file, overwrite = TRUE)
  
  # Then, convert the h5 file into an AnnData object
  SeuratDisk::Convert(intermediate_h5_file, dest = "h5ad", overwrite = TRUE)
  
  # Remove Seurat h5 intermediate object
  base::system2('rm', args = intermediate_h5_file)
}




convert_genes_human_to_mouse <- function(hsapiens_gene_list){
  ###
  # Converts human genes to mouse genes
  #
  # Parameters
  # hsapiens_gene_list: a list of human genes to be converted
  #
  # Return
  # It returns a list of mouse genes
  ###
  require(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol",
                   values = hsapiens_gene_list,
                   mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows = T)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  row.names(genesV2) <- genesV2$HGNC.symbol
  mmusculus_gene_list <- genesV2[x, 2 ]
  
  return(mmusculus_gene_list)
}


convert_genes_mouse_to_human <- function(mmusculus_gene_list){
  ###
  # Converts mouse genes to human genes
  #
  # Parameters
  # mmusculus_gene_list: a list of mouse genes to be converted
  #
  # Return
  # It returns a list of human genes
  ###
  require(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol",
                   values = mmusculus_gene_list,
                   mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human,
                   uniqueRows = T)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  row.names(genesV2) <- genesV2$HGNC.symbol
  hsapiens_gene_list <- genesV2[x, 2 ]
  
  return(hsapiens_gene_list)
}



# 
#' This function takes in a Seurat 3.0 object and returns a named list with 2
#' objects formated to be loaded in the ShinyApp:
#' https://singlecellgenomics-cnag-crg.shinyapps.io/Annotation/
#' 
#' 1st the metadata + the coordinates of the 2D embeding, 
#' the later with names coord_x and coord_y, and second the expression matrix selected.
#'
#' @param se_obj Object of class Seurat from which we want to extract the information.
#' @param assay Object of class Character indicating from which assay to extract expression data.
#' @param slot Object of class Character indicating from which slot to extract expression data, by default data.
#' @param reduction Object of class Character indicating from which dimensionality reduction we want to extract the coordinates.
#' @return This function returns a named list, the first position contains the joint metadata + 2D embeding and the second contains the expression data matrix. 2D dimensions are labelles as coord_x and coord_y.
#' @export
#' @examples
#'
#'
prepare_se_obj <- function(se_obj,
                           assay,
                           slot = "data",
                           reduction = "umap") {
  suppressMessages(library(dplyr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tibble))
  suppressMessages(library(SummarizedExperiment))
  
  # Input check
  if (! is(se_obj, "Seurat")) stop("Object passed is not a Seurat object;")
  if (! assay %in% Seurat::Assays(se_obj)) stop("assay not in the Seurat object's available assays.")
  if (! slot %in% c("counts", "data", "scale.data")) warning("slot not in the Seurat object's assays.")
  if (! reduction %in% names(se_obj@reductions)) stop("reduction not in the Seurat object's available reductions.")
  
  # Extract 2D cordinates
  embed_df <- se_obj@reductions[[reduction]]@cell.embeddings %>%
    data.frame() %>%
    tibble::rownames_to_column("barcode")
  
  # Change 2D coordinate names to coord_x and coord_y
  colnames(embed_df) <- c("barcode", "coord_x", "coord_y")
  
  # Join metadata with coordinates
  metadata_df <- se_obj@meta.data %>%
    data.frame() %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(embed_df, by = "barcode")
  
  # Extract expression data
  assay_data <- Seurat::GetAssayData(se_obj,
                                     slot = slot,
                                     assay = assay)
  
  return(list(metadata = metadata_df,
              expr_data = assay_data))
}



################################################################################
############################### General Variables ############################## 
################################################################################

