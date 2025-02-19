# boxplot -----------------------------------------------------------------

create_boxplot <- function(df, metadata) {

    fill_colors <- c("CA1" = "#fa7876",
                     "CA3" = "#edd9a3",
                     "DG" = "#cc79a7")


    df_long <- df %>%
        tidyr::pivot_longer(cols = starts_with("S"),
                            values_to = "log2_intensities",
                            names_to = "sample_id") %>%
        dplyr::inner_join(metadata) %>%
        mutate(sample_id = factor(sample_id, levels = unique(sample_id[order(area)]))) %>%
        dplyr::select(sample_id, log2_intensities, intervention, area)

    # Generate the boxplot
    p <- ggplot(df_long, aes(x = sample_id, y = log2_intensities, fill = area)) +
        geom_boxplot() +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
        scale_fill_manual(values = fill_colors) +
        labs(x = "sample_id", y = "log2_intensities") +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8))

    return(p)
}



# eulerr ------------------------------------------------------------------

calculate_venn_counts <- function(df1, df2, df3, name1 = "Set1", name2 = "Set2", name3 = "Set3") {
  # Extract row names from data frames
  set1 <- rownames(df1)
  set2 <- rownames(df2)
  set3 <- rownames(df3)
  
  # Compute intersections
  n_set1_set2 <- length(intersect(set1, set2))
  n_set1_set3 <- length(intersect(set1, set3))
  n_set2_set3 <- length(intersect(set2, set3))
  
  # Compute triple intersection
  n_set1_set2_set3 <- length(intersect(intersect(set1, set2), set3))
  
  # Adjust for double counting
  n_set1_set2_only <- n_set1_set2 - n_set1_set2_set3
  n_set1_set3_only <- n_set1_set3 - n_set1_set2_set3
  n_set2_set3_only <- n_set2_set3 - n_set1_set2_set3
  
  # Compute exclusive elements
  n_set1_only <- length(set1) - (n_set1_set2_only + n_set1_set3_only + n_set1_set2_set3)
  n_set2_only <- length(set2) - (n_set1_set2_only + n_set2_set3_only + n_set1_set2_set3)
  n_set3_only <- length(set3) - (n_set1_set3_only + n_set2_set3_only + n_set1_set2_set3)
  
  # Return named vector
  return(setNames(c(
    n_set1_only, n_set2_only, n_set3_only,
    n_set1_set2_only, n_set1_set3_only, n_set2_set3_only,
    n_set1_set2_set3
  ), c(
    name1, name2, name3,
    paste(name1, name2, sep = "&"), paste(name1, name3, sep = "&"), paste(name2, name3, sep = "&"),
    paste(name1, name2, name3, sep = "&")
  )))
}



# violin plot of area-specific proteins  -----------------------------------------------------------------

plot_hipp_proteins <- function(df, metadata,
                           protein_list = c("Prox1","Dkk3","Calb2", "Pdzd2",
                                            "Ociad2", "Cacng5", "Fibcd1", "Iyd", "Wfs1", "Dcn"),
                           areaspecific_protein = list(
                               pyramidal = c("Ociad2", "Wfs1", "Fibcd1", "Iyd", "Dcn"),
                               mossy = c("Calb2"),
                               non_granule = c("Dkk3"),
                               granule = c("Prox1")),
                           fill_colors = c("CA1" = "#fa7876",
                                         "CA3" = "#edd9a3",
                                         "DG" = "#cc79a7")) {

    #extract df - specific proteins
    df_long <- df[row.names(df) %in% protein_list, ] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("proteins") %>%
        tidyr::pivot_longer(cols = !proteins,
                            values_to = "log2_intensities",
                            names_to = "sample_id") %>%
        #dplyr::inner_join(metadata_m) %>%
        dplyr::inner_join(metadata_f) %>%
        mutate(proteins_by_tissue = as.factor(case_when(
            proteins %in% areaspecific_protein$pyramidal ~ "pyramidal proteins",
            proteins %in% areaspecific_protein$mossy ~ "mossy proteins",
            proteins %in% areaspecific_protein$non_granule ~ "non_granule proteins",
            proteins %in% areaspecific_protein$granule ~ "granule proteins",
            TRUE ~ "Other"
        ))) %>%
        mutate(proteins = factor(proteins, levels = unique(proteins[order(proteins_by_tissue)]))) %>%
        dplyr::select(proteins, sample_id, log2_intensities, area, intervention, proteins_by_tissue)

 
    # plot
    p <- ggplot(df_long, aes(x = area, y = log2_intensities, fill = area)) +
        geom_violin(trim = FALSE) +  # Set trim to FALSE to display the entire violin plot
        geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
        scale_fill_manual(values = fill_colors) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
        labs(x = "Tissue", y = "log2 protein abundance", title = "known proteins and their abundances in our data") +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8)) +
        facet_wrap(proteins~proteins_by_tissue, scales = "free_y", nrow =2)

    list(data = df_long, plot = p)
}



# mds plot ----------------------------------------------------------------

create_mds_plots <- function(mds_data, metadata) {
  metadata_variables <- c("area", "intervention", "batch", "termination_date","cryo_date", "lcm_date",
                          "duration_frozen_tissue", "duration_frozen_section", "duration_frozen_postlcm"
  )
  
  # Define the color palette for the area variable
  area_colors <- c(
    "CA1" = "#fa7876",
    "CA3" = "#edd9a3",
    "DG" = "#cc79a7")
  
  # Define the shapes
  shape_annotations <- list(
    "intervention" = as.factor(c("FFMD" = 21, "LFD" = 22)),
    "batch" = as.factor(c("1" = 23, "3" = 24)),
    "termination_date" = as.factor(c("2024-04-19", "2024-04-25", "2024-04-26")),
    "cryo_date" = as.factor(c("2024-06-12" = 0, "2024-06-13" = 15, "2024-06-14" = 2)),
    "lcm_date" = as.factor(c("2024-06-23" = 0, "2024-06-26" = 1, "2024-06-24" = 2, "2024-07-10" = 3,
                             "2024-07-14" = 4, "2024-07-16" = 5, "2024-07-17" = 6, "2024-07-29" = 7)),
    "duration_frozen_tissue" = as.factor(c("47" = 0, "48" = 1, "50" = 2, "54" = 3, "55" = 4)),
    "duration_frozen_section" = as.factor(c("10" = 0, "14" = 1, "15" = 2, "27" = 3, "31" = 4,
                                            "32" = 5, "34" = 6, "35" = 7, "46" = 8)),
    "duration_frozen_postlcm" = as.factor(c("135" = 0, "147" = 1, "148" = 2, "150" = 3, "154" = 4,
                                            "167" = 5, "168" = 6, "171" = 7)))
  
  
  
  plot_data_area <- mds_data[, c("dim1", "dim2", "area")]
  
  area_plot <- ggplot(plot_data_area,
                      aes(x = dim1, y = dim2, colour = area)) +
    geom_point(size = 4) +
    labs(title = "df_fz_diann__area", colour = "area") +
    theme_minimal() +
    scale_colour_manual(values = area_colors)
  
  # Save the area plot
  ggsave(filename = "doc/df_fz_diann_mds_dim1_dim2.pdf", plot = area_plot, width = 6, height = 6, units = "in", dpi = 300)
  print(area_plot)
  
  # Generate plots for each metadata variable with area as fill and variable as shape
  for (variable in metadata_variables) {
    plot_data <- mds_data[, c("dim1", "dim2", "area", variable)]
    
    plot <- ggplot(plot_data,
                   aes(x = dim1, y = dim2, colour = area, shape = as.factor(mds_data[[variable]]))) +
      geom_point(size = 4) +
      labs(title = paste("df_fz_diann_", variable), colour = "area", shape = variable) +
      theme_minimal() +
      scale_colour_manual(values = area_colors)
    
    # Apply shape mapping for the current variable
    if (variable %in% names(shape_annotations)) {
      plot <- plot + scale_shape_manual(values = shape_annotations[[variable]])
    }
    
    # Construct the file path and save the plot
    filename <- file.path("doc/", paste0("df_fz_diann_mds_dim1_dim2_", variable, ".pdf"))
    ggsave(filename = filename, plot = plot, width = 6, height = 6, units = "in", dpi = 300)
    print(plot)
  }
}

## female 
create_mds_plots_fem <- function(mds_data, metadata) {
  metadata_variables <- c("area", "intervention", "termination_date")
  
  # Define the color palette for the area variable
  area_colors <- c(
    "CA1" = "#fa7876",
    "CA3" = "#edd9a3",
    "DG" = "#cc79a7")
  
  # Define the shapes
  shape_annotations <- list(
    "intervention" = as.factor(c("FFMD" = 21, "LFD" = 22)),
    "termination_date" = as.factor(c("2024-10-19", "2024-10-20")),
    )
  
  
  
  plot_data_area <- mds_data[, c("dim1", "dim2", "area")]
  
  area_plot <- ggplot(plot_data_area,
                      aes(x = dim1, y = dim2, colour = area)) +
    geom_point(size = 4) +
    labs(title = "diann_df_fem_cleaned_fz__area", colour = "area") +
    theme_minimal() +
    scale_colour_manual(values = area_colors)
  
  # Save the area plot
  ggsave(filename = "doc/female/diann_df_fem_cleaned_fz_mds_dim1_dim2.pdf", plot = area_plot, width = 6, height = 6, units = "in", dpi = 300)
  print(area_plot)
  
  # Generate plots for each metadata variable with area as fill and variable as shape
  for (variable in metadata_variables) {
    plot_data <- mds_data[, c("dim1", "dim2", "area", variable)]
    
    plot <- ggplot(plot_data,
                   aes(x = dim1, y = dim2, colour = area, shape = as.factor(mds_data[[variable]]))) +
      geom_point(size = 4) +
      labs(title = paste("diann_df_fem_cleaned_fz", variable), colour = "area", shape = variable) +
      theme_minimal() +
      scale_colour_manual(values = area_colors)
    
    # Apply shape mapping for the current variable
    if (variable %in% names(shape_annotations)) {
      plot <- plot + scale_shape_manual(values = shape_annotations[[variable]])
    }
    
    # Construct the file path and save the plot
    filename <- file.path("doc/female", paste0("diann_df_fem_cleaned_fz_dim1_dim2_", variable, ".pdf"))
    ggsave(filename = filename, plot = plot, width = 6, height = 6, units = "in", dpi = 300)
    print(plot)
  }
}


# mds plot for batchcorr ---------------------------------------------------------------

create_mds_plots_bc <- function(mds_data, metadata) {
  metadata_variables <- c("area", "intervention", "batch", "cryo_date", "lcm_date",
                          "duration_frozen_tissue", "duration_frozen_section", "duration_frozen_postlcm"
  )
  
  # Define the color palette for the area variable
  area_colors <- c(
    "CA1" = "#fa7876",
    "CA3" = "#edd9a3",
    "DG" = "#cc79a7")
  
  # Define the shapes
  shape_annotations <- list(
    "intervention" = as.factor(c("FFMD" = 21, "LFD" = 22)),
    "batch" = as.factor(c("1" = 23, "3" = 24)),
    "cryo_date" = as.factor(c("2024-06-12" = 0, "2024-06-13" = 15, "2024-06-14" = 2)),
    "lcm_date" = as.factor(c("2024-06-23" = 0, "2024-06-26" = 1, "2024-06-24" = 2, "2024-07-10" = 3,
                             "2024-07-14" = 4, "2024-07-16" = 5, "2024-07-17" = 6, "2024-07-29" = 7)),
    "duration_frozen_tissue" = as.factor(c("47" = 0, "48" = 1, "50" = 2, "54" = 3, "55" = 4)),
    "duration_frozen_section" = as.factor(c("10" = 0, "14" = 1, "15" = 2, "27" = 3, "31" = 4,
                                            "32" = 5, "34" = 6, "35" = 7, "46" = 8)),
    "duration_frozen_postlcm" = as.factor(c("135" = 0, "147" = 1, "148" = 2, "150" = 3, "154" = 4,
                                            "167" = 5, "168" = 6, "171" = 7)))
  
  
  
  plot_data_area <- mds_data[, c("dim1", "dim2", "area")]
  
  area_plot <- ggplot(plot_data_area,
                      aes(x = dim1, y = dim2, colour = area)) +
    geom_point(size = 4) +
    labs(title = "df_fz_diann_batchcorr__area", colour = "area") +
    theme_minimal() +
    scale_colour_manual(values = area_colors)
  
  # Save the area plot
  ggsave(filename = "doc/df_fz_mds_dim1_dim2_diann_batchcorr.pdf", plot = area_plot, width = 6, height = 6, units = "in", dpi = 300)
  print(area_plot)
  
  # Generate plots for each metadata variable with area as fill and variable as shape
  for (variable in metadata_variables) {
    plot_data <- mds_data[, c("dim1", "dim2", "area", variable)]
    
    plot <- ggplot(plot_data,
                   aes(x = dim1, y = dim2, colour = area, shape = as.factor(mds_data[[variable]]))) +
      geom_point(size = 4) +
      labs(title = paste("df_fz_bc", variable), colour = "area", shape = variable) +
      theme_minimal() +
      scale_colour_manual(values = area_colors)
    
    # Apply shape mapping for the current variable
    if (variable %in% names(shape_annotations)) {
      plot <- plot + scale_shape_manual(values = shape_annotations[[variable]])
    }
    
    # Construct the file path and save the plot
    filename <- file.path("doc", paste0("df_fz_mds_dim1_dim2_diann_batchcorr", variable, ".pdf"))
    ggsave(filename = filename, plot = plot, width = 6, height = 6, units = "in", dpi = 300)
    print(plot)
  }
}


# mds plot 2 WORK ON HERE--------------------------------------------------------------
create_mds_plots2 <- function(data, metadata, area_label, 
                              dim_x = "dim1", dim_y = "dim2", 
                              color_by = "area", shape_by = "batch", 
                              title = "MDS Plot") {
  # Perform MDS using plotMDS
  mds_data <- plotMDS(data, ndim = 4, plot = FALSE)$cmdscale.out
  
  # Combine MDS results with metadata
  mds_coords <- data.frame(
    mds_data,
    samples = metadata$sample_id,
    area = metadata$area,
    intervention = metadata$intervention,
    batch = metadata$batch
  )
  
  # Rename MDS dimensions
  colnames(mds_coords)[1:4] <- c("dim1", "dim2", "dim3", "dim4")
  
  # Add area label
  mds_coords$area_label <- area_label
  
  # Convert column names to symbols for tidy evaluation
  dim_x_sym <- sym(dim_x)
  dim_y_sym <- sym(dim_y)
  color_by_sym <- sym(color_by)
  shape_by_sym <- sym(shape_by)
  
  # Plot MDS
  plot <- ggplot(mds_coords, aes(x = !!dim_x_sym, y = !!dim_y_sym, color = !!color_by_sym, shape = !!shape_by_sym)) +
    geom_point(size = 4) +
    labs(title = title, x = dim_x, y = dim_y, color = color_by, shape = shape_by) +
    theme_minimal() +
    scale_colour_manual(values = c("CA1" = "#fa7876", "CA3" = "#edd9a3", "DG" = "#cc79a7")) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  return(plot)
}



process_mds <- function(data, metadata, area_label) {
  # Perform MDS
  mds_data <- cmdscale(dist(t(data)), k = 4)
  
  # Combine MDS results with metadata
  mds_coords <- data.frame(
    mds_data,
    samples = metadata$sample_id,
    area = metadata$area,
    intervention = metadata$intervention,
    batch = metadata$batch
  )
  
  # Rename MDS dimensions
  colnames(mds_coords)[1:4] <- c("dim1", "dim2", "dim3", "dim4")
  
  # Add a label for the area if needed
  mds_coords$area_label <- area_label
  
  return(mds_coords)
}


plot_mds <- function(mds_data, dim_x, dim_y, color_by = "area", shape_by = "batch", title = "MDS Plot") {
  ggplot(mds_data, aes_string(x = dim_x, y = dim_y, color = color_by, shape = shape_by)) +
    geom_point(size = 4) +
    labs(title = title, x = dim_x, y = dim_y, color = color_by, shape = shape_by) +
    theme_minimal() +
    scale_colour_manual(values = c("CA1" = "#fa7876", "CA3" = "#edd9a3", "DG" = "#cc79a7")) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
}




create_mds_plots2 <- function(df, metadata, dimensions_list, output_dir = "doc", prefix = "df_mds") {
  mds_data <- plotMDS(df, ndim = max(unlist(dimensions_list)))
  
  mds_data <- data.frame(
    mds_data$cmdscale.out[, seq_len(max(unlist(dimensions_list)))],
    samples = metadata$sample_id,
    area = metadata$area,
    intervention = metadata$intervention,
    batch = metadata$batch
  )
  
  colnames(mds_data)[1:max(unlist(dimensions_list))] <- paste0("dim", 1:max(unlist(dimensions_list)))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (dims in dimensions_list) {
    dim_x <- paste0("dim", dims[1])
    dim_y <- paste0("dim", dims[2])
    
    plot <- ggplot(mds_data, aes_string(x = dim_x, y = dim_y, colour = "area", shape = "batch")) +
      geom_point(size = 4) +
      labs(title = paste(prefix, paste(dim_x, dim_y, sep = "_"), sep = "_"), colour = "area") +
      theme_minimal() +
      scale_colour_manual(values = c("CA1" = "#fa7876",
                                     "CA3" = "#edd9a3",
                                     "DG" = "#cc79a7"))
    
    # Save plot
    output_path <- file.path(output_dir, paste0(prefix, "_", dim_x, "_", dim_y, "_batch.pdf"))
    ggsave(output_path, plot = plot, device = "pdf")
  }
  
  # Additional plots for intervention
  for (dims in dimensions_list) {
    dim_x <- paste0("dim", dims[1])
    dim_y <- paste0("dim", dims[2])
    
    plot <- ggplot(mds_data, aes_string(x = dim_x, y = dim_y, colour = "area", shape = "intervention")) +
      geom_point(size = 4) +
      labs(title = paste(prefix, paste(dim_x, dim_y, "intervention", sep = "_"), sep = "_"), colour = "area") +
      theme_minimal() +
      scale_colour_manual(values = c("CA1" = "#fa7876",
                                     "CA3" = "#edd9a3",
                                     "DG" = "#cc79a7"))
    
    # Save plot
    output_path <- file.path(output_dir, paste0(prefix, "_", dim_x, "_", dim_y, "_intervention.pdf"))
    ggsave(output_path, plot = plot, device = "pdf")
  }
}

# z-score  ----------------------------------------------------------------

z_score_by_median <- function(df) {
  df %>%
    as.data.frame() %>%
    mutate(across(everything(), ~ scale(.x, center = median(.x, na.rm = TRUE), scale = T) %>%
                    as.numeric())) %>% 
    mutate(protein = rownames(.))
}



# valid count - bar plot --------------------------------------------------
barplot_validcounts <- function(df, metadata) {
  # Calculate number of valid (non-NA) values per column
  counts_valid <- colSums(!is.na(df))
  
  # Create a dataframe with valid counts and metadata
  df_valid <- tibble(
    counts_valid = counts_valid,
    sampleid = metadata$sample_id,
    area = metadata$area
  ) %>%
    mutate(
      sampleid = factor(sampleid, levels = unique(sampleid[order(metadata$area)])),
      group = factor(area, levels = unique(area))
    )
  
  # Define cutoff values for each area
  thresholds <- tibble(
    area = c("CA1", "CA3", "DG"),
    #cutoff = c(3982, 4028, 4065) #male
    cutoff = c(3471, 3570, 3997) #female
  )
  
  # Create the bar plot
  ggplot(df_valid, aes(x = sampleid, y = counts_valid, fill = area)) +
    geom_col() +
    scale_fill_manual(values = c("CA1" = "#fa7876", 
                                 "CA3" = "#edd9a3", 
                                 "DG" = "#ca3c97")) +
    scale_x_discrete(labels = df_valid$sampleid) +
    theme_minimal() +
    geom_hline(data = thresholds, aes(yintercept = cutoff, color = area), 
               linetype = "dashed", size = 0.5) +
    scale_color_manual(values = c("CA1" = "black", 
                                  "CA3" = "black", 
                                  "DG" = "black")) +
    labs(x = "Sample", y = "Valid Counts", fill = "area", color = "Cutoff") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    facet_grid(~ area, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = function(x) ifelse(df_valid$area[match(x, df_valid$sampleid)] == unique(df_valid$area), x, ""))
}

