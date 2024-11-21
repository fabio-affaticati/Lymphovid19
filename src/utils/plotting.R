pvalue_to_asterisk <- function(pvalue) {
  
  if (pvalue < 0.001) {
    pvalue <- "***"
  } else if (pvalue < 0.01) {
    pvalue <- "**"
  } else if (pvalue < 0.05) {
    pvalue <- "*"
  } else {
    pvalue <- ""
  }
  return (pvalue)

}



plot_clonal_depth <- function(data, path){

    require(ggplot2)
    require(tidyverse)
    require(scales)
    
    # Define colors for conditions
    colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")

    p <-ggplot(data, aes(x=TIMEPOINTS,y=cloneFraction, fill=CONDITION)) +
    geom_boxplot(notch = FALSE, width = 0.5, outlier.shape = NA) +
    geom_point() +
    annotation_logticks(sides = "l") +
    scale_fill_manual(values = colors) +  # Set fill colors manualsly
    scale_colour_manual(values = colors) +  # Set line colors manually
    geom_line(aes(group=SAMPLE, colour = CONDITION)) +
    theme_bw(base_size = 10) +
    ylim(0.00005, 0.0008) +
    labs(title = "Repertoire clonal fraction of Covid-specific TCRs",
         y = "Clonal fraction",
         x = "Timepoint") +
    theme(
      axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",family = "JetMono Brains"),
    ) + 
    scale_y_continuous(
      labels = scales::percent_format()  # Format y-axis labels as percentages
    )

    ggsave(paste0(path,"scatteratio_depth.png"))
}



plot_scatteratio_breadth <- function(data, path){

  require(ggplot2)
  require(tidyverse)
  require(scales)

  # Define colors for conditions
  colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")
  
  p <-ggplot(data, aes(x=TIMEPOINTS,y=fraction_sequences, fill=CONDITION)) +
    geom_boxplot(notch = FALSE, width = 0.5, outlier.shape = NA) +
    geom_point() +
    scale_y_continuous(
      labels = scales::percent_format()  # Format y-axis labels as percentages
    ) +
    scale_fill_manual(values = colors) +  # Set fill colors manually
    scale_colour_manual(values = colors) +  # Set line colors manually
    geom_line(aes(group=SAMPLE, colour = CONDITION)) +
    theme_bw(base_size = 10) +
    #ylim(0.1, 0.3) +
    labs(title = "",
         y = "Covid Specific fraction (breadth)",
         x = "Timepoint") +
    theme(
      axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",family = "JetMono Brains"),
    )
  ggsave(paste0(path,"scatteratio_breadth_R.png"))

}


plot_betweenstats_breadth <- function(df, path){

  require(ggstatsplot)
  require(ggplot2)
  require(ggrepel)
  require(ggsignif)
  require(tidyverse)
  
  
  colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")

  p1 <- ggstatsplot::ggwithinstats(data = df,
                                   x = TIMEPOINTS,
                                   y = fraction_sequences,
                                   ylab = 'Breadth',
                                   #title = 'Global Comparison of Breadth Across Timepoints',
                                   pairwise.display = 'all',
                                   #pairwise.display = 's',
                                   p.adjust.method = 'BH',
                                   digits = 4L,
                                   outlier.tagging = TRUE,
                                   ggtheme = theme_bw(base_size = 10) +
                                     theme(
                                       text = element_text(size = 12),  # Adjust text size
                                       axis.text = element_text(size = 12),  # Adjust axis text size
                                       axis.title = element_text(size = 10)  # Adjust axis title size
                                     ),
                                   boxplot.args = list(width = 0.3, alpha = 0.2, color = "gray"),
                                   points.args = list(size = 30, alpha = 1),
                                   violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                                   plot.type = "boxplot", type = "nonparametric",
                                   geom_signif_args = list(textsize = 30),
                                   #geom_signif(list(map_signif_level = TRUE)),
                                   centrality.point.args = list(size = 4,color = "#8a0f00"))
  
  p1 <- p1 + geom_line(data = df, aes(group = SAMPLE), colour = "lightgray", alpha = 0.5)

  ggsave(paste0(path,"testing_breadth_global.png"), p1, width=8)
  
  p2 <- ggstatsplot::grouped_ggwithinstats(data = df,
                                   x = TIMEPOINTS,
                                   y = fraction_sequences,
                                   grouping.var = CONDITION,
                                   ylab = 'Breadth',
                                   #title = 'Global Comparison of Breadth Across Timepoints',
                                   pairwise.display = 'all',
                                   #pairwise.display = 's',
                                   p.adjust.method = 'BH',
                                   digits = 4L,
                                   outlier.tagging = TRUE,
                                   ggtheme = theme_bw(base_size = 10) +
                                     theme(
                                       text = element_text(size = 12),  # Adjust text size
                                       axis.text = element_text(size = 12),  # Adjust axis text size
                                       axis.title = element_text(size = 10)  # Adjust axis title size
                                     ),
                                   boxplot.args = list(width = 0.3, alpha = 0.2, color = "gray"),
                                   points.args = list(size = 30, alpha = 1),
                                   violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                                   plot.type = "boxplot", type = "nonparametric",
                                   geom_signif_args = list(textsize = 30),
                                   #geom_signif(list(map_signif_level = TRUE)),
                                   centrality.point.args = list(size = 4,color = "#8a0f00"))

  ggsave(paste0(path,"testing_breadth_grouped.png"), p2, width=8)



  p3 <- ggplot(df, aes(TIMEPOINTS, fraction_sequences, group = interaction(TIMEPOINTS, CONDITION))) +
    geom_point(aes(color = CONDITION, fill = after_scale(alpha(colour, 0.5))), 
               position = position_jitterdodge(dodge.width = 0.9, 0.1),
               size = 2, shape = 21) +
    geom_boxplot(fill = NA, color = "gray", width = 0.2, linewidth = 0.4,
                 position = position_dodge(0.9)) +
    geom_point(stat = "summary", size = 3, color = "#8a0f00",
               position = position_dodge(0.9), fun = median) +
    geom_label_repel(stat = "summary", fun = median, size = 3,
                     aes(label = paste0("hat(mu)*scriptstyle(median)==", 
                                        round(after_stat(y), 5))),
                     parse = TRUE, position = position_dodge(0.9)) +
    geom_signif(y_position = 0.0019, xmin = 1:3 - 0.22, xmax = 1:3 + 0.22,
                annotations = scales::pvalue(sapply(split(df, df$TIMEPOINTS), 
                                                    \(x) wilcox.test(fraction_sequences~CONDITION, x, exact = FALSE)$p.value),
                                             add_p = TRUE)) +
    scale_y_continuous(sec.axis = sec_axis(~., 
                                           bquote(Pairwise~Test~paste(":")~bold(Wilcoxon~Test)))) +
    scale_y_continuous(
      labels = scales::percent_format()  # Format y-axis labels as percentages
    ) +
    #scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(axis.title = element_text(face = 2),
          legend.position = "bottom",
          axis.text.y.right = element_blank()) +
    labs(y = "Breadth",
         x = "Timepoint")

  ggsave(paste0(path,"testing_breadth_betweenconditions.png"), p3, width=8)



  p <- combine_plots(
    plotlist = list(p1, p2, p3),
    #plotgrid.args = list(nrow = 2),
    annotation.args = list(
      title = "",
      caption = ""
    )
  )


  ggsave(paste0(path,"composite_breadth.png"), p, width=30)

}



plot_betweenstats_clone_fraction <- function(df, path){
  
  require(ggstatsplot)
  require(ggplot2)
  require(ggrepel)
  require(ggsignif)
  require(tidyverse)
  require(rstatix)
  require(ggpubr)
  
  colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")

  p1 <- ggstatsplot::ggwithinstats(data = df,
                                   x = TIMEPOINTS,
                                   y = cloneFraction,
                                   ylab = 'Clone Fraction',
                                   #title = 'Global Comparison of Breadth Across Timepoints',
                                   pairwise.display = 'all',
                                   #pairwise.display = 's',
                                   p.adjust.method = 'BH',
                                   digits = 4L,
                                   outlier.tagging = TRUE,
                                   ggtheme = theme_bw(base_size = 10) +
                                     theme(
                                       text = element_text(size = 12),  # Adjust text size
                                       axis.text = element_text(size = 12),  # Adjust axis text size
                                       axis.title = element_text(size = 10)  # Adjust axis title size
                                     ),
                                   #ggplot.component = list(scale_y_log10()),
                                   boxplot.args = list(width = 0.3, alpha = 0.2, color = "gray"),
                                   points.args = list(size = 30, alpha = 1),
                                   violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                                   plot.type = "boxplot", type = "nonparametric",
                                   geom_signif_args = list(textsize = 30,
                                                    y_position = c(0.023, 0.024, 0.025),),
                                   #geom_signif(list(map_signif_level = TRUE)),
                                   centrality.point.args = list(size = 4,color = "#8a0f00"))
  p1 <- p1 + geom_line(data = df, aes(group = SAMPLE), colour = "lightgray", alpha = 0.5)

  ggsave(paste0(path,"testing_clonefraction_global.png"), p1, width=8)


  
  p2 <- ggstatsplot::grouped_ggwithinstats(data = df,
                                   x = TIMEPOINTS,
                                   y = cloneFraction,
                                   grouping.var = CONDITION,
                                   ylab = 'Clone Fraction',
                                   #title = 'Global Comparison of Breadth Across Timepoints',
                                   pairwise.display = 'all',
                                   #pairwise.display = 's',
                                   p.adjust.method = 'BH',
                                   digits = 4L,
                                   outlier.tagging = TRUE,
                                   ggtheme = theme_bw(base_size = 10) +
                                     theme(
                                       text = element_text(size = 12),  # Adjust text size
                                       axis.text = element_text(size = 12),  # Adjust axis text size
                                       axis.title = element_text(size = 10)  # Adjust axis title size
                                     ),
                                   #ggplot.component = list(scale_y_log10()),
                                   boxplot.args = list(width = 0.3, alpha = 0.2, color = "gray"),
                                   points.args = list(size = 30, alpha = 1),
                                   violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                                   plot.type = "boxplot", type = "nonparametric",
                                   geom_signif_args = list(textsize = 30),
                                   #geom_signif(list(map_signif_level = TRUE)),
                                   centrality.point.args = list(size = 4,color = "#8a0f00"))

  ggsave(paste0(path,"testing_clonefraction_grouped.png"), p2, width=8)



  p3 <- ggplot(df, aes(TIMEPOINTS, cloneFraction, group = interaction(TIMEPOINTS, CONDITION))) +
    geom_point(aes(color = CONDITION, fill = after_scale(alpha(colour, 0.5))), 
               position = position_jitterdodge(dodge.width = 0.9, 0.1),
               size = 2, shape = 21) +
    geom_boxplot(fill = NA, color = "gray", width = 0.2, linewidth = 0.4,
                 position = position_dodge(0.9)) +
    geom_point(stat = "summary", size = 3, color = "#8a0f00",
               position = position_dodge(0.9), fun = median) +
    #geom_label_repel(stat = "summary", fun = median, size = 3,
    #                 aes(label = paste0("hat(mu)*scriptstyle(median)==", 
    #                                    round(after_stat(y), 5))),
    #                 parse = TRUE, position = position_dodge(0.9)) +
    geom_signif(y_position = 0.0008, xmin = 1:3 - 0.22, xmax = 1:3 + 0.22,
                annotations = scales::pvalue(sapply(split(df, df$TIMEPOINTS), 
                                                    \(x) wilcox.test(cloneFraction~CONDITION, x, exact = FALSE)$p.value),
                                             add_p = TRUE)) +
    scale_y_continuous(sec.axis = sec_axis(~., 
                                           bquote(Pairwise~Test~paste(":")~bold(Wilcoxon~Test)))) +
    scale_y_continuous(
      labels = scales::percent_format()  # Format y-axis labels as percentages
    ) +
    #scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(axis.title = element_text(face = 2),
          legend.position = "bottom",
          axis.text.y.right = element_blank()) +
    labs(y = "Clone Fraction",
         x = "Timepoint") #+ scale_y_continuous(trans='log10')

  ggsave(paste0(path,"testing_clonefraction_betweenconditions.png"), p3, width=8)


  p <- combine_plots(
    plotlist = list(p1, p2, p3),
    #plotgrid.args = list(nrow = 2),
    annotation.args = list(
      title = "",
      caption = ""
    )
  )


  ggsave(paste0(path,"composite_clonefraction.png"), p, width=30)

}



linear_mixed_model <- function(df, model_path){
  
  require(lme4)
  require(lmerTest)
  library(emmeans)
  require(ggplot2)


  df$SAMPLE <- as.factor(df$SAMPLE)
  df$TIMEPOINTS <- as.factor(df$TIMEPOINTS)
  df$CONDITION <- as.factor(df$CONDITION)
  #df$TIMEPOINTS <- relevel(df$TIMEPOINTS, ref = "V3")
  #df$CONDITION <- relevel(df$CONDITION, ref = "Lymphomas")

  # if fraction_sequences is in the columns of the dataframe
  if("fraction_sequences" %in% colnames(df)){
    model <- lmer(fraction_sequences ~ CONDITION * TIMEPOINTS + (1 | SAMPLE), data = df)
    saveRDS(model, file = paste0(model_path, "lmm_model_breadth.rds"))
    emms <- emmeans(model, specs = ~ TIMEPOINTS | CONDITION, type = "response", adjust = "tukey")
    print(contrast(emms, method = "pairwise"))
    ggsave(paste0(model_path, "lmm_model_breadth_emms.png"), plot(emms, comparisons = TRUE), width = 10)
    ggsave(paste0(model_path, "lmm_model_breadth_pwpp.png"), pwpp(emms, method = "pairwise",
            sort = F, values = T, type = "response") + 
            geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.5), width = 10)
  }

  # if cloneFraction is in the columns of the dataframe
  else if("cloneFraction" %in% colnames(df)){
    model <- lmer(cloneFraction ~ CONDITION * TIMEPOINTS + (1 | SAMPLE), data = df)
    saveRDS(model, file = paste0(model_path, "lmm_model_depth.rds"))
    emms <- emmeans(model, specs = ~ TIMEPOINTS | CONDITION, type = "response", adjust = "tukey")
    
    print(contrast(emms, method = "pairwise"))

    ggsave(paste0(model_path, "lmm_model_depth_emms.png"), plot(emms, comparisons = TRUE), width = 10)
    ggsave(paste0(model_path, "lmm_model_depth_pwpp.png"), pwpp(emms, method = "pairwise",
            sort = F, values = T, type = "response") + 
            geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.5), width = 10)
  }
  else{
    print("No suitable column found in the dataframe")
  }
  
  #p <- emmip(model,~ TIMEPOINTS | CONDITION, CIs = TRUE, type = "response")
  #ggsave(paste0(model_path, "test.png"), p, width = 10)

  
}




export_lmm_results <- function(model_path){

  require(jtools)

  # read RDS files
  model_breadth <- readRDS(paste0(model_path, 'lmm_model_breadth.rds'))
  model_depth <- readRDS(paste0(model_path, 'lmm_model_depth.rds'))

  export_summs(model_breadth,  model_depth, to.file = 'pdf', file.name = paste0(model_path, "lmm_model_reports.pdf"),
              model.names = c("LMM_Breadth", "LMM_Depth"), digits = 6, to.numeric = TRUE, error_format = "[{conf.low}, {conf.high}]")

  
  # Rename coefficients
coef_names <- c(
  'Lymphoma' = 'CONDITIONLymphomas',  
  'Time V1' = 'TIMEPOINTSV1',          
  'Time V3' = 'TIMEPOINTSV3',          
  'Lymphoma at V1' = 'CONDITIONLymphomas:TIMEPOINTSV1',  
  'Lymphoma at V3' = 'CONDITIONLymphomas:TIMEPOINTSV3' 
)

  p <- plot_summs(model_breadth, model_depth, ci_level = .95, model.names = c("LMM_Breadth", "LMM_Depth"), coefs = coef_names)

  summary_model_breadth <- summary(model_breadth)$coefficients[, "Pr(>|t|)"]
  summary_model_depth <- summary(model_depth)$coefficients[, "Pr(>|t|)"]

  # annotate the plot with the p values taken from the models
  p <- p + annotate("text", x = 0.00055, y = 5.15, label = pvalue_to_asterisk(summary_model_breadth['CONDITIONLymphomas']), size = 5) +
          annotate("text", x = 0.00055, y = 4.9, label = pvalue_to_asterisk(summary_model_depth['CONDITIONLymphomas']), size = 5) +
          annotate("text", x = 0.00055, y = 4.15, label = pvalue_to_asterisk(summary_model_breadth['TIMEPOINTSV1']), size = 5) +
          annotate("text", x = 0.00055, y = 3.9, label = pvalue_to_asterisk(summary_model_depth['TIMEPOINTSV1']), size = 5) +
          annotate("text", x = 0.00055, y = 3.15, label = pvalue_to_asterisk(summary_model_breadth['TIMEPOINTSV3']), size = 5) +
          annotate("text", x = 0.00055, y = 2.9, label = pvalue_to_asterisk(summary_model_depth['TIMEPOINTSV3']), size = 5) +
          annotate("text", x = 0.00055, y = 2.15, label = pvalue_to_asterisk(summary_model_breadth['CONDITIONLymphomas:TIMEPOINTSV1']), size = 5) +
          annotate("text", x = 0.00055, y = 1.9, label = pvalue_to_asterisk(summary_model_depth['CONDITIONLymphomas:TIMEPOINTSV1']), size = 5) +
          annotate("text", x = 0.00055, y = 1.15, label = pvalue_to_asterisk(summary_model_breadth['CONDITIONLymphomas:TIMEPOINTSV3']), size = 5) +
          annotate("text", x = 0.00055, y = 0.9, label = pvalue_to_asterisk(summary_model_depth['CONDITIONLymphomas:TIMEPOINTSV3']), size = 5) +
    scale_x_continuous(
      labels = scales::percent_format()  # Format y-axis labels as percentages
    )

  ggsave(paste0(model_path, "models_summs.png"), p, width = 10)
  
}


mean_clone_dist <- function(df, selected_resampling_val, plotPath) {

  require(ggplot2)

  bar_plot <- ggplot(df, aes(x = Value, y = Frequency)) +
    geom_point() +
    geom_smooth(stat = "identity") +   # Use 'identity' to plot actual values
    scale_fill_manual(values = c("Below" = "gray", "Above" = "red")) +
    labs(title = "", x = "Mean Covid-specific Clone fraction", y = "Frequency") +
    theme_bw() +
    geom_vline(xintercept = selected_resampling_val, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
    annotate("text", x = selected_resampling_val + 0.000005, y = 22, label = "Selected Resampling",
              angle = 0, vjust = 1.5, color = "red", size = 4) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 14),  # Increase x-axis title size
      axis.title.y = element_text(size = 14),  # Increase y-axis title size
      axis.text.x = element_text(size = 14),   # Increase x-axis text size
      axis.text.y = element_text(size = 12),    # Increase y-axis text size
      panel.grid.major = element_blank(),      # Remove major gridlines
      panel.grid.minor = element_blank(),      # Remove minor gridlines
    )

  ggsave(paste0(plotPath,'dist_permutations.png'), bar_plot, scale = 1, width = 10, height = 8)
}


correlation_heatmap <- function(data, plotsDir){

  require(pheatmap)
  require(ggplot2)
  require(glue)
  require(Hmisc)
  require(RColorBrewer)

  clusters_n = 8

  # Get numeric data for correlation
  numeric_data <- data[, sapply(data, is.numeric)]

  # Create an annotation dataframe using the CONDITION and TIMEPOINTS columns
  annotation_data <- data[, c("CONDITION", "TIMEPOINTS", "SAMPLE", "HLA_A", "HLA_B", "HLA_C")]
  data$unique_id <- paste(data$SAMPLE, data$CONDITION, data$TIMEPOINTS, sep="_")

  # Compute correlation matrix and p-value matrix using Hmisc's rcorr function
  corr_results <- rcorr(as.matrix(t(numeric_data)), type = "pearson")
  mat <- corr_results$r       # Correlation matrix
  pvals <- corr_results$P     # P-value matrix
  
  unique_samples <- unique(data$SAMPLE)
  sample_colors <- hcl.colors(length(unique_samples), palette = "Dark3")
  names(sample_colors) <- unique_samples  # Name colors by sample

  unique_hlasa <- unique(data$HLA_A)
  hlaa_colors <- hcl.colors(length(unique_hlasa), palette = "Grays")
  names(hlaa_colors) <- unique_hlasa

  unique_hlasb <- unique(data$HLA_B)
  hlab_colors <- hcl.colors(length(unique_hlasb), palette = "Grays")
  names(hlab_colors) <- unique_hlasb

  unique_hlasc <- unique(data$HLA_C)
  hlac_colors <- hcl.colors(length(unique_hlasc), palette = "Grays")
  names(hlac_colors) <- unique_hlasc




  # Create a matrix of asterisks based on p-value thresholds
  signif_matrix <- ifelse(pvals < 0.0001, "****",
                    ifelse(pvals < 0.001, "***",
                    ifelse(pvals < 0.01, "**",
                    ifelse(pvals < 0.05, "*", ""))))

  # Set diagonal of the significance matrix to empty strings
  diag(signif_matrix) <- ""

  # Define colors for the annotation bars (you can customize these)
  annotation_colors <- list(
    CONDITION = c("Healthy" = "blue", "Lymphomas" = "red"),  # Replace with your actual conditions
    TIMEPOINTS = c("baseline" = "yellow", "V1" = "orange", "V3" = "purple"),  # Replace with your actual timepoints
    HLA_A = hlaa_colors,
    HLA_B = hlab_colors,
    HLA_C = hlac_colors,
    SAMPLE = sample_colors  # Add sample colors to the annotation colors
  )

  rownames(mat) <- data$unique_id  # Set row names to column names
  colnames(mat) <- data$unique_id  # Set column names to column names
  rownames(annotation_data) <- data$unique_id  # Set unique identifiers as row names for annotation

# Create heatmap with pheatmap
  pheatmap(as.matrix(mat),
    color = colorRampPalette(c("white", "red"))(100),
    annotation_row = annotation_data,   # Add annotations for rows
    annotation_col = annotation_data,   # Add annotations for columns
    annotation_colors = annotation_colors,  # Add colors for the annotations
    main = "",  # Title
    display_numbers = signif_matrix, #TRUE,  # Show correlation values inside the cells
    number_format = "%.2f",  # Format numbers to 2 decimal places
    fontsize = 10,  # Font size for general text
    fontsize_number = 10,  # Font size for numbers in cells
    cutree_rows = clusters_n,  # Number of clusters for rows
    cutree_cols = clusters_n,  # Number of clusters for columns
    cluster_rows = F,  # Cluster rows
    cluster_cols = TRUE,  # Cluster columns
    filename = paste0(plotsDir, glue("vgene_corrplot_euclidean_distance.png")),
    width = 20, height = 20)

  return

}


correlation_heatmap_pertime <- function(data, time, plotsDir){

  clusters_n = 4

  require(pheatmap)
  require(ggplot2)
  require(glue)
  require(Hmisc)
  require(RColorBrewer)

  # Generate correlation matrix using only numeric columns
  # Previous use, now to be removed
  mat <- cor(t(data[, sapply(data, is.numeric)]))

  # Get numeric data for correlation (transpose to   correlate rows)
  numeric_data <- data[, sapply(data, is.numeric)]

  # Create an annotation dataframe using the CONDITION and TIMEPOINTS columns
  annotation_data <- data[, c("CONDITION", "SAMPLE", "HLA_A", "HLA_B", "HLA_C")]
  data$unique_id <- paste(data$SAMPLE, data$CONDITION, sep="_")

  # Compute correlation matrix and p-value matrix using Hmisc's rcorr function
  corr_results <- rcorr(as.matrix(t(numeric_data)), type = "pearson")
  mat <- corr_results$r       # Correlation matrix
  pvals <- corr_results$P     # P-value matrix

  unique_samples <- unique(data$SAMPLE)
  sample_colors <- hcl.colors(length(unique_samples), palette = "Dark3")
  names(sample_colors) <- unique_samples  # Name colors by sample

  unique_hlasa <- unique(data$HLA_A)
  hlaa_colors <- hcl.colors(length(unique_hlasa), palette = "Grays")
  names(hlaa_colors) <- unique_hlasa

  unique_hlasb <- unique(data$HLA_B)
  hlab_colors <- hcl.colors(length(unique_hlasb), palette = "Grays")
  names(hlab_colors) <- unique_hlasb

  unique_hlasc <- unique(data$HLA_C)
  hlac_colors <- hcl.colors(length(unique_hlasc), palette = "Grays")
  names(hlac_colors) <- unique_hlasc


  # Create a matrix of asterisks based on p-value thresholds
  signif_matrix <- ifelse(pvals < 0.0001, "****",
                    ifelse(pvals < 0.001, "***",
                    ifelse(pvals < 0.01, "**",
                    ifelse(pvals < 0.05, "*", ""))))

  # Set diagonal of the significance matrix to empty strings
  diag(signif_matrix) <- ""

  # Define colors for the annotation bars (you can customize these)
  annotation_colors <- list(
    CONDITION = c("Healthy" = "blue", "Lymphomas" = "red"),  # Replace with your actual conditions
    HLA_A = hlaa_colors,
    HLA_B = hlab_colors,
    HLA_C = hlac_colors,
    SAMPLE = sample_colors  # Add sample colors to the annotation colors
  )

  rownames(mat) <- data$unique_id  # Set row names to column names
  colnames(mat) <- data$unique_id  # Set column names to column names
  rownames(annotation_data) <- data$unique_id  # Set unique identifiers as row names for annotation

# Create heatmap with pheatmap
  pheatmap_result <- pheatmap(as.matrix(mat),
    color = colorRampPalette(c("white", "red"))(100),
    annotation_row = annotation_data,   # Add annotations for rows
    annotation_col = annotation_data,   # Add annotations for columns
    annotation_colors = annotation_colors,  # Add colors for the annotations
    main = glue("Correlation plot {time}"),  # Title
    display_numbers = signif_matrix, #TRUE,  # Show correlation values inside the cells
    number_format = "%.2f",  # Format numbers to 2 decimal places
    fontsize = 15,  # Font size for general text
    fontsize_number = 10,  # Font size for numbers in cells
    cutree_rows = clusters_n,  # Number of clusters for rows
    cutree_cols = clusters_n,  # Number of clusters for columns
    cluster_rows = T,  # Cluster rows
    cluster_cols = T,  # Cluster columns
    filename = paste0(plotsDir, glue("vgene_corrplot_euclidean_distance_{time}.png")),
    width = 30, height = 30)



  row_clusters <- cutree(pheatmap_result$tree_row, k = clusters_n)  # Change 'k' to the number of desired clusters
  # Create a data frame with sample names and their corresponding clusters
  cluster_data <- data.frame(SAMPLE = names(row_clusters), CLUSTER = row_clusters)

  metadata <- data.frame(SAMPLE = data$unique_id, CONDITION = data$CONDITION)

  # Merge with metadata to include conditions
  combined_data <- merge(cluster_data, metadata, by = "SAMPLE")

  print(glue("Cluster Analysis for Timepoint {time}:"))
  # Create a contingency table for each cluster and perform the fisher test
  p_values <- c()
  for (cluster in unique(combined_data$CLUSTER)) {
    contingency_table <- table(combined_data$CLUSTER == cluster, combined_data$CONDITION)
    
    print(glue("Contingency Table for Cluster {cluster}:"))
    print(contingency_table)
    # Perform the fisher test
    fisher_test <- fisher.test(contingency_table)
    
    # Collect the p-value
    p_values <- c(p_values, fisher_test$p.value)
  }

  # Print the collected p-values
  print("P-Values from Fisher Tests:")
  print(p_values)

  p_adjusted <- p.adjust(p_values, method = "fdr")

  # Print the adjusted p-values
  print("Adjusted P-Values:")
  print(p_adjusted)
  return


}


pca_biplot_vgenes <- function(data, metadata, plotsDir){

  metadata$CovidBreadth <- metadata$fraction_sequences

  # Load necessary libraries
  require(ggplot2)
  require(gridExtra)
  require(dplyr)
  require(broom)

  # Keep only columns with non-zero variance
  data <- data[, apply(data, 2, var) != 0]  

  # Ensure metadata factors are correctly formatted
  metadata$TIMEPOINTS <- as.factor(metadata$TIMEPOINTS)
  metadata$CONDITION <- as.factor(metadata$CONDITION)
  metadata$SAMPLE <- as.factor(metadata$SAMPLE)
  metadata$HLA_A <- as.factor(metadata$HLA_A)
  metadata$HLA_B <- as.factor(metadata$HLA_B)
  metadata$HLA_C <- as.factor(metadata$HLA_C)
  rownames(metadata) <- rownames(metadata$sample_id)

  # Perform PCA
  pca_result <- prcomp(data, center = TRUE, scale. = TRUE)

  # Extract PCA loadings
  loadings <- as.data.frame(pca_result$rotation)
  loadings$variable <- rownames(loadings)

  # Calculate explained variance
  explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100


  loading_df <- data.frame(
    variable = rownames(loadings),
    PC1 = loadings[, 1],
    PC2 = loadings[, 2]
  )

  # Define quadrants based on PC1 and PC2 signs
  loading_df$quadrant <- with(loading_df, ifelse(PC1 > 0 & PC2 > 0, "Q1",
                                    ifelse(PC1 < 0 & PC2 > 0, "Q2",
                                    ifelse(PC1 < 0 & PC2 < 0, "Q3", "Q4"))))

  # Select the top 5 variables for each quadrant
  top_features <- loading_df %>%
    group_by(quadrant) %>%
    top_n(5, wt = abs(PC1) + abs(PC2))  # Sum of absolute loadings for PC1 and PC2
  
  # Filter loadings to keep only top features
  top_loadings <- loadings %>%
    filter(variable %in% top_features$variable)


  # Prepare data for ggplot
  pca_scores <- as.data.frame(pca_result$x)  # Scores for individuals
  pca_scores$CONDITION <- metadata$CONDITION
  pca_scores$TIMEPOINTS <- metadata$TIMEPOINTS
  pca_scores$SAMPLE <- metadata$SAMPLE
  pca_scores$CovidBreadth <- metadata$CovidBreadth

  # Correlations with PCs (for annotations)
  cor_with_pc1 <- cor(pca_scores$CovidBreadth, pca_scores$PC1, method = "spearman", use = "complete.obs")
  cor_with_pc2 <- cor(pca_scores$CovidBreadth, pca_scores$PC2, method = "spearman", use = "complete.obs")

  # Find common limits for PC1 and PC2
  x_limits <- range(pca_scores$PC1, na.rm = TRUE)
  y_limits <- range(pca_scores$PC2, na.rm = TRUE)

  # PCA biplot
  pca_plot <- ggplot() +
    geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = CONDITION, shape = TIMEPOINTS), size = 2.5) +
    geom_text_repel(data = pca_scores, aes(x = PC1, y = PC2, color = CONDITION, label = SAMPLE), nudge_x = 0.3, size = 3, segment.color = NA) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    # Plot top feature loadings as arrows
    geom_segment(data = top_loadings, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10), 
                arrow = arrow(length = unit(0.2, "cm")), color = "darkgray", alpha = 0.6) +
    # Annotate the features
    geom_text(data = top_loadings, aes(x = PC1 * 10, y = PC2 * 10, label = variable), size = 2.5, hjust = 1.2) +
    # Customize axis labels
    labs(x = paste0("PC1 (", round(explained_var[1], 2), "%)"), 
        y = paste0("PC2 (", round(explained_var[2], 2), "%)")) +
    # Customize colors
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(text = element_text(size = 10),
        legend.position = "none") +  # Remove legend here
    coord_cartesian(xlim = x_limits, ylim = y_limits)

  # Scatter plot of the variable against PC1
  scatter_pc1 <- ggplot(pca_scores, aes(x = PC1, y = CovidBreadth, color = CONDITION, shape = TIMEPOINTS)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(y = "CovidBreadth", x = '', subtitle = paste("Spearman correlation with PC1:", round(cor_with_pc1, 5))) +
    coord_cartesian(xlim = x_limits)

  # Scatter plot of the variable against PC2
  scatter_pc2 <- ggplot(pca_scores, aes(x = CovidBreadth, y = PC2, color = CONDITION, shape = TIMEPOINTS)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(x = "CovidBreadth", y = '', subtitle = paste("Spearman correlation with PC2:", round(cor_with_pc2, 5))) +
    coord_cartesian(ylim = y_limits)

  # Blank plot
  blankPlot <- ggplot() + geom_blank() +
    theme_void()

  # Extract the legend from the PCA plot
  legend_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = CONDITION, shape = TIMEPOINTS)) +
    geom_point() +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(legend.position = "right")  # Keep the legend here

  # Extract the legend as a separate object
  legend <- get_legend(legend_plot)


  combined_plot <- ggarrange(
    scatter_pc1, blankPlot, pca_plot, scatter_pc2, 
    ncol = 2, nrow = 2, widths = c(3, 1.4), heights = c(1.4, 3),
    legend.grob = legend,  # Add the legend as a grob
    legend = "right"       # Position the legend to the right
  )


  # Save the plot
  ggsave(paste0(plotsDir, "pca_biplot_12.png"), plot = combined_plot, width = 15, height = 10, dpi = 600)



 loading_df <- data.frame(
    variable = rownames(loadings),
    PC3 = loadings[, 3],
    PC4 = loadings[, 4]
  )

  # Define quadrants based on PC1 and PC2 signs
  loading_df$quadrant <- with(loading_df, ifelse(PC3 > 0 & PC4 > 0, "Q1",
                                    ifelse(PC3 < 0 & PC4 > 0, "Q2",
                                    ifelse(PC3 < 0 & PC4 < 0, "Q3", "Q4"))))

  # Select the top 5 variables for each quadrant
  top_features <- loading_df %>%
    group_by(quadrant) %>%
    top_n(5, wt = abs(PC3) + abs(PC4))  # Sum of absolute loadings for PC3 and PC4
  
  # Filter loadings to keep only top features
  top_loadings <- loadings %>%
    filter(variable %in% top_features$variable)


  # Correlations with PCs (for annotations)
  cor_with_pc1 <- cor(pca_scores$CovidBreadth, pca_scores$PC3, method = "spearman", use = "complete.obs")
  cor_with_pc2 <- cor(pca_scores$CovidBreadth, pca_scores$PC4, method = "spearman", use = "complete.obs")

  x_limits <- range(pca_scores$PC3, na.rm = TRUE)
  y_limits <- range(pca_scores$PC4, na.rm = TRUE)

  # Create a ggplot biplot
  pca_plot <- ggplot() +
    # Plot individuals
    geom_point(data = pca_scores, aes(x = PC3, y = PC4, color = CONDITION, shape = TIMEPOINTS), size = 2.5) +
    geom_text_repel(data = pca_scores, aes(x = PC3, y = PC4, color = CONDITION, label = SAMPLE), nudge_x = 0.3, size = 3, segment.color = NA) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    # Plot top feature loadings as arrows
    geom_segment(data = top_loadings, aes(x = 0, y = 0, xend = PC3 * 10, yend = PC4 * 10), 
                arrow = arrow(length = unit(0.2, "cm")), color = "darkgray", alpha = 0.6) +
    # Annotate the features
    geom_text(data = top_loadings, aes(x = PC3 * 10, y = PC4 * 10, label = variable), size = 2.5, hjust = 1.2) +
    # Customize axis labels
    labs(x = paste0("PC3 (", round(explained_var[3], 2), "%)"), 
        y = paste0("PC4 (", round(explained_var[4], 2), "%)")) +
    # Customize colors
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(text = element_text(size = 8),
        legend.position = "none") +  # Remove legend here
    coord_cartesian(xlim = x_limits, ylim = y_limits)

    # Scatter plot of the variable against PC3
  scatter_pc3 <- ggplot(pca_scores, aes(x = PC3, y = CovidBreadth, color = CONDITION, shape = TIMEPOINTS)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(y = "CovidBreadth", x = '', subtitle = paste("Spearman correlation with PC3:", round(cor_with_pc1, 5))) +
    coord_cartesian(xlim = x_limits)

  # Scatter plot of the variable against PC2
  scatter_pc4 <- ggplot(pca_scores, aes(x = CovidBreadth, y = PC4, color = CONDITION, shape = TIMEPOINTS)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(x = "CovidBreadth", y = '', subtitle = paste("Spearman correlation with PC4:", round(cor_with_pc2, 5))) +
    coord_cartesian(ylim = y_limits)

  # Blank plot
  blankPlot <- ggplot() + geom_blank() +
    theme_void()

  # Extract the legend from the PCA plot
  legend_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = CONDITION, shape = TIMEPOINTS)) +
    geom_point() +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(legend.position = "right")  # Keep the legend here

  # Extract the legend as a separate object
  legend <- get_legend(legend_plot)


  combined_plot <- ggarrange(
    scatter_pc3, blankPlot, pca_plot, scatter_pc4, 
    ncol = 2, nrow = 2, widths = c(3, 1.4), heights = c(1.4, 3),
    legend.grob = legend,  # Add the legend as a grob
    legend = "right"       # Position the legend to the right
  )

  # Save the plot
  ggsave(paste0(plotsDir, "pca_biplot_34.png"), plot = combined_plot, width = 15, height = 10, dpi = 600)








  library(variancePartition)
  require(tidyr)

  pca_scores <- as.data.frame(pca_result$x)
  rownames(pca_scores) <- rownames(metadata)

  # Step 2: Print the dimensions of variance_input and metadata for verification
  print("Dimensions of variance_input:")
  print(dim(data))
  
  print("Dimensions of metadata:")
  print(dim(metadata))
  
  # Step 3: Check and print summary of metadata
  print("Metadata summary:")
  print(summary(metadata))


  form <- ~ (1 | TIMEPOINTS) + (1 | CONDITION) + (1 | HLA_A) + (1 | HLA_B) + (1 | HLA_C) + (1 | SAMPLE)



 # Step 5: Fit the model, check for errors
  print("Fitting the model...")
  varPart <- tryCatch({
    fitExtractVarPartModel(t(pca_scores), form, metadata)
  }, error = function(e) {
    print(paste("Error in model fitting:", e))
    return(NULL)
  })

  if (is.null(varPart)) {
    stop("Model fitting failed. Check the input data and model formula.")
  }

  # Step 6: Print and check variance partitioning results
  print("Variance partitioning results:")
  print(summary(varPart))
  
  # Step 7: Sort the variance results
  vp <- sortCols(varPart)
  print("Sorted variance partitioning results:")
  print(vp)

   # Step 8: Plot percentage variance explained
  png(paste0(plotsDir, 'variance_explained.png'), width = 800, height = 500, res = 100)
  p <- plotPercentBars(vp[1:5, ])
  print(p)
  dev.off()

  # Step 9: Plot global variance explained
  png(paste0(plotsDir, 'global_variance_explained.png'), width = 800, height = 800, res = 100)
  p <- plotVarPart(vp)
  print(p)
  dev.off()


}


pca_biplot_vgenes_pertime <- function(data, metadata, time, plotsDir){

  metadata$CovidBreadth <- metadata$fraction_sequences
  
  # Load necessary libraries
  require(ggplot2)
  require(gridExtra)
  require(dplyr)
  require(broom)

  # Keep only columns with non-zero variance
  data <- data[, apply(data, 2, var) != 0]  

  # Ensure metadata factors are correctly formatted
  metadata$CONDITION <- as.factor(metadata$CONDITION)
  metadata$SAMPLE <- as.factor(metadata$SAMPLE)
  metadata$HLA_A <- as.factor(metadata$HLA_A)
  metadata$HLA_B <- as.factor(metadata$HLA_B)
  metadata$HLA_C <- as.factor(metadata$HLA_C)
  rownames(metadata) <- rownames(metadata$sample_id)

  # Perform PCA
  pca_result <- prcomp(data, center = TRUE, scale. = TRUE)

  # Extract PCA loadings
  loadings <- as.data.frame(pca_result$rotation)
  loadings$variable <- rownames(loadings)

  # Calculate explained variance
  explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100


  loading_df <- data.frame(
    variable = rownames(loadings),
    PC1 = loadings[, 1],
    PC2 = loadings[, 2]
  )

  # Define quadrants based on PC1 and PC2 signs
  loading_df$quadrant <- with(loading_df, ifelse(PC1 > 0 & PC2 > 0, "Q1",
                                    ifelse(PC1 < 0 & PC2 > 0, "Q2",
                                    ifelse(PC1 < 0 & PC2 < 0, "Q3", "Q4"))))

  # Select the top 5 variables for each quadrant
  top_features <- loading_df %>%
    group_by(quadrant) %>%
    top_n(5, wt = abs(PC1) + abs(PC2))  # Sum of absolute loadings for PC1 and PC2
  
  # Filter loadings to keep only top features
  top_loadings <- loadings %>%
    filter(variable %in% top_features$variable)


  # Prepare data for ggplot
  pca_scores <- as.data.frame(pca_result$x)  # Scores for individuals
  pca_scores$CONDITION <- metadata$CONDITION
  pca_scores$SAMPLE <- metadata$SAMPLE
  pca_scores$CovidBreadth <- metadata$CovidBreadth

  # Correlations with PCs (for annotations)
  cor_with_pc1 <- cor(pca_scores$CovidBreadth, pca_scores$PC1, method = "spearman")
  cor_with_pc2 <- cor(pca_scores$CovidBreadth, pca_scores$PC2, method = "spearman")

  # Find common limits for PC1 and PC2
  x_limits <- range(pca_scores$PC1, na.rm = TRUE)
  y_limits <- range(pca_scores$PC2, na.rm = TRUE)


  # PCA biplot
  pca_plot <- ggplot() +
    geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = CONDITION), size = 2.5) +
    geom_text_repel(data = pca_scores, aes(x = PC1, y = PC2, color = CONDITION, label = SAMPLE), nudge_x = 0.3, size = 3, segment.color = NA) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    # Plot top feature loadings as arrows
    geom_segment(data = top_loadings, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10), 
                arrow = arrow(length = unit(0.2, "cm")), color = "darkgray", alpha = 0.6) +
    # Annotate the features
    geom_text(data = top_loadings, aes(x = PC1 * 10, y = PC2 * 10, label = variable), size = 2.5, hjust = 1.2) +
    # Customize axis labels
    labs(x = paste0("PC1 (", round(explained_var[1], 2), "%)"), 
        y = paste0("PC2 (", round(explained_var[2], 2), "%)")) +
    # Customize colors
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(text = element_text(size = 8),
        legend.position = "none") +  # Remove legend here
    coord_cartesian(xlim = x_limits, ylim = y_limits)

  # Scatter plot of the variable against PC1
  scatter_pc1 <- ggplot(pca_scores, aes(x = PC1, y = CovidBreadth, color = CONDITION)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(y = "CovidBreadth", x = '', subtitle = paste("Spearman correlation with PC1:", round(cor_with_pc1, 5))) +
    coord_cartesian(xlim = x_limits)

  # Scatter plot of the variable against PC2
  scatter_pc2 <- ggplot(pca_scores, aes(x = CovidBreadth, y = PC2, color = CONDITION)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(x = "CovidBreadth", y = '', subtitle = paste("Spearman correlation with PC2:", round(cor_with_pc2, 5))) +
    coord_cartesian(ylim = y_limits)

  # Blank plot
  blankPlot <- ggplot() + geom_blank() +
    theme_void()

  # Extract the legend from the PCA plot
  legend_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = CONDITION)) +
    geom_point() +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(legend.position = "right")  # Keep the legend here

  # Extract the legend as a separate object
  legend <- get_legend(legend_plot)


  combined_plot <- ggarrange(
    scatter_pc1, blankPlot, pca_plot, scatter_pc2, 
    ncol = 2, nrow = 2, widths = c(3, 1.4), heights = c(1.4, 3),
    legend.grob = legend,  # Add the legend as a grob
    legend = "right"       # Position the legend to the right
  )


  # Save the plot
  ggsave(paste0(plotsDir, glue("pca_biplot_{time}_12.png")), plot = combined_plot, width = 15, height = 10, dpi = 600)



 
  loading_df <- data.frame(
    variable = rownames(loadings),
    PC3 = loadings[, 3],
    PC4 = loadings[, 4]
  )

  # Define quadrants based on PC1 and PC2 signs
  loading_df$quadrant <- with(loading_df, ifelse(PC3 > 0 & PC4 > 0, "Q1",
                                    ifelse(PC3 < 0 & PC4 > 0, "Q2",
                                    ifelse(PC3 < 0 & PC4 < 0, "Q3", "Q4"))))

  # Select the top 5 variables for each quadrant
  top_features <- loading_df %>%
    group_by(quadrant) %>%
    top_n(5, wt = abs(PC3) + abs(PC4))  # Sum of absolute loadings for PC3 and PC4
  
  # Filter loadings to keep only top features
  top_loadings <- loadings %>%
    filter(variable %in% top_features$variable)


  # Correlations with PCs (for annotations)
  cor_with_pc1 <- cor(pca_scores$CovidBreadth, pca_scores$PC3, method = "spearman", use = "complete.obs")
  cor_with_pc2 <- cor(pca_scores$CovidBreadth, pca_scores$PC4, method = "spearman", use = "complete.obs")

  x_limits <- range(pca_scores$PC3, na.rm = TRUE)
  y_limits <- range(pca_scores$PC4, na.rm = TRUE)



  # Create a ggplot biplot
  pca_plot <- ggplot() +
    # Plot individuals
    geom_point(data = pca_scores, aes(x = PC3, y = PC4, color = CONDITION), size = 2.5) +
    geom_text_repel(data = pca_scores, aes(x = PC3, y = PC4, color = CONDITION, label = SAMPLE), nudge_x = 0.3, size = 3, segment.color = NA) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = .2) +
    # Plot top feature loadings as arrows
    geom_segment(data = top_loadings, aes(x = 0, y = 0, xend = PC3 * 10, yend = PC4 * 10), 
                arrow = arrow(length = unit(0.2, "cm")), color = "darkgray", alpha = 0.6) +
    # Annotate the features
    geom_text(data = top_loadings, aes(x = PC3 * 10, y = PC4 * 10, label = variable), size = 2.5, hjust = 1.2) +
    # Customize axis labels
    labs(x = paste0("PC3 (", round(explained_var[3], 2), "%)"), 
        y = paste0("PC4 (", round(explained_var[4], 2), "%)")) +
    # Customize colors
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(text = element_text(size = 8),
        legend.position = "none") +  # Remove legend here
    coord_cartesian(xlim = x_limits, ylim = y_limits)

    # Scatter plot of the variable against PC3
  scatter_pc3 <- ggplot(pca_scores, aes(x = PC3, y = CovidBreadth, color = CONDITION)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(y = "CovidBreadth", x = '', subtitle = paste("Spearman correlation with PC3:", round(cor_with_pc1, 5))) +
    coord_cartesian(xlim = x_limits)


  # Scatter plot of the variable against PC2
  scatter_pc4 <- ggplot(pca_scores, aes(x = CovidBreadth, y = PC4, color = CONDITION)) +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "none",
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank()) +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    labs(x = "CovidBreadth", y = '', subtitle = paste("Spearman correlation with PC4:", round(cor_with_pc2, 5))) +
    coord_cartesian(ylim = y_limits)

  # Blank plot
  blankPlot <- ggplot() + geom_blank() +
    theme_void()

  # Extract the legend from the PCA plot
  legend_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = CONDITION)) +
    geom_point() +
    scale_color_manual(values = c("Healthy" = "blue", "Lymphomas" = "red")) +
    theme_bw() +
    theme(legend.position = "right")  # Keep the legend here

  # Extract the legend as a separate object
  legend <- get_legend(legend_plot)


  combined_plot <- ggarrange(
    scatter_pc3, blankPlot, pca_plot, scatter_pc4, 
    ncol = 2, nrow = 2, widths = c(3, 1.4), heights = c(1.4, 3),
    legend.grob = legend,  # Add the legend as a grob
    legend = "right"       # Position the legend to the right
  )

  # Save the plot
  ggsave(paste0(plotsDir, glue("pca_biplot_{time}_34.png")), plot = combined_plot, width = 15, height = 10, dpi = 600)




  library(variancePartition)
  require(tidyr)

  pca_scores <- as.data.frame(pca_result$x)
  rownames(pca_scores) <- rownames(metadata)

  # Step 2: Print the dimensions of variance_input and metadata for verification
  print("Dimensions of variance_input:")
  print(dim(data))
  
  print("Dimensions of metadata:")
  print(dim(metadata))
  
  # Step 3: Check and print summary of metadata
  print("Metadata summary:")
  print(summary(metadata))


  form <- ~ (1 | CONDITION) + (1 | HLA_A) + (1 | HLA_B) + (1 | HLA_C)



 # Step 5: Fit the model, check for errors
  print("Fitting the model...")
  varPart <- tryCatch({
    fitExtractVarPartModel(t(pca_scores), form, metadata)
  }, error = function(e) {
    print(paste("Error in model fitting:", e))
    return(NULL)
  })

  if (is.null(varPart)) {
    stop("Model fitting failed. Check the input data and model formula.")
  }

  # Step 6: Print and check variance partitioning results
  print("Variance partitioning results:")
  print(summary(varPart))
  
  # Step 7: Sort the variance results
  vp <- sortCols(varPart)
  print("Sorted variance partitioning results:")
  print(vp)

   # Step 8: Plot percentage variance explained
  png(paste0(plotsDir, glue('variance_explained_{time}.png')), width = 800, height = 500, res = 100)
  p <- plotPercentBars(vp[1:5, ])
  print(p)
  dev.off()

  # Step 9: Plot global variance explained
  png(paste0(plotsDir, glue('global_variance_explained_{time}.png')), width = 800, height = 800, res = 100)
  p <- plotVarPart(vp)
  print(p)
  dev.off()


}



test_cluster_vgenes <- function(trbv_gene_usage, metadata_test, cluster, time, cluster_n, plotsDir){
  # Required libraries
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(reshape2)

  # Transform list to array
  cluster <- unlist(cluster)

  # Merge the TRBV gene usage data with metadata for cluster1 and not_cluster1
  cluster_data <- trbv_gene_usage %>%
    filter(TIMEPOINTS == time) %>%
    inner_join(metadata_test %>% select(sample_id, SAMPLE), by = "sample_id") %>%
    filter(SAMPLE %in% cluster) %>%
    select(-SAMPLE, -sample_id, -CONDITION, -TIMEPOINTS)

  not_cluster_data <- trbv_gene_usage %>%
    filter(TIMEPOINTS == time) %>%
    inner_join(metadata_test %>% select(sample_id, SAMPLE), by = "sample_id") %>%
    filter(!SAMPLE %in% cluster) %>%
    select(-SAMPLE, -sample_id, -CONDITION, -TIMEPOINTS)

  cluster_data$Group <- "Cluster"
  not_cluster_data$Group <- "NotCluster"

  # Combine the two data frames
  plot_data <- bind_rows(
    mutate(cluster_data, Group = "Cluster"),
    mutate(not_cluster_data, Group = "NotCluster")
  )

  plot_data <- melt(plot_data, id.vars = "Group", variable.name = "VGene", value.name = "Value")


  # Initialize a list to store p-values
  p_values <- numeric(length(unique(plot_data$VGene)))

  # Perform Wilcoxon test for each V gene
  for (v_gene in unique(plot_data$VGene)) {
    cluster_values <- plot_data$Value[plot_data$VGene == v_gene & plot_data$Group == "Cluster"]
    not_cluster_values <- plot_data$Value[plot_data$VGene == v_gene & plot_data$Group == "NotCluster"]
    
    # Perform the Wilcoxon test
    test_result <- wilcox.test(cluster_values, not_cluster_values, exact = FALSE)
    p_values[which(unique(plot_data$VGene) == v_gene)] <- test_result$p.value
  }

  # Create a data frame for p-values
  p_value_df <- data.frame(
    VGene = unique(plot_data$VGene),
    PValue = p_values
  )
  # Adjust p-values for multiple testing (optional)
  p_value_df$PValue[is.nan(p_value_df$PValue)] <- 1
  p_value_df$AdjustedPValue <- p.adjust(p_value_df$PValue, method = "fdr")

  p_value_df$Significance <- cut(p_value_df$PValue,
                               breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 1),
                               labels = c("****", "***", "**", "*", ""))

  
  # Plot with larger text and annotations
  p <- ggplot(plot_data, aes(x = VGene, y = Value, color = Group)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 2) +
    theme_bw() +
    labs(y = "V gene clonal fraction", x = "V Gene") +
    theme(
      legend.position = "bottom",
      text = element_text( size = 14),  # General text size
      axis.title = element_text(size = 16),  # Axis title size
      axis.text = element_text(size = 14),   # Axis label size
      legend.text = element_text(size = 12),  # Legend text size
      axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
    ) +
    coord_cartesian(ylim = c(min(plot_data$Value), max(plot_data$Value) + 0.000035)) +
    scale_color_manual(values = c("Cluster" = "black", "NotCluster" = "darkgray")) +
    geom_text(data = p_value_df, aes(x = VGene, y = max(plot_data$Value) + 0.0000001, label = Significance),
                   color = "black", size = 9, angle = 0, vjust = 0) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))



  ggsave(paste0(plotsDir, glue("Vgene_spec_cluster{cluster_n}.png")), plot = p, width = 20, height = 6, dpi = 300)

}


vgene_lmm <- function(df, plotsDir) {

  require(glue)
  require(ggplot2)
  require(lme4)
  require(lmerTest)

  df$CONDITION <- as.factor(df$CONDITION)
  df$SAMPLE <- as.factor(df$SAMPLE)
  df$v_call <- as.factor(df$v_call)
  df$concave <- as.factor(df$concave)

  ###################################################################### Linear Mixed Model with the modules
  model <- lmer(cloneFraction ~ (TIMEPOINTS * v_call * CONDITION) + (1 | SAMPLE), data = df)

  llm_pvalues <- summary(model)$coefficients[, "Pr(>|t|)"]
  llm_pvalues <- llm_pvalues[llm_pvalues < 0.05]


  p <- plot_summs(model, ci_level = .95, model.names = c(glue("Model LLM")), coefs = names(llm_pvalues))

  # annotate the plot with the p values taken from the models
  i = 1
  for (name in names(llm_pvalues)) {
    value <- llm_pvalues[name]
    plot_data <- ggplot_build(p)$data[[1]]
    p <- p + annotate("text", x = plot_data$xmax[i], y = plot_data$y[i],
                      label = pvalue_to_asterisk(summary(model)$coefficients[name,]['Pr(>|t|)']), size = 5)
    i = i + 1
  }
  ggsave(paste0(plotsDir, glue("models_summs.png")), p, width = 10, height = 10)


}