# Custom breaks and labels function
custom_breaks <- function(x) {
        10^seq(floor(log10(min(x))), ceiling(log10(max(x))), by = 1)
}




plot_clonal_depth = function(data, path){

    require(ggplot2)
    require(tidyverse)
    require(scales)
    
    # Define colors for conditions
    colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")

    p <-ggplot(data, aes(x=TIMEPOINTS,y=cloneFraction, fill=CONDITION)) +
    geom_boxplot(notch = FALSE, width = 0.5, outlier.shape = NA) +
    geom_point() +
    scale_y_log10(breaks = custom_breaks, labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Logarithmic scale with custom labels
    annotation_logticks(sides = "l") +
    scale_fill_manual(values = colors) +  # Set fill colors manualsly
    scale_colour_manual(values = colors) +  # Set line colors manually
    geom_line(aes(group=SAMPLE, colour = CONDITION)) +
    theme_bw(base_size = 10) +
    #ylim(0.1, 0.3) +
    labs(title = "Repertoire clonal fraction of Covid-specific TCRs",
         y = "Clonal fraction",
         x = "Timepoint") +
    theme(
      axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",family = "JetMono Brains"),
    )
    ggsave(paste0(path,"clonal_preds_R.png"))
}



plot_scatteratio_breadth = function(data, path){

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


plot_heatmap_predictions = function(data, path){

  require(ggplot2)
  require(pheatmap)
  require(RColorBrewer)
  require(reshape2)
  
  color_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
  
  # Reshape data to wide format
  data_wide <- dcast(data, epitope ~ sample_id, value.var = "cloneFraction", fun.aggregate = sum)
  
  # Set rownames to epitope and remove the epitope column
  rownames(data_wide) <- data_wide$epitope
  data_wide$epitope <- NULL
  
  # Create annotation data frames for columns
  sample_annotations <- data[!duplicated(data$sample_id), c("sample_id", "TIMEPOINTS", "CONDITION")]
  rownames(sample_annotations) <- sample_annotations$sample_id
  sample_annotations$sample_id <- NULL
  
  # Sort the annotations by condition
  sorted_annotations <- sample_annotations[order(sample_annotations$TIMEPOINTS), ]
  
  # Reorder the columns of data_wide to match the sorted annotations
  data_wide <- data_wide[, rownames(sorted_annotations)]
  
  
  # Choose a color palette
  color_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
  
  # Define annotation colors
  anno_colors <- list(
    timepoints = c('baseline' = "#E41A1C", 'V1' = "#377EB8", 'V3' = "#4DAF4A"),
    condition = c('Healthy' = "blue", 'Lymphomas' = "red")
  )
  
  # Create the heatmap
  p <- pheatmap(data_wide,
           color = color_palette,
           cluster_rows = F,
           cluster_cols = F,
           #scale = "column",
           #scale = "none",
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize = 10,
           fontsize_row = 10,
           fontsize_col = 10,
           main = "Epitope immunodominance",
           annotation_col = sample_annotations,
           annotation_colors = anno_colors,
           border_color = NA)
  
  ggsave(paste0(path,"heatmap.png"), p)
  
}


plot_betweenstats_breadth = function(df, path){

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



plot_betweenstats_clone_fraction = function(df, path){
  
  require(ggstatsplot)
  require(ggplot2)
  require(ggrepel)
  require(ggsignif)
  require(tidyverse)
  require(rstatix)
  require(ggpubr)
  
  colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")

  df <- df %>% convert_as_factor(SAMPLE, TIMEPOINTS)
  res.fried <- df %>% friedman_test(cloneFraction ~ TIMEPOINTS | SAMPLE)
  pwc <- df %>% wilcox_test(cloneFraction ~ TIMEPOINTS, paired = TRUE, p.adjust.method = "BH")
  pwc <- pwc %>% add_xy_position(x = "TIMEPOINTS")
  t <- ggboxplot(df, x = "TIMEPOINTS", y = "cloneFraction", add = "point") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(
      subtitle = get_test_label(res.fried,  detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  ggsave(paste0(path,"test.png"), t, width=8)

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
                                   boxplot.args = list(width = 0.3, alpha = 0.2, color = "gray"),
                                   points.args = list(size = 30, alpha = 1),
                                   violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                                   plot.type = "boxplot", type = "nonparametric",
                                   geom_signif_args = list(textsize = 30),
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
    geom_label_repel(stat = "summary", fun = median, size = 3,
                     aes(label = paste0("hat(mu)*scriptstyle(median)==", 
                                        round(after_stat(y), 5))),
                     parse = TRUE, position = position_dodge(0.9)) +
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
         x = "Timepoint")
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



plot_gene_usage = function(df, path){

  require(ggstatsplot)
  require(ggplot2)
  require(ggrepel)
  require(ggsignif)
  require(tidyverse)
  require(rstatix)
  require(ggpubr)
  
  colors <- c("Healthy" = "#1f77b4", "Lymphomas" = "#d62728")


  p1 <- ggstatsplot:: grouped_ggbetweenstats(data = df,
                                   x = CONDITION,
                                   y = cloneFraction,
                                   ylab = 'Clone Fraction',
                                   grouping.var = v_call,
                                   #title = 'Global Comparison of Breadth Across Timepoints',
                                   #pairwise.display = 'all',
                                   pairwise.display = 's',
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


  options(ragg.max_dim = 1000000)

  ggsave(paste0(path,"test_gene_usage.png"), p1, width=1000, limitsize = FALSE)

}











### TEST
workingDir = "/Users/fabioaffaticati/Desktop/Work/lymphoma_covid/"
setwd(workingDir) 

library(lme4)
library(lmerTest)
library(ggplot2)
df <- read.csv("scatteratio_data.csv", row.names = 1)


df$SAMPLE <- as.factor(df$SAMPLE)
df$TIMEPOINTS <- as.factor(df$TIMEPOINTS)
df$CONDITION <- as.factor(df$CONDITION)

model <- lmer(fraction_sequences ~ CONDITION * TIMEPOINTS + (1 | SAMPLE), data = df)
summary(model)


df <- read.csv("clonal_fractions_scatteratio_data.csv", row.names = 1)

df$SAMPLE <- as.factor(df$SAMPLE)
df$TIMEPOINTS <- as.factor(df$TIMEPOINTS)
df$CONDITION <- as.factor(df$CONDITION)

model <- lmer(cloneFraction ~ CONDITION * TIMEPOINTS + (1 | SAMPLE), data = df)
summary(model)
