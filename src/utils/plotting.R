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
  require(ggplot2)


  df$SAMPLE <- as.factor(df$SAMPLE)
  df$TIMEPOINTS <- as.factor(df$TIMEPOINTS)
  df$CONDITION <- as.factor(df$CONDITION)

  # if fraction_sequences is in the columns of the dataframe
  if("fraction_sequences" %in% colnames(df)){
    model <- lmer(fraction_sequences ~ CONDITION * TIMEPOINTS + (1 | SAMPLE), data = df)
    saveRDS(model, file = paste0(model_path, "lmm_model_breadth.rds"))
  }
  # if cloneFraction is in the columns of the dataframe
  else if("cloneFraction" %in% colnames(df)){
    model <- lmer(cloneFraction ~ CONDITION * TIMEPOINTS + (1 | SAMPLE), data = df)
    saveRDS(model, file = paste0(model_path, "lmm_model_depth.rds"))
  }
  else{
    print("No suitable column found in the dataframe")
  }


}

export_lmm_results <- function(model_path){

  require(jtools)

  # read RDS files
  model_breadth <- readRDS(paste0(model_path, 'lmm_model_breadth.rds'))
  model_depth <- readRDS(paste0(model_path, 'lmm_model_depth.rds'))

  #export_summs(model_breadth,  model_depth, to.file = 'pdf', file.name = paste0(model_path, "lmm_model_reports.pdf"),
  #            model.names = c("LMM_Breadth", "LMM_Depth"), digits = 6, to.numeric = TRUE, error_format = "[{conf.low}, {conf.high}]")

  
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