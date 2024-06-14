plot_predictions = function(data, path){

    require(ggplot2)
    require(tidyverse)

    p <-ggplot(data, aes(x=TIMEPOINTS,y=cloneFractionCovid, fill=CONDITION)) +
    geom_boxplot(notch = FALSE, width = 0.5) +
    geom_point() +
    geom_line(aes(group=SAMPLE, colour = CONDITION)) +
    theme_bw(base_size = 15) +
    #ylim(0.1, 0.3) +
    labs(title = "Repertoire clonal fraction of Covid-specific TCRs",
         y = "Clonal fraction",
         x = "Timepoint") +
    theme(
      axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",family = "JetMono Brains"),
    )
    ggsave(paste0(path,"tcrexpreds_R.png"), scale = 1, width = 8, height = 8)
}



plot_scatteratio_breadth = function(data, path){

  require(ggplot2)
  require(tidyverse)
  
  p <-ggplot(data, aes(x=TIMEPOINTS,y=fraction_betas, fill=CONDITION)) +
    geom_boxplot(notch = FALSE, width = 0.5) +
    geom_point() +
    geom_line(aes(group=SAMPLE, colour = CONDITION)) +
    theme_bw(base_size = 15) +
    #ylim(0.1, 0.3) +
    labs(title = "",
         y = "Covid Specific fraction betas (breadth)",
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





### TEST
workingDir = "/Users/fabioaffaticati/Desktop/Work/lymphoma_covid/"
setwd(workingDir) 
data <- read.csv("plotteddata.csv", row.names=1)

data <- read.csv("heatmap.csv", row.names = 1)





