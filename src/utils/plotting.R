plot_tcrex_predictions = function(data, path){

    require(ggplot2)
    require(tidyverse)

    p <-ggplot(data, aes(x=TIMEPOINTS,y=cloneFractionCovid, fill=CONDITION)) +
    geom_boxplot(notch = FALSE, width = 0.5) +
    geom_point() +
    geom_line(aes(group=SAMPLE)) +
    theme_bw(base_size = 15) +
    ylim(0.1, 0.3) +
    labs(title = "Repertoire clonal fraction of Covid-specific TCRs",
         y = "Clonal fraction",
         x = "Timepoint") +
    theme(
      axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",family = "JetMono Brains"),
    )
    ggsave(paste0(path,"tcrexpreds_R.png"), scale = 1, width = 8, height = 8)
}



plot_tcrex_predictions2 = function(data, path){
  # Box plot facetted by "dose"
  require(ggpubr)
  p <- ggboxplot(data, x = "CONDITION", y = "cloneFractionCovid",
                 color = "cloneFractionCovid", #palette = "jco",
                 add = "jitter",
                 facet.by = "TIMEPOINTS", short.panel.labs = FALSE)
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5)
  ggsave(paste0(path,"tcrexpreds_R.png"))

}








### TEST
workingDir = "/Users/fabioaffaticati/Desktop/Work/lymphoma_covid/"
setwd(workingDir) 
data <- read.csv("plotteddata.csv", row.names=1)

plot_tcrex_predictions2(data, '/Users/fabioaffaticati/Desktop/Work/lymphoma_covid/results/plots/')
