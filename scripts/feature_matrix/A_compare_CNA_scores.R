# rm(list=ls())
gc(full=T)

ml_merged <- readRDS("~/mountHD/ml_models/data/ml_merged.rds")

tumor_type <- 'BRCA'

pdf(file = 'compare_CNA_scores.pdf', width = 15)
for(i in names(ml_merged)){
  toPlot <- ml_merged[[i]]
  toPlot <- toPlot[toPlot$type == tumor_type,]
  
  toPlot[,c(1:10,12:24)] <- apply(toPlot[,c(1:10,12:24)], 2, as.numeric)
  
  print(ggarrange(ggplot(toPlot, aes(x = del_score, y = del_score_sg)) +
                    geom_point() +
                    ggtitle(i),
                  ggplot(toPlot, aes(x = ampl_score, y = ampl_score_sg)) +
                    geom_point() +
                    ggtitle(i)))
}
dev.off()