# useful notes/data for shiny app

setwd('Desktop/Data shiny/') # set this folder as directory

# SHAP values of the models are
load('Avg_shap_values_InteractomeINSIDER.RData')
shap.df <- models.shap.df$`Mid-length::Amplification model`

# remove features that are always 0 across columns
shap.df <- shap.df[,apply(shap.df, 2, function(x){all(x!=0)})]

shap.df$type <- do.call(rbind, str_split(shap.df$labels, pattern = '-'))[,2]
shap.df$chr <- do.call(rbind, str_split(shap.df$labels, '_'))[,1]
shap.df$bin <- do.call(rbind, str_split(do.call(rbind, str_split(shap.df$labels, '_'))[,2], '-'))[,1]

# remove all non-segment features and useless columns
shap.df <- shap.df %>% select(-c("Chromosome_Length", "Centromere_Length", "Centromere_Type", 'BIAS'))
shap.df <- shap.df %>% select("labels", "type", "chr", "bin",
                              "dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                              "distance.to.centromere", "distance.to.telomere", "mutations_norm",
                              "genes.bin", "Length_Counts.E17", "Length_Counts.E19", "Ess.distance_pancancer",
                              "Length_Counts.E25", "Length_Counts.E1")

View(shap.df) # this is the table you should use for plotting sum(abs(SHAP)) for each point

# example of the plot for tt <- 'BRCA' and chromosome 1 and 2
shap.abs.sum <- data.frame(apply(shap.df %>% filter(type %in% 'BRCA' & chr %in% c(1,2)) %>% select(-c('labels','type','chr','bin')), 2, 
                                 function(x){sum(abs(x))}))
colnames(shap.abs.sum) <- 'value'
shap.abs.sum$feature <- rownames(shap.abs.sum)

ggplot(shap.abs.sum, aes(x = reorder(feature, value), y = value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "SHAP Value Contribution per Feature",
       x = "Feature",
       y = "Sum of Absolute SHAP Values") +
  theme_minimal()

# Quindi questo plot ideamente è quello che vorrei sotto la schermata
# uno per modello (uno ampl e uno del)



# per il plot grande, assumento che questa sia la tabella finale 
toplot.plot <- readRDS('tabella_finale?.rds')

# Compute chromosome ranges for background shading
chr_bounds <- toplot.plot %>%
  group_by(chr) %>%
  summarize(start = min(pos), end = max(pos)) %>%
  mutate(fill = rep(c("lightgrey", "darkgrey"), length.out = n()))

base_plot <- ggplot() +
  geom_rect(data = chr_bounds,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.3) +
  scale_fill_manual(values = c("#eeeeee", "#cccccc")) +
  geom_line(data = toplot.plot, aes(x = pos, y = ampl_score), color = "red") +
  geom_line(data = toplot.plot, aes(x = pos, y = -del_score), color = "blue") +
  geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
  labs(x = "Genomic Position", y = "Amplification Score (Mid-length)") +
  theme_classic() +
  theme(legend.position = 'none')

# Add cluster lines above plot
cluster_ticks <- toplot.plot %>%
  mutate(cluster_ymin = max(toplot.plot$ampl_score) * 1.07 + clusters15 * 0.015,
         cluster_ymax = cluster_ymin + 0.01)

base_plot +
  geom_segment(data = cluster_ticks,
               aes(x = pos, xend = pos, y = cluster_ymin, yend = cluster_ymax),
               color = "black", size = 0.3)
