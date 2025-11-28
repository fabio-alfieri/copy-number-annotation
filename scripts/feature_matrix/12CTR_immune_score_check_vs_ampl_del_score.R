tt <- 'LUSC'
i <- '1Mbp'
hyper <- F
source('/mnt/fabiogokce/fabiohd/ml_models/code_ML/0_load_data.R', local = T)

# SECOND ----
x <- as.data.frame(cor(ml_tables_df[,c(features_background,
                                       features_selection,
                                       features_immune[grep('kruskall', features_immune)],
                                       features_immune[grep('wil.ampl', features_immune)],
                                       'ampl_score','del_score')],
                       use="complete.obs"))

y <- x$ampl_score
names(y) <- colnames(x)
y <- as.data.frame(y)
y$cell_type <- colnames(x)

y %>% filter(cell_type != 'ampl_score' & cell_type != 'del_score') %>%
  ggplot(aes(x = cell_type, y = y)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0.3, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = -0.3, col = 'red', linetype = 'dashed')

y <- x$del_score
names(y) <- colnames(x)
y <- as.data.frame(y)
y$cell_type <- colnames(x)

y %>% filter(cell_type != 'ampl_score' & cell_type != 'del_score') %>%
  ggplot(aes(x = cell_type, y = y)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0.3, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = -0.3, col = 'red', linetype = 'dashed')


heatmap(cor(ml_tables_df[,c(features_background,features_selection,
                            features_immune[grep('kruskall', features_immune)],
                            features_immune[grep('wil.ampl', features_immune)],
                            'ampl_score')], use="complete.obs"), na.rm = T)


# FIRST ----
x <- as.data.frame(cor(ml_tables_df[,c(features_background,features_selection,
                                       features_immune[grep('kruskall', features_immune)],
                                       # features_immune[grep('wil.ampl', features_immune)],
                                       'ampl_score','del_score')],
                       use="complete.obs"))

y <- x$ampl_score
names(y) <- colnames(x)
y <- as.data.frame(y)
y$cell_type <- colnames(x)

y %>% filter(cell_type != 'ampl_score' & cell_type != 'del_score') %>%
  ggplot(aes(x = cell_type, y = y)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0.3, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = -0.3, col = 'red', linetype = 'dashed')

y <- x$del_score
names(y) <- colnames(x)
y <- as.data.frame(y)
y$cell_type <- colnames(x)

y %>% filter(cell_type != 'ampl_score' & cell_type != 'del_score') %>%
  ggplot(aes(x = cell_type, y = y)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 0.3, col = 'red', linetype = 'dashed') +
  geom_hline(yintercept = -0.3, col = 'red', linetype = 'dashed')

heatmap(cor(ml_tables_df[,c(features_background,features_selection,
                            features_immune[grep('kruskall', features_immune)],
                            # features_immune[grep('wil.ampl', features_immune)],
                            'ampl_score')], use="complete.obs"), na.rm = T)



# check ML correlations ----
df <- read.table(file = '/mnt/fabiogokce/fabiohd/US_intern/comparison_feature_sets_immune.txt',
                 header = T)

ggpubr::ggarrange(df %>% filter(model == 'ampl') %>%
  ggplot(aes(y = cor, x = tumor_type, fill = feature_set)) +
  geom_bar(position = 'dodge', stat = 'identity'),
  df %>% filter(model == 'del') %>%
  ggplot(aes(y = cor, x = tumor_type, fill = feature_set)) +
  geom_bar(position = 'dodge', stat = 'identity'))
