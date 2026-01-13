wd <- 'path/to/GitHub/copy-number-annotation/'

gains <- read.table(wd, "data/gain/gain_table.tsv", header = T, sep = "\t")
gains <- gains[gains$model != "Output.regressor_NoCluster_no_cluster",]

gains_ampl <- gains[gains$scna == "Amplification model",]
gains_del <- gains[gains$scna == "Deletion model",]

occ <- c("Centromere_Length", "Chromosome_Length", "dist.to.closest.FGS", "distance.to.telomere", "distance.to.centromere", "Length_Counts.E25", "Length_Counts.E22", "Length_Counts.E21", "Length_Counts.E19", "Length_Counts.E17", "Length_Counts.E13")
gains_ampl_occ <- gains_ampl[gains_ampl$Feature %in% occ,]
gains_del_occ <- gains_del[gains_del$Feature %in% occ,]

sel <- c("dist.to.closest.TSG", "dist.to.closest.OG", "all.int.trans", "Ess.distance_pancancer", "mutations_norm", "genes.bin", "partners.trans")
gains_ampl_sel <- gains_ampl[gains_ampl$Feature %in% sel,]
gains_del_sel <- gains_del[gains_del$Feature %in% sel,]

model_order <- c("Small-scale", "Mid-length", "Arm-level", "Chromosome-level")

gains_ampl_occ$model <- factor(gains_ampl_occ$model, levels = model_order)
gains_del_occ$model  <- factor(gains_del_occ$model,  levels = model_order)

gains_ampl_sel$model <- factor(gains_ampl_sel$model, levels = model_order)
gains_del_sel$model  <- factor(gains_del_sel$model,  levels = model_order)

p1 <- ggplot(gains_ampl_sel, aes(x = model, y = Gain, color = Feature, group = Feature)) +
  geom_line(size = 1) +
  geom_point() +
  scale_x_discrete(drop = FALSE) +
  theme_minimal()

ggsave(paste0(wd, "data/plots/gains_ampl_sel.pdf"), p1, width = 7, height = 5)


p2 <- ggplot(gains_del_occ, aes(x = model, y = Gain, color = Feature, group = Feature)) +
  geom_line(size = 1) +
  geom_point() +
  scale_x_discrete(drop = FALSE) +
  theme_minimal()

ggsave(paste0(wd, "data/plots/gains_del_occ.pdf"), p2, width = 7, height = 5)


p3 <- ggplot(gains_del_sel, aes(x = model, y = Gain, color = Feature, group = Feature)) +
  geom_line(size = 1) +
  geom_point() +
  scale_x_discrete(drop = FALSE) +
  theme_minimal()

ggsave(paste0(wd, "data/plots/gains_del_sel.pdf"), p3, width = 7, height = 5)


p4 <- ggplot(gains_ampl_occ, aes(x = model, y = Gain, color = Feature, group = Feature)) +
  geom_line(size = 1) +
  geom_point() +
  scale_x_discrete(drop = FALSE) +
  theme_minimal()

ggsave(paste0(wd, "data/plots/gains_ampl_occ.pdf"), p4, width = 7, height = 5)


