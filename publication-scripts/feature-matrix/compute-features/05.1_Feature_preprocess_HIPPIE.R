# Date: August 8, 2024 - Updated on October 3, 2024
# Adding PPI related features based on PPIs from HIPPIE to backbone files

rm(list = ls())

packages <- c(
  "ggplot2", "ggpubr"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

# 1. Calculate new features 

#Most of the non Tissue-specific features are gene based. Therefore, we need genes located on the corresponding bin.
#Genes per bins for focal CNA

backbone_path <- "./Data/All_levels_backbonetables.RData"
genes_per_bins_table <- "./Data/All_levels_genes_per_bins.RData"
ensembl_corum_path <- "./Data/Ensembl_CORUM_features_tissue_general.RData"
hippie_path <- "./Data/HIPPIE_PPIs_tissue_general.RData"
output_path <- "./Data/Backbone_tables_with_non_tissue_specific_features_CORUM_HIPPIE.RData"

load(genes_per_bins_table)

#All backbone tables for focal CNAs and aneuploidies
load(backbone_path)

#Features: Ensembl and CORUM related features (Just to use CORUM)
load(ensembl_corum_path)

gene.table <- gene.table[,c(1,2,17:21)]

#Features: PPIs from HIPPIE
load(hippie_path)
colnames(gene.table.hippie) <- c("Gene.name","PPIs-HIPPIE","n_PPIs-HIPPIE",
                                 "PPIs_063-HIPPIE","n_PPIs_063-HIPPIE",
                                 "PPIs_073-HIPPIE","n_PPIs_073-HIPPIE")
gene.table.hippie[,c(3,5,7)] <- sapply(gene.table.hippie[,c(3,5,7)], as.numeric)

#Merge two gene tables
gene.table <- merge(gene.table, gene.table.hippie, by = "Gene.name", all.x = T)

# Correlation check
ggplot(gene.table, aes(x=n_PPIs, y=`n_PPIs-HIPPIE`)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) +
  xlab("Interactome INSIDER") + ylab("HIPPIE") +
  theme_bw() + stat_cor(method = "pearson")

ggplot(gene.table, aes(x=n_PPIs, y=`n_PPIs_063-HIPPIE`)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) +
  xlab("Interactome INSIDER") + ylab("HIPPIE - Cut-off = .63") +
  theme_bw() + stat_cor(method = "pearson")

ggplot(gene.table, aes(x=n_PPIs, y=`n_PPIs_073-HIPPIE`)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE) +
  xlab("Interactome INSIDER") + ylab("HIPPIE - Cut-off = .73") +
  theme_bw() + stat_cor(method = "pearson")

levels <- names(chr_backbone_namesfixed)
backbone_tables_w_features <- list()

for(level in levels){
  backbone <- do.call(rbind,chr_backbone_namesfixed[[level]])
  if(!level %in% c("Arm","Chromosome")){backbone$bin <- paste(backbone$chr,backbone$bin,sep = "_")} #To match the bin names between backbone and genes_on_bins tables
  features.df <- c()
  for(bin in backbone$bin){
    bin.genes <- genes_on_bins[[level]][[bin]]
    #Some bins do not have genes, and those are not included in genes_on_bins data
    if(length(bin.genes) == 0){row <- c(bin,"",length(bin.genes),rep(NA,17))}
    else{
      df.genes <- gene.table[gene.table$Gene.name %in% bin.genes,] #Gene table only has protein coding genes!
      #Find complex partners in trans
      partners <- unique(unlist(lapply(df.genes$Partners, function(x) strsplit(x,"::")[[1]])))
      partners <- partners[complete.cases(partners)]
      partners.trans <- setdiff(partners,bin.genes)
      #Find interacting proteins in trans based on HIPPIE (no-cutoff)
      ppis <- unique(unlist(lapply(df.genes$`PPIs-HIPPIE`, function(x) strsplit(x,"::")[[1]])))
      ppis <- ppis[complete.cases(ppis)]
      ppis.trans <- setdiff(ppis,bin.genes)
      #Find interacting proteins in trans based on HIPPIE (.63)
      ppis.063 <- unique(unlist(lapply(df.genes$`PPIs_063-HIPPIE`, function(x) strsplit(x,"::")[[1]])))
      ppis.063 <- ppis.063[complete.cases(ppis.063)]
      ppis.063.trans <- setdiff(ppis.063,bin.genes)
      #Find interacting proteins in trans based on HIPPIE (.73)
      ppis.073 <- unique(unlist(lapply(df.genes$`PPIs_073-HIPPIE`, function(x) strsplit(x,"::")[[1]])))
      ppis.073 <- ppis.073[complete.cases(ppis.073)]
      ppis.073.trans <- setdiff(ppis.073,bin.genes)
      #Create bin row with features
      row <- c(bin,
               paste(bin.genes,collapse = "::"), #Genes on the bin
               length(bin.genes), #Number of genes on the bin
               #Number of protein complex subunit on the bin
               length(which(df.genes$ComplexSubunit == "True")),
               #Partners and number of partners
               paste(partners, collapse = "::"),length(partners), 
               #Partners in trans and number of partners in trans
               paste(partners.trans, collapse = "::"),length(partners.trans),
               #Interacting proteins and number of them (no cut-off)
               paste(ppis, collapse = "::"),length(ppis),
               #Interacting proteins in trans and number of them (no cut-off)
               paste(ppis.trans, collapse = "::"), length(ppis.trans),
               #Interacting proteins and number of them (.63)
               paste(ppis.063, collapse = "::"),length(ppis.063),
               #Interacting proteins in trans and number of them (.63)
               paste(ppis.063.trans, collapse = "::"), length(ppis.063.trans),
               #Interacting proteins and number of them (.73)
               paste(ppis.073, collapse = "::"),length(ppis.073),
               #Interacting proteins in trans and number of them (.73)
               paste(ppis.073.trans, collapse = "::"), length(ppis.073.trans)
               )
      }
    features.df <- rbind(features.df,row)}
  features.df <- as.data.frame(features.df)
  colnames(features.df) <- c("bin",
                             "genes",
                             "n_genes",
                             "total.complex.subunit",
                             "Partners","total_n_partners",
                             "Partners.trans","total_n_partners.trans",
                             "PPIs-HIPPIE","total_n_PPIs-HIPPIE",
                             "PPIs.trans-HIPPIE","total_n_PPIs.trans-HIPPIE",
                             "PPIs_063-HIPPIE","total_n_PPIs_063-HIPPIE",
                             "PPIs_063.trans-HIPPIE","total_n_PPIs_063.trans-HIPPIE",
                             "PPIs_073-HIPPIE","total_n_PPIs_073-HIPPIE",
                             "PPIs_073.trans-HIPPIE","total_n_PPIs_073.trans-HIPPIE")
  features.df[,c(3,4,6,8,10,12,14,16,18,20)] <- sapply(features.df[,c(3,4,6,8,10,12,14,16,18,20)], as.numeric)
  backbone <- merge(backbone,features.df, by="bin")
  backbone_tables_w_features[[level]] <- backbone}

save(backbone_tables_w_features, file = output_path)

# After this, use the output to calculate mRNA and protein mean/medians in 03.1_Feature_preprocess_mrna_protein_levels.R

