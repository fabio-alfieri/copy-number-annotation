# Pre-processing of feature data
# Ensembl Biomart and CORUM features
# Tissue-general

rm(list = ls())

packages <- c(
  "stringr", "tidyr"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

wd <- 'path/to/GitHub/copy-number-annotation/'
setwd(wd)

# Data from Ensembl Biomart (downloaded on 14 February 2024 - Genome assembly GRCh37.p13)
# List of genes with the largest transcript (autosomal protein coding genes)

ensembl_biomart_path <- "feature_matrix_data/hg19_14022024.txt"
ccds_path <- "feature_matrix_data/CCDS.current_hg19.txt"
subunits_path <- ".feature_matrix_data/HumanPCgenes_Subunits_allchr.RData"
PPIs_intINSIDER <- "feature_matrix_data/H_sapiens_interfacesHQ_genename.RData"
hippie_tissue_general_path <- "feature_matrix_data/HIPPIE_PPIs_tissue_general.RData"
ensembl_corum_path <- "feature_matrix_data/Ensembl_CORUM_features_tissue_general.RData" 

hg19 <- read.delim(ensembl_biomart_path) # 215404
genes <- unique(hg19$Gene.name)
hg19 <- hg19[hg19$Gene.type == "protein_coding",] # 160238
hg19 <- hg19[hg19$Chromosome.scaffold.name %in% as.character(seq(1,22,1)),] # 141152 
hg19 <- hg19[order(hg19$Transcript.length..including.UTRs.and.CDS.,decreasing = T),]
hg19 <- hg19[!duplicated(hg19$Gene.name),] #19349 - The unique gene set with the longest transcript
hg19$Gene.length <- hg19$Gene.end..bp. - hg19$Gene.start..bp.
hg19 <- hg19[,c(1,9,3,8,10,11,14,15)]

#Data from CCDC - Keep the one with the longest CCDS
ccds <- read.delim(ccds_path) #27809 genes
ccds <- ccds[ccds$chromosome %in% as.character(seq(1,22,1)),] #26473 genes
ccds$CCDS.ID <- gsub("\\..*","",ccds$ccds_id)
ccds <- ccds[ccds$ccds_status == "Public",]
ccds[,c(8,9)] <- sapply(ccds[,c(8,9)],as.numeric)
#Keep the gene with the longest coding region
ccds$cds_length <- ccds$cds_to - ccds$cds_from
ccds <- ccds[order(ccds$cds_length, decreasing = T),]
ccds <- ccds[!duplicated(ccds$gene),] # 17559 - the longest cdc
ccds <- ccds[,c(1:5,8,9,12,13)]
#Merging two tables
gene.table <- merge(hg19,ccds,by.x="Gene.name",by.y = "gene",all = TRUE) #19621 genes (union)

rm(list = setdiff(ls(),"gene.table"))

#Adding if a gene is subunit or not, and if so, the number of partners & partners
#CORUM complexes 
load(subunits_path) #partner info for subunits
rm(list = setdiff(ls(),c("human_proteincomplex","gene.table")))
extended.complexes <- separate_rows(human_proteincomplex,subunits.Gene.name.,sep = ";")
extended.complexes <- extended.complexes[extended.complexes$subunits.Gene.name. != "",]
extended.complexes$subunits.Gene.name. <- str_trim(extended.complexes$subunits.Gene.name.) #remove white spaces from start/end of proteins
corum.proteins <- unique(extended.complexes$subunits.Gene.name.) #3662 proteins
#Add subunit info to the gene.table
gene.table$ComplexSubunit <- ifelse(gene.table$Gene.name %in% corum.proteins, "True","False")
#Add partners other than gene itself information
for(i in 1:nrow(gene.table)){
  gene <- gene.table$Gene.name[i]
  if(gene.table$ComplexSubunit[i] == "False"){
    gene.table$Partners[i] <- NA
    gene.table$n_partner[i] <- NA}
  else{
    partners <- extended.complexes[extended.complexes$ComplexName %in% 
                                     extended.complexes[extended.complexes$subunits.Gene.name. == gene,]$ComplexName,]$subunits.Gene.name.
    partners <- unique(setdiff(partners,gene))
    gene.table$Partners[i] <- paste(partners,collapse = "::")
    gene.table$n_partner[i] <- length(partners)}}

rm(list = setdiff(ls(),"gene.table"))
gc()

#Adding PPIs interactions from Interactome INSIDER
#Curated PPIs from Interactome INSIDER
load(PPIs_intINSIDER)
HS_interfaces_genenames$Homodimer <- ifelse(HS_interfaces_genenames$P1 == HS_interfaces_genenames$P2, "Yes","No")
HS_interfaces_genenames <- HS_interfaces_genenames[HS_interfaces_genenames$Homodimer != "Yes",] #114287
interactome.proteins <- union(HS_interfaces_genenames$P1, HS_interfaces_genenames$P2)
for(i in 1:nrow(gene.table)){
  gene <- gene.table$Gene.name[i]
  if(gene %in% interactome.proteins){
    bi <- HS_interfaces_genenames[HS_interfaces_genenames$P1 == gene |
                                    HS_interfaces_genenames$P2 == gene,]
    bi <- unique(setdiff(union(bi$P1,bi$P2),gene))
    gene.table$PPIs_IntINSIDER[i] <- paste(bi,collapse = "::")
    gene.table$n_PPIs_IntINSIDER[i] <- length(bi)}
  else{
    gene.table$PPIs_IntINSIDER[i] <- NA
    gene.table$n_PPIs_IntINSIDER[i] <- NA}}

rm(list = setdiff(ls(),"gene.table"))
gc()

# PPIs from HIPPIE
# # url for the data
# #url <- "https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt" # v2.3 - 831932 rows
# url <- "https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-current.mitab.txt" #v2.3 - 783182 rows
# download.file(url, destfile = "./Downloads/HIPPIE_PPIs.txt")
# ppis <- read.delim("./Downloads/HIPPIE_PPIs.txt", header = T)
# 
# all.genes <- union(ppis$Gene.Name.Interactor.A, ppis$Gene.Name.Interactor.B)
# 
# gene.table.hippie <- c()
# for(gene in all.genes){
#   bi <- ppis[ppis$Gene.Name.Interactor.A == gene | ppis$Gene.Name.Interactor.B == gene,]
#   # Remove homodimers
#   bi <- bi[bi$Gene.Name.Interactor.A != bi$Gene.Name.Interactor.B,]
#   bi <- bi[complete.cases(bi),]
#   # Confidence interval cutoff
#   bi_063 <- bi[bi$Confidence.Value >= 0.63,]
#   bi_073 <- bi[bi$Confidence.Value >= 0.73,]
#   # PPIs
#   bi <- unique(setdiff(union(bi$Gene.Name.Interactor.A,bi$Gene.Name.Interactor.B),gene))
#   PPIs <- paste(bi,collapse = "::") # without confidence cutoff
#   n_PPIs <- length(bi)
#   # PPIs with confidence interval cutoff
#   # Cutoff 0.63
#   PPIs_063 <- unique(setdiff(union(bi_063$Gene.Name.Interactor.A,bi_063$Gene.Name.Interactor.B),gene))
#   n_PPIs_063 <- length(PPIs_063)
#   PPIs_063 <- paste(PPIs_063,collapse = "::")
#   # Cutoff 0.73
#   PPIs_073 <- unique(setdiff(union(bi_073$Gene.Name.Interactor.A,bi_073$Gene.Name.Interactor.B),gene))
#   n_PPIs_073 <- length(PPIs_073)
#   PPIs_073 <- paste(PPIs_073,collapse = "::")
#   
#   gene.table.hippie <- rbind(gene.table.hippie, 
#                              c(gene,PPIs,n_PPIs,PPIs_063,n_PPIs_063,PPIs_073,n_PPIs_073))
# }
# gene.table.hippie <- as.data.frame(gene.table.hippie)
# colnames(gene.table.hippie) <- c("Gene","PPIs","n_PPIs","PPIs_063","n_PPIs_063","PPIs_073","n_PPIs_073")
# 
# save(gene.table.hippie, file = "./Data/HIPPIE_PPIs_tissue_general.RData") # Updated on October 2, 2024

load(hippie_tissue_general_path)
colnames(gene.table.hippie)[2:7] <- paste(colnames(gene.table.hippie)[2:7],"HIPPIE",sep = "_")
gene.table <- merge(gene.table, gene.table.hippie, by.x = "Gene.name", by.y = "Gene", all.x = TRUE)

save(gene.table, file = ensembl_corum_path)
